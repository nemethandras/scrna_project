# scripts/demux.py
"""
Demultiplex single-cell RNA-seq data from a CellRanger BAM by matching
each cell barcode's genotype profile against known cell line genotypes
stored in the database.

Requires cellSNP-lite output (run via the cellsnp_pileup Snakemake rule).

Usage:
    python scripts/demux.py \
        --cellsnp-dir results/demux/my_run/cellsnp \
        --db          results/variants.db \
        --run-ids     HeLa_run1 HEK293_run1 K562_run1 \
        --output      results/demux/my_run/assignments.tsv \
        --demux-run-id my_run

Each cell is assigned to the cell line with the highest binomial
log-likelihood, or flagged as doublet/no_match.
"""

import sys
import sqlite3
import argparse
import numpy as np
import scipy.io
import scipy.sparse
from collections import Counter
from pathlib import Path
from scipy.special import gammaln
from cyvcf2 import VCF


# Sequencing error rate — prevents log(0) and models base-call errors.
# 0/0 is scored as p=ERR, 0/1 as p=0.5, 1/1 as p=1-ERR.
_ERR = 0.02

GT_DOSAGE = {
    "0/0": _ERR,       "0|0": _ERR,
    "0/1": 0.5,        "1/0": 0.5,
    "0|1": 0.5,        "1|0": 0.5,
    "1/1": 1.0 - _ERR, "1|1": 1.0 - _ERR,
}


# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Demultiplex cells by genotype against the variant database"
    )
    p.add_argument("--cellsnp-dir",   required=True,
                   help="cellSNP-lite output directory")
    p.add_argument("--db",            required=True,
                   help="SQLite database path")
    p.add_argument("--run-ids",       nargs="+", required=True,
                   help="run_ids of reference cell lines in the DB, or 'all'")
    p.add_argument("--output",        required=True,
                   help="output TSV path for cell assignments")
    p.add_argument("--demux-run-id",  default=None,
                   help="unique label for this demux run (required with --load-db)")
    p.add_argument("--load-db",       action="store_true",
                   help="write assignments into the DB cell_assignments table")
    p.add_argument("--min-depth",     type=int,   default=10,
                   help="min read depth at a position in a cell to include it (default: 10)")
    p.add_argument("--min-positions", type=int,   default=50,
                   help="min positions with coverage to attempt assignment (default: 50)")
    p.add_argument("--doublet-gap",   type=float, default=2.0,
                   help="min log-likelihood gap between top two donors for a singlet call (default: 2.0)")
    p.add_argument("--no-match-threshold", type=float, default=-0.5,
                   help="mean LL per position below which the call is no_match (default: -0.5)")
    return p.parse_args()


# ─────────────────────────────────────────────
# Load cellSNP-lite output
# ─────────────────────────────────────────────
def load_cellsnp_output(cellsnp_dir):
    """
    Load cellSNP-lite output directory.

    Returns
    -------
    positions : list of (chrom, pos, ref, alt)
    barcodes  : list of cell barcode strings
    ad        : scipy CSC sparse matrix, shape (n_positions, n_cells), alt read counts
    dp        : scipy CSC sparse matrix, shape (n_positions, n_cells), total depth
    """
    d = Path(cellsnp_dir)

    barcodes = (d / "cellSNP.samples.tsv").read_text().splitlines()

    vcf_path = d / "cellSNP.base.vcf"
    if not vcf_path.exists():
        vcf_path = d / "cellSNP.base.vcf.gz"
    positions = []
    vcf = VCF(str(vcf_path))
    for v in vcf:
        alt = str(v.ALT[0]) if v.ALT else "."
        positions.append((v.CHROM, v.POS, v.REF, alt))
    vcf.close()

    # CSC format: efficient column (per-cell) slicing
    ad = scipy.io.mmread(str(d / "cellSNP.tag.AD.mtx")).tocsc()
    dp = scipy.io.mmread(str(d / "cellSNP.tag.DP.mtx")).tocsc()

    return positions, barcodes, ad, dp


# ─────────────────────────────────────────────
# Query DB for reference genotypes
# ─────────────────────────────────────────────
def query_db_genotypes(db_path, run_ids, positions, min_depth=10):
    """
    Fetch genotypes for the given run_ids at the given positions.

    Returns
    -------
    donors  : list of run_id strings (column order)
    dosage  : float array, shape (n_positions, n_donors)
              NaN where a donor has no call at that position.
    """
    conn = sqlite3.connect(db_path)

    if run_ids == ["all"]:
        run_ids = [r[0] for r in conn.execute("SELECT run_id FROM runs").fetchall()]

    missing = [r for r in run_ids
               if not conn.execute("SELECT 1 FROM runs WHERE run_id=?", (r,)).fetchone()]
    if missing:
        print(
            f"\nWARNING: the following run_ids were not found in the database and will be skipped:\n"
            + "\n".join(f"  {r}" for r in missing)
            + "\nRun 'sqlite3 results/variants.db \"SELECT run_id FROM runs;\"' "
              "to list available run_ids.\n",
            file=sys.stderr,
        )
        run_ids = [r for r in run_ids if r not in missing]

    if not run_ids:
        conn.close()
        raise SystemExit(
            "\nERROR: no reference run_ids found in the database.\n"
            "Run bulk samples first to build the reference genotype DB,\n"
            "or check 'sqlite3 results/variants.db \"SELECT run_id FROM runs;\"'.\n"
        )

    pos_index = {p: i for i, p in enumerate(positions)}
    donor_idx = {r: i for i, r in enumerate(run_ids)}
    dosage    = np.full((len(positions), len(run_ids)), np.nan)

    placeholders = ",".join("?" * len(run_ids))
    rows = conn.execute(f"""
        SELECT r.run_id,
               v.chromosome, v.position, v.ref_allele, v.alt_allele,
               gc.genotype
        FROM   genotype_calls gc
        JOIN   variants v ON gc.variant_id = v.variant_id
        JOIN   runs     r ON gc.run_id     = r.run_id
        WHERE  r.run_id IN ({placeholders})
          AND  gc.genotype IS NOT NULL
          AND  gc.depth >= ?
    """, (*run_ids, min_depth)).fetchall()
    conn.close()

    for run_id, chrom, pos, ref, alt, gt in rows:
        key = (chrom, pos, ref, alt)
        if key not in pos_index or gt not in GT_DOSAGE:
            continue
        dosage[pos_index[key], donor_idx[run_id]] = GT_DOSAGE[gt]

    return run_ids, dosage


# ─────────────────────────────────────────────
# Scoring
# ─────────────────────────────────────────────
def _binom_logpmf(k, n, p):
    """Binomial log PMF, vectorised. Safe against p=0 or p=1."""
    p = np.clip(p, 1e-6, 1.0 - 1e-6)
    return (
        gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        + k * np.log(p)
        + (n - k) * np.log(1.0 - p)
    )


def score_cells(ad, dp, dosage, min_depth):
    """
    Compute per-(cell, donor) log-likelihood scores.

    Iterates over donors (not cells) so the sparse column structure of
    the CSC matrices is exploited efficiently without loading everything
    into a dense array.

    Returns
    -------
    scores  : float array (n_cells, n_donors)  — summed log-likelihoods
    n_used  : int array   (n_cells, n_donors)  — positions used per pair
    """
    n_pos, n_cells = ad.shape
    n_donors = dosage.shape[1]

    scores = np.full((n_cells, n_donors), np.nan)
    n_used = np.zeros((n_cells, n_donors), dtype=np.int32)

    for di in range(n_donors):
        p_all = dosage[:, di]          # (n_pos,) — NaN where no GT
        has_gt = ~np.isnan(p_all)
        if not has_gt.any():
            continue

        p_sub = p_all[has_gt]          # positions this donor has calls at
        pos_idx = np.where(has_gt)[0]

        # Slice the sparse matrices to only the rows with GT for this donor,
        # then convert to dense — much smaller than the full matrix.
        ad_sub = np.asarray(ad[pos_idx, :].todense())   # (n_gt_pos, n_cells)
        dp_sub = np.asarray(dp[pos_idx, :].todense())

        depth_ok = dp_sub >= min_depth                   # (n_gt_pos, n_cells)

        ll = _binom_logpmf(
            ad_sub,
            dp_sub,
            p_sub[:, np.newaxis],      # broadcast p over cells
        )
        ll[~depth_ok] = 0.0            # exclude low-depth positions

        scores[:, di] = ll.sum(axis=0)
        n_used[:, di] = depth_ok.sum(axis=0)

    return scores, n_used


# ─────────────────────────────────────────────
# Cell assignment
# ─────────────────────────────────────────────
def assign_cells(scores, n_used, donors, doublet_gap, no_match_threshold, min_positions):
    """
    Assign each cell to a donor or flag as doublet / no_match.

    Returns a list of dicts, one per cell.
    """
    results = []
    n_cells = scores.shape[0]

    for ci in range(n_cells):
        row  = scores[ci]
        used = n_used[ci]

        valid = ~np.isnan(row)
        max_positions = int(used[valid].max()) if valid.any() else 0

        if not valid.any() or max_positions < min_positions:
            results.append({
                "cell_idx":       ci,
                "assigned_line":  "no_match",
                "status":         "insufficient_coverage",
                "score":          None,
                "second_score":   None,
                "n_positions":    max_positions,
                "doublet":        False,
                "mean_ll_per_pos": None,
            })
            continue

        order        = np.argsort(row[valid])[::-1]
        valid_donors = np.where(valid)[0]
        best_di      = valid_donors[order[0]]
        best_score   = float(row[best_di])
        n_pos        = int(used[best_di])
        mean_ll      = best_score / n_pos if n_pos > 0 else None

        second_score = None
        doublet      = False
        second_name  = None

        if len(order) > 1:
            second_di    = valid_donors[order[1]]
            second_score = float(row[second_di])
            second_name  = donors[second_di]
            doublet      = (best_score - second_score) < doublet_gap

        if mean_ll is not None and mean_ll < no_match_threshold:
            status         = "no_match"
            assigned_line  = "no_match"
        elif doublet:
            status         = "doublet"
            assigned_line  = f"doublet:{donors[best_di]}+{second_name}"
        else:
            status         = "assigned"
            assigned_line  = donors[best_di]

        results.append({
            "cell_idx":        ci,
            "assigned_line":   assigned_line,
            "status":          status,
            "score":           round(best_score, 3),
            "second_score":    round(second_score, 3) if second_score is not None else None,
            "n_positions":     n_pos,
            "doublet":         doublet,
            "mean_ll_per_pos": round(mean_ll, 4) if mean_ll is not None else None,
        })

    return results


# ─────────────────────────────────────────────
# Output
# ─────────────────────────────────────────────
def write_results(results, barcodes, output_path):
    with open(output_path, "w") as f:
        f.write(
            "barcode\tassigned_line\tstatus\t"
            "score\tsecond_score\tn_positions\tdoublet\tmean_ll_per_pos\n"
        )
        for r in results:
            bc = barcodes[r["cell_idx"]]
            f.write(
                f"{bc}\t{r['assigned_line']}\t{r['status']}\t"
                f"{r['score'] if r['score'] is not None else ''}\t"
                f"{r['second_score'] if r['second_score'] is not None else ''}\t"
                f"{r['n_positions']}\t{r['doublet']}\t"
                f"{r['mean_ll_per_pos'] if r['mean_ll_per_pos'] is not None else ''}\n"
            )


def init_cell_assignments_table(conn):
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cell_assignments (
            assignment_id   INTEGER PRIMARY KEY AUTOINCREMENT,
            demux_run_id    TEXT    NOT NULL,
            cell_barcode    TEXT    NOT NULL,
            assigned_line   TEXT,
            status          TEXT,
            score           REAL,
            second_score    REAL,
            n_positions     INTEGER,
            doublet         INTEGER,
            mean_ll_per_pos REAL,
            UNIQUE(demux_run_id, cell_barcode)
        )
    """)
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_assignments_demux_run
            ON cell_assignments(demux_run_id)
    """)
    conn.commit()


def load_assignments_to_db(results, barcodes, db_path, demux_run_id):
    conn = sqlite3.connect(db_path)
    init_cell_assignments_table(conn)

    # Check for existing demux_run_id
    existing = conn.execute(
        "SELECT COUNT(*) FROM cell_assignments WHERE demux_run_id=?", (demux_run_id,)
    ).fetchone()[0]
    if existing:
        conn.close()
        raise SystemExit(
            f"\nERROR: demux_run_id '{demux_run_id}' already has {existing} rows "
            f"in cell_assignments.\nChoose a different --demux-run-id.\n"
        )

    rows = [
        (
            demux_run_id,
            barcodes[r["cell_idx"]],
            r["assigned_line"],
            r["status"],
            r["score"],
            r["second_score"],
            r["n_positions"],
            int(r["doublet"]),
            r["mean_ll_per_pos"],
        )
        for r in results
    ]
    conn.executemany("""
        INSERT INTO cell_assignments
            (demux_run_id, cell_barcode, assigned_line, status,
             score, second_score, n_positions, doublet, mean_ll_per_pos)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, rows)
    conn.commit()
    conn.close()
    print(f"  {len(rows):,} cell assignments loaded (demux_run_id={demux_run_id})")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    args = parse_args()

    if args.load_db and not args.demux_run_id:
        raise SystemExit("ERROR: --demux-run-id is required when using --load-db\n")

    print(f"\nLoading cellSNP-lite output from: {args.cellsnp_dir}")
    positions, barcodes, ad, dp = load_cellsnp_output(args.cellsnp_dir)
    print(f"  {len(positions):,} positions   {len(barcodes):,} cells")

    print(f"Querying DB for reference genotypes...")
    donors, dosage = query_db_genotypes(args.db, args.run_ids, positions, args.min_depth)

    print(f"  Using {len(donors)} reference line(s):")
    for d, n in zip(donors, (~np.isnan(dosage)).sum(axis=0)):
        print(f"    {d}: {n:,} positions with genotype call")

    n_shared = (~np.isnan(dosage)).any(axis=1).sum()
    print(f"  {n_shared:,} positions covered by at least one reference line")

    print(f"\nScoring {len(barcodes):,} cells against {len(donors)} reference line(s)...")
    scores, n_used = score_cells(ad, dp, dosage, args.min_depth)

    print("Assigning cells...")
    results = assign_cells(
        scores, n_used, donors,
        args.doublet_gap,
        args.no_match_threshold,
        args.min_positions,
    )

    # Summary
    status_counts = Counter(r["status"] for r in results)
    print(f"\nResults ({len(results):,} cells total):")
    for status, count in sorted(status_counts.items()):
        pct = 100 * count / len(results)
        print(f"  {status:<25} {count:>6,}  ({pct:.1f}%)")

    assigned = [r for r in results if r["status"] == "assigned"]
    if assigned:
        print("\n  Per cell line:")
        for line, count in Counter(r["assigned_line"] for r in assigned).most_common():
            print(f"    {line:<35} {count:>6,}")

    print(f"\nWriting {args.output}...")
    write_results(results, barcodes, args.output)

    if args.load_db:
        print(f"Loading assignments into DB ({args.db})...")
        load_assignments_to_db(results, barcodes, args.db, args.demux_run_id)

    print("\nDone.\n")


if __name__ == "__main__":
    main()
