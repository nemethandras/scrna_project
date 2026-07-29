"""
Match Vireo's anonymous donor genotypes to cell line identities in the DB.

Vireo clusters cells into donors (donor0, donor1, ...) and infers a genotype
for each donor in GT_donors.vcf.gz. This script compares those inferred
genotypes against the reference genotypes in the DB to identify which cell
line each donor is.

Usage:
    python scripts/match_vireo.py \
        --vireo-dir results/demux/pool_ctr_001/vireo \
        --db results/variants.db \
        --output results/demux/pool_ctr_001/donor_matches.tsv
"""

import sys
import argparse
import sqlite3
import numpy as np
from pathlib import Path
from cyvcf2 import VCF
from scipy.optimize import linear_sum_assignment


def parse_args():
    p = argparse.ArgumentParser(
        description="Match Vireo donor genotypes to cell lines in the DB"
    )
    p.add_argument("--vireo-dir",  required=True,
                   help="Vireo output directory (must contain GT_donors.vcf.gz)")
    p.add_argument("--db",         required=True,
                   help="SQLite database path")
    p.add_argument("--run-ids",    nargs="+", default=["all"],
                   help="run_ids to compare against, or 'all' (default: all)")
    p.add_argument("--output",     required=True,
                   help="output TSV path for donor→cell_line matches")
    p.add_argument("--min-concordance", type=float, default=0.65, metavar="F",
                   help="min genotype concordance to accept a match (default: 0.65)")
    p.add_argument("--min-gap", type=float, default=0.10, metavar="F",
                   help="min gap between best and second-best concordance to accept a match "
                        "(default: 0.10)")
    p.add_argument("--min-positions", type=int, default=200,
                   help="min shared positions to attempt matching (default: 200)")
    p.add_argument("--unique", action="store_true",
                   help="enforce one-to-one donor↔reference matching via Hungarian algorithm "
                        "(prevents multiple donors from matching the same cell line)")
    return p.parse_args()


def load_donor_genotypes(vireo_dir):
    """
    Parse GT_donors.vcf.gz from Vireo output.

    Returns
    -------
    donors    : list of donor name strings
    positions : list of (chrom, pos) tuples
    dosage    : int array (n_positions, n_donors) — 0/1/2, -1 = no call
    """
    d = Path(vireo_dir)
    for candidate in ["GT_donors.vireo.vcf.gz", "GT_donors.vcf.gz", "GT_donors.vcf"]:
        vcf_path = d / candidate
        if vcf_path.exists():
            break
    else:
        raise SystemExit(
            f"\nERROR: GT_donors VCF not found in {vireo_dir}\n"
            "Re-run Vireo without -d (genotype-free mode) to produce inferred donor genotypes.\n"
        )

    vcf = VCF(str(vcf_path))
    donors = list(vcf.samples)
    positions = []
    rows = []

    GT_MAP = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 2}

    for v in vcf:
        positions.append((v.CHROM, v.POS))
        row = []
        for gt in v.genotypes:
            a1, a2 = gt[0], gt[1]
            row.append(GT_MAP.get((a1, a2), -1))
        rows.append(row)

    vcf.close()

    dosage = np.array(rows, dtype=np.int8)  # (n_positions, n_donors)
    print(f"  {len(donors)} Vireo donors, {len(positions):,} positions")
    return donors, positions, dosage


def load_cell_line_names(db_path, run_ids):
    """Return dict run_id → cell_line (falls back to sample name, then run_id)."""
    conn = sqlite3.connect(db_path)
    names = {}
    for run_id in run_ids:
        row = conn.execute("""
            SELECT s.cell_line, s.name
            FROM   runs r JOIN samples s ON r.sample_id = s.sample_id
            WHERE  r.run_id = ?
        """, (run_id,)).fetchone()
        if row:
            cl, sname = row
            names[run_id] = cl if cl else sname
        else:
            names[run_id] = run_id
    conn.close()
    return names


def load_db_genotypes(db_path, run_ids, positions):
    """
    Fetch reference genotypes from the DB for the given positions.

    Returns
    -------
    run_ids : list of run_id strings (filtered to those present in DB)
    dosage  : int array (n_positions, n_runs) — 0/1/2, -1 = no call
    """
    conn = sqlite3.connect(db_path)

    if run_ids == ["all"]:
        run_ids = [r[0] for r in conn.execute("SELECT run_id FROM runs").fetchall()]

    missing = [r for r in run_ids
               if not conn.execute("SELECT 1 FROM runs WHERE run_id=?", (r,)).fetchone()]
    if missing:
        print(f"WARNING: {len(missing)} run_ids not found in DB, skipping: {missing}",
              file=sys.stderr)
        run_ids = [r for r in run_ids if r not in missing]
    if not run_ids:
        conn.close()
        raise SystemExit("\nERROR: no reference run_ids found in the database.\n")

    # Key on (chrom, pos) only — not (chrom, pos, ref, alt) — so that 0/0 DB calls
    # (stored with alt_allele='.') match Vireo positions that carry a real ALT allele
    # (because another donor at that site has the ALT). The SNP panel has unique positions.
    pos_index = {(p[0], p[1]): i for i, p in enumerate(positions)}
    GT_MAP = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
    dosage = np.full((len(positions), len(run_ids)), -1, dtype=np.int8)

    placeholders = ",".join("?" * len(run_ids))
    rows = conn.execute(f"""
        SELECT r.run_id,
               v.chromosome, v.position,
               gc.genotype
        FROM   genotype_calls gc
        JOIN   variants v ON gc.variant_id = v.variant_id
        JOIN   runs     r ON gc.run_id     = r.run_id
        WHERE  r.run_id IN ({placeholders})
          AND  gc.genotype IS NOT NULL
    """, run_ids).fetchall()
    conn.close()

    run_idx = {r: i for i, r in enumerate(run_ids)}
    for run_id, chrom, pos, gt in rows:
        key = (chrom, pos)
        if key in pos_index and gt in GT_MAP:
            idx = pos_index[key]
            ri  = run_idx[run_id]
            # Keep the ALT call if one exists (prefer non-ref over 0/0 for multiallelic)
            if dosage[idx, ri] < 0 or GT_MAP[gt] > 0:
                dosage[idx, ri] = GT_MAP[gt]

    print(f"  {len(run_ids)} reference lines in DB")
    return run_ids, dosage


def compute_concordance(vireo_dosage, db_dosage):
    """
    Compute binary ALT-presence concordance between each Vireo donor and each DB reference.

    Metric: among all positions where BOTH the Vireo donor and the DB reference have a
    valid call (dosage ≥ 0), what fraction agree on whether there is ANY alt allele?
    That is: (dosage > 0) in Vireo == (dosage > 0) in DB.

    This is variant-count-agnostic: the denominator includes 0/0 positions (for DB
    entries that store them), so a reference with fewer variants is no longer
    artificially preferred. At the same time, ignoring the het/hom distinction
    (0/1 vs 1/1 both count as "ALT present") avoids penalising Vireo's inference noise
    from sparse scRNA coverage, where truly homozygous-alt positions are often called
    heterozygous due to read-count variability.

    Graceful degradation: for DB entries that only store ALT calls (no 0/0 rows, e.g.
    samples loaded under the old variant-only pipeline), the shared positions are purely
    ALT positions and the metric reduces to the original ALT-recall.

    Returns
    -------
    concordance : float array (n_donors, n_runs)  — binary ALT-presence concordance
    n_shared    : int array   (n_donors, n_runs)  — positions with valid calls in both
    """
    n_pos, n_donors = vireo_dosage.shape
    n_runs = db_dosage.shape[1]

    concordance = np.full((n_donors, n_runs), np.nan)
    n_shared    = np.zeros((n_donors, n_runs), dtype=np.int32)

    for di in range(n_donors):
        vd = vireo_dosage[:, di]
        v_valid = vd >= 0

        for ri in range(n_runs):
            rd = db_dosage[:, ri]
            r_valid = rd >= 0

            # All positions where both have a valid call
            shared = v_valid & r_valid
            n = int(shared.sum())
            if n == 0:
                continue
            # Binary agreement: does each position carry an ALT allele?
            v_has_alt = vd[shared] > 0
            r_has_alt = rd[shared] > 0
            matching = int((v_has_alt == r_has_alt).sum())
            concordance[di, ri] = matching / n
            n_shared[di, ri]    = n

    return concordance, n_shared


def unique_match(donors, run_ids, concordance, n_shared, min_positions):
    """
    Hungarian 1:1 matching: each reference can be assigned to at most one donor.
    Returns assigned_ri[di] = reference index for each donor, or -1 if no valid positions.
    """
    n_donors = len(donors)
    n_refs   = len(run_ids)

    # Mask positions below min_positions as -inf so they can't be chosen
    cost = concordance.copy()
    invalid = np.isnan(cost) | (n_shared < min_positions)
    cost[invalid] = 0.0  # 0 = worst (we maximize)

    # linear_sum_assignment minimizes; negate to maximize concordance
    if n_donors <= n_refs:
        row_ind, col_ind = linear_sum_assignment(-cost)
    else:
        # More donors than refs: transpose and solve
        row_ind, col_ind = linear_sum_assignment(-cost.T)
        row_ind, col_ind = col_ind, row_ind

    assigned_ri = np.full(n_donors, -1, dtype=int)
    for di, ri in zip(row_ind, col_ind):
        if not invalid[di, ri]:
            assigned_ri[di] = ri
    return assigned_ri


def write_matches(donors, run_ids, concordance, n_shared, output_path,
                  min_positions, min_concordance, min_gap, use_unique=False,
                  cell_line_names=None):
    # Compute 1:1 assignments if requested
    unique_ri = unique_match(donors, run_ids, concordance, n_shared, min_positions) \
                if use_unique else None

    cl = cell_line_names or {}

    with open(output_path, "w") as f:
        f.write(
            "vireo_donor\tassigned_line\tcell_line\tconcordance\tn_positions"
            "\tsecond_line\tsecond_concordance\tgap\tconfidence\n"
        )
        for di, donor in enumerate(donors):
            row = concordance[di]
            ns  = n_shared[di]

            valid = ~np.isnan(row) & (ns >= min_positions)

            if not valid.any():
                f.write(f"{donor}\tno_match\t\t\t0\t\t\t\tno_data\n")
                continue

            # In unique mode use the Hungarian-assigned reference, else pick best globally
            if use_unique and unique_ri is not None:
                best_ri = unique_ri[di]
                if best_ri == -1:
                    f.write(f"{donor}\tno_match\t\t\t0\t\t\t\tno_data\n")
                    continue
            else:
                order   = np.argsort(row[valid])[::-1]
                vi      = np.where(valid)[0]
                best_ri = vi[order[0]]

            best_c   = float(row[best_ri])
            best_n   = int(ns[best_ri])
            best_run = run_ids[best_ri]

            # Second-best: best valid index that is NOT best_ri
            other_valid = np.where(valid & (np.arange(len(run_ids)) != best_ri))[0]
            if len(other_valid):
                second_ri  = other_valid[np.argmax(row[other_valid])]
                second_run = run_ids[second_ri]
                second_c   = float(row[second_ri])
            else:
                second_run = ""
                second_c   = np.nan

            gap = best_c - second_c if not np.isnan(second_c) else best_c

            # Determine assignment and confidence
            if best_c >= min_concordance and gap >= min_gap:
                assigned   = best_run
                confidence = "high"
            elif best_c >= min_concordance * 0.7:
                # Below threshold but still best available — flag as low_confidence
                assigned   = best_run
                confidence = "low"
            else:
                assigned   = "no_match"
                confidence = "no_match"

            second_c_str = f"{second_c:.4f}" if not np.isnan(second_c) else ""
            assigned_cl  = cl.get(assigned, assigned) if assigned != "no_match" else ""
            f.write(
                f"{donor}\t{assigned}\t{assigned_cl}\t{best_c:.4f}\t{best_n}"
                f"\t{second_run}\t{second_c_str}\t{gap:.4f}\t{confidence}\n"
            )


def main():
    args = parse_args()

    print(f"\nLoading Vireo donor genotypes from: {args.vireo_dir}")
    donors, positions, vireo_dosage = load_donor_genotypes(args.vireo_dir)

    print(f"Loading DB reference genotypes...")
    run_ids, db_dosage = load_db_genotypes(args.db, args.run_ids, positions)

    n_shared_total = ((vireo_dosage >= 0).any(axis=1) & (db_dosage >= 0).any(axis=1)).sum()
    print(f"  {n_shared_total:,} positions present in both Vireo and DB")

    print(f"Computing concordance ({len(donors)} donors × {len(run_ids)} references)...")
    concordance, n_shared = compute_concordance(vireo_dosage, db_dosage)

    print(f"\nDonor → cell line matches:")
    for di, donor in enumerate(donors):
        row  = concordance[di]
        ns   = n_shared[di]
        valid = ~np.isnan(row) & (ns >= args.min_positions)
        if not valid.any():
            print(f"  {donor:<12}  no_match  (no shared positions)")
            continue
        best_ri = int(np.nanargmax(row * valid))
        print(f"  {donor:<12}  {run_ids[best_ri]:<35}  concordance={row[best_ri]:.3f}  n={ns[best_ri]:,}")

    cell_line_names = load_cell_line_names(args.db, run_ids)

    if args.unique:
        print("Using 1:1 Hungarian matching (--unique)")

    print(f"\nWriting {args.output}...")
    write_matches(donors, run_ids, concordance, n_shared, args.output,
                  args.min_positions, args.min_concordance, args.min_gap,
                  use_unique=args.unique, cell_line_names=cell_line_names)
    print("Done.\n")


if __name__ == "__main__":
    main()
