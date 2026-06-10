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
    p.add_argument("--min-concordance", type=float, default=0.4, metavar="F",
                   help="min ALT-recall to accept a match (default: 0.4)")
    p.add_argument("--min-gap", type=float, default=0.15, metavar="F",
                   help="min gap between best and second-best concordance to accept a match "
                        "(default: 0.15)")
    p.add_argument("--min-positions", type=int, default=50,
                   help="min shared positions to attempt matching (default: 100)")
    return p.parse_args()


def load_donor_genotypes(vireo_dir):
    """
    Parse GT_donors.vcf.gz from Vireo output.

    Returns
    -------
    donors    : list of donor name strings
    positions : list of (chrom, pos, ref, alt) tuples
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
        alt = str(v.ALT[0]) if v.ALT else "."
        positions.append((v.CHROM, v.POS, v.REF, alt))
        row = []
        for gt in v.genotypes:
            a1, a2 = gt[0], gt[1]
            row.append(GT_MAP.get((a1, a2), -1))
        rows.append(row)

    vcf.close()

    dosage = np.array(rows, dtype=np.int8)  # (n_positions, n_donors)
    print(f"  {len(donors)} Vireo donors, {len(positions):,} positions")
    return donors, positions, dosage


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

    pos_index = {p: i for i, p in enumerate(positions)}
    GT_MAP = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
    dosage = np.full((len(positions), len(run_ids)), -1, dtype=np.int8)

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
    """, run_ids).fetchall()
    conn.close()

    run_idx = {r: i for i, r in enumerate(run_ids)}
    for run_id, chrom, pos, ref, alt, gt in rows:
        key = (chrom, pos, ref, alt)
        if key in pos_index and gt in GT_MAP:
            dosage[pos_index[key], run_idx[run_id]] = GT_MAP[gt]

    print(f"  {len(run_ids)} reference lines in DB")
    return run_ids, dosage


def compute_concordance(vireo_dosage, db_dosage):
    """
    Compute ALT-recall between each Vireo donor and each DB reference.

    Metric: among positions where the DB has a non-ref call (0/1 or 1/1),
    what fraction does Vireo also call as non-ref?

    This is robust to the old DB samples that only contain ALT calls
    (no 0/0 entries), and handles sparse scRNA-seq coverage — Vireo will
    default to 0/0 at uncovered positions, but the signal comes from the
    positions where it does infer a variant.

    Returns
    -------
    concordance : float array (n_donors, n_runs)  — ALT recall
    n_shared    : int array   (n_donors, n_runs)  — DB ALT positions used
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

            # Focus on positions where the DB has a non-ref call
            db_alt = r_valid & (rd > 0)
            shared = v_valid & db_alt
            n = int(shared.sum())
            if n == 0:
                continue
            # Fraction of DB ALT positions where Vireo also calls non-ref
            vireo_also_alt = int((vd[shared] > 0).sum())
            concordance[di, ri] = vireo_also_alt / n
            n_shared[di, ri]    = n

    return concordance, n_shared


def write_matches(donors, run_ids, concordance, n_shared, output_path,
                  min_positions, min_concordance, min_gap):
    with open(output_path, "w") as f:
        f.write(
            "vireo_donor\tassigned_line\tconcordance\tn_positions"
            "\tsecond_line\tsecond_concordance\tgap\n"
        )
        for di, donor in enumerate(donors):
            row = concordance[di]
            ns  = n_shared[di]

            valid = ~np.isnan(row) & (ns >= min_positions)

            if not valid.any():
                f.write(f"{donor}\tno_match\t\t0\t\t\t\n")
                continue

            order    = np.argsort(row[valid])[::-1]
            vi       = np.where(valid)[0]
            best_ri  = vi[order[0]]
            best_c   = float(row[best_ri])
            best_n   = int(ns[best_ri])
            best_run = run_ids[best_ri]

            second_run = ""
            second_c   = np.nan
            if len(order) > 1:
                second_ri  = vi[order[1]]
                second_run = run_ids[second_ri]
                second_c   = float(row[second_ri])

            gap = best_c - second_c if not np.isnan(second_c) else best_c

            # Reject if below absolute threshold or gap is too small
            if best_c < min_concordance or gap < min_gap:
                assigned = "no_match"
            else:
                assigned = best_run

            second_c_str = f"{second_c:.4f}" if not np.isnan(second_c) else ""
            f.write(
                f"{donor}\t{assigned}\t{best_c:.4f}\t{best_n}"
                f"\t{second_run}\t{second_c_str}\t{gap:.4f}\n"
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

    print(f"\nWriting {args.output}...")
    write_matches(donors, run_ids, concordance, n_shared, args.output,
                  args.min_positions, args.min_concordance, args.min_gap)
    print("Done.\n")


if __name__ == "__main__":
    main()
