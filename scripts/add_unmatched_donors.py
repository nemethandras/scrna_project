"""
Add unmatched Vireo donors to the variants DB as new reference samples.

For each donor in donor_matches.tsv whose confidence is 'no_match' or 'no_data',
reads the inferred genotype from GT_donors.vireo.vcf.gz and inserts it into the DB
as a new sample so future pools can match against it.

New entries:
  samples.name        = unknown_<demux_run_id>_<donor_id>
  samples.cell_line   = NULL  (fill in once you identify the cell line)
  runs.run_id         = vireo_<demux_run_id>_<donor_id>
  runs.reference      = vireo_inferred
  genotype_calls      = non-ref (0/1, 1/1) calls only — consistent with ALT-only DB entries

Re-running is safe: INSERT OR IGNORE prevents duplicates.

Usage:
    python scripts/add_unmatched_donors.py \\
        --vireo-dir    results/demux/pool_ctr_v4/vireo \\
        --matches      results/demux/pool_ctr_v4/donor_matches.tsv \\
        --db           results/variants.db \\
        --demux-run-id pool_ctr_v4

After identifying the cell line, update the DB:
    UPDATE samples SET cell_line='IS3'
    WHERE name = 'unknown_pool_ctr_v4_donor3';
"""

import argparse
import csv
import datetime
import sqlite3
import sys
from pathlib import Path

from cyvcf2 import VCF


def parse_args():
    p = argparse.ArgumentParser(
        description="Insert unmatched Vireo donors into the DB as new reference samples"
    )
    p.add_argument("--vireo-dir",    required=True, metavar="DIR",
                   help="Vireo output directory (must contain GT_donors.vireo.vcf.gz)")
    p.add_argument("--matches",      required=True, metavar="TSV",
                   help="donor_matches.tsv from match_vireo.py")
    p.add_argument("--db",           required=True, metavar="PATH",
                   help="SQLite database path")
    p.add_argument("--demux-run-id", required=True, metavar="ID",
                   help="demux run label used to build the new sample/run names")
    return p.parse_args()


def load_unmatched_donors(matches_path):
    """Return list of donor names with confidence 'no_match' or 'no_data'."""
    unmatched = []
    with open(matches_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("confidence", "") in ("no_match", "no_data"):
                unmatched.append(row["vireo_donor"])
    return unmatched


def load_vireo_vcf(vireo_dir):
    """
    Parse GT_donors.vireo.vcf.gz.

    Returns
    -------
    donors  : list of donor name strings
    records : list of (chrom, pos, ref, alt, {donor: gt_str})
              Only positions where at least one donor has a non-ref call.
    """
    d = Path(vireo_dir)
    for candidate in ["GT_donors.vireo.vcf.gz", "GT_donors.vcf.gz", "GT_donors.vcf"]:
        vcf_path = d / candidate
        if vcf_path.exists():
            break
    else:
        raise SystemExit(
            f"\nERROR: GT_donors VCF not found in {vireo_dir}\n"
            "Vireo must be run in genotype-free mode to produce inferred donor genotypes.\n"
        )

    vcf = VCF(str(vcf_path))
    donors = list(vcf.samples)
    GT_STR = {(0, 0): "0/0", (0, 1): "0/1", (1, 0): "0/1", (1, 1): "1/1"}

    records = []
    for v in vcf:
        alt = str(v.ALT[0]) if v.ALT else "."
        gts = {}
        for donor, gt in zip(donors, v.genotypes):
            a1, a2 = gt[0], gt[1]
            s = GT_STR.get((a1, a2))
            if s:
                gts[donor] = s
        if gts:
            records.append((v.CHROM, v.POS, v.REF, alt, gts))

    vcf.close()
    print(f"  {len(donors)} donors, {len(records):,} positions in VCF")
    return donors, records


def _variant_type(ref, alt):
    return "SNP" if len(ref) == 1 and len(alt) == 1 else "INDEL"


def insert_donor(conn, demux_run_id, donor, records):
    """
    Insert one unmatched donor as a new sample + run + genotype_calls.
    Only non-ref calls (0/1, 1/1) are inserted.
    Returns (sample_name, run_id, n_calls_inserted).
    """
    sample_name = f"unknown_{demux_run_id}_{donor}"
    run_id      = f"vireo_{demux_run_id}_{donor}"
    today       = datetime.date.today().isoformat()

    conn.execute(
        "INSERT OR IGNORE INTO samples (name, cell_line, date_added) VALUES (?, NULL, ?)",
        (sample_name, today),
    )
    sample_id = conn.execute(
        "SELECT sample_id FROM samples WHERE name=?", (sample_name,)
    ).fetchone()[0]

    conn.execute(
        """INSERT OR IGNORE INTO runs
           (run_id, sample_id, reference, sequencing, run_date,
            mapping_rate, total_raw, total_filt)
           VALUES (?, ?, 'vireo_inferred', '10x_scrna', ?, NULL, NULL, NULL)""",
        (run_id, sample_id, today),
    )

    n_inserted = 0
    for chrom, pos, ref, alt, gts in records:
        gt_str = gts.get(donor)
        if not gt_str or gt_str == "0/0":
            continue

        conn.execute(
            """INSERT OR IGNORE INTO variants
               (chromosome, position, ref_allele, alt_allele, variant_type)
               VALUES (?, ?, ?, ?, ?)""",
            (chrom, pos, ref, alt, _variant_type(ref, alt)),
        )
        row = conn.execute(
            """SELECT variant_id FROM variants
               WHERE chromosome=? AND position=? AND ref_allele=? AND alt_allele=?""",
            (chrom, pos, ref, alt),
        ).fetchone()
        if row is None:
            continue
        variant_id = row[0]

        af = 0.5 if gt_str == "0/1" else 1.0
        conn.execute(
            """INSERT OR IGNORE INTO genotype_calls
               (run_id, variant_id, genotype, quality, depth,
                ref_depth, alt_depth, allele_freq, filter_status)
               VALUES (?, ?, ?, NULL, NULL, NULL, NULL, ?, 'PASS')""",
            (run_id, variant_id, gt_str, af),
        )
        n_inserted += 1

    return sample_name, run_id, n_inserted


def main():
    args = parse_args()

    print(f"Loading donor matches from: {args.matches}")
    unmatched = load_unmatched_donors(args.matches)
    if not unmatched:
        print("No unmatched donors found — nothing to add.")
        return
    print(f"  Unmatched donors: {unmatched}")

    print(f"Loading Vireo VCF from: {args.vireo_dir}")
    vcf_donors, records = load_vireo_vcf(args.vireo_dir)

    missing = [d for d in unmatched if d not in vcf_donors]
    if missing:
        print(f"  WARNING: {missing} not found in VCF sample list — skipping")
        unmatched = [d for d in unmatched if d in vcf_donors]
    if not unmatched:
        print("Nothing left to insert.")
        return

    conn = sqlite3.connect(args.db)
    try:
        for donor in unmatched:
            sample_name, run_id, n = insert_donor(conn, args.demux_run_id, donor, records)
            print(f"  {donor:<12} → {sample_name}  ({n:,} ALT calls, run_id={run_id})")
        conn.commit()
    finally:
        conn.close()

    print(f"\nDone. To label the new entries once identified:")
    print(f"  UPDATE samples SET cell_line='<NAME>'")
    print(f"  WHERE name LIKE 'unknown_{args.demux_run_id}_%';")


if __name__ == "__main__":
    main()
