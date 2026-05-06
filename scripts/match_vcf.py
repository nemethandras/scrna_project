"""
Match an unknown VCF against all samples in the database.

For each sample, calculates what percentage of the query VCF's variants
are also present in that sample at the same position with the same allele.
Reports all matches above a minimum threshold.

Usage:
    python scripts/match_vcf.py --vcf path/to/unknown.vcf
    python scripts/match_vcf.py --vcf path/to/unknown.vcf --min-match 30
    python scripts/match_vcf.py --vcf path/to/unknown.vcf --metric jaccard
"""

import argparse
import sqlite3
from cyvcf2 import VCF


def parse_args():
    p = argparse.ArgumentParser(description="Match unknown VCF against database samples")
    p.add_argument("--vcf",       required=True, help="path to the unknown VCF file")
    p.add_argument("--db",        default="results/variants.db", help="SQLite database path")
    p.add_argument("--min-match", type=float, default=40.0,
                   help="minimum match %% to report (default: 40)")
    p.add_argument("--metric",    choices=["overlap", "jaccard"], default="overlap",
                   help="overlap: shared/query total; jaccard: shared/union (default: overlap)")
    return p.parse_args()


def load_query_variants(vcf_path):
    """Parse VCF and return a set of (chrom, pos, ref, alt) tuples."""
    variants = set()
    for v in VCF(vcf_path):
        if v.ALT:
            variants.add((v.CHROM, v.POS, v.REF, str(v.ALT[0])))
    return variants


def load_db_variants(conn):
    """Return a dict of run_id -> set of (chrom, pos, ref, alt) tuples."""
    db_variants = {}
    rows = conn.execute("""
        SELECT gc.run_id, v.chromosome, v.position, v.ref_allele, v.alt_allele
        FROM genotype_calls gc
        JOIN variants v ON gc.variant_id = v.variant_id
    """)
    for run_id, chrom, pos, ref, alt in rows:
        db_variants.setdefault(run_id, set()).add((chrom, pos, ref, alt))
    return db_variants


def match_score(query, sample, metric):
    shared = len(query & sample)
    if metric == "overlap":
        return round(shared / len(query) * 100, 2) if query else 0.0
    else:
        union = len(query | sample)
        return round(shared / union * 100, 2) if union else 0.0


def main():
    args = parse_args()

    print(f"\nLoading query VCF: {args.vcf}")
    query = load_query_variants(args.vcf)
    print(f"  {len(query):,} variants in query")

    conn = sqlite3.connect(args.db)

    # Fetch sample name lookup
    sample_names = {
        r[0]: r[1] for r in conn.execute(
            "SELECT r.run_id, s.name FROM runs r JOIN samples s ON r.sample_id = s.sample_id"
        )
    }

    print(f"  Loading database variants...")
    db_variants = load_db_variants(conn)
    conn.close()

    print(f"\nMatching against {len(db_variants)} samples "
          f"(metric: {args.metric}, min: {args.min_match}%)\n")

    results = []
    for run_id, sample_variants in db_variants.items():
        score = match_score(query, sample_variants, args.metric)
        shared = len(query & sample_variants)
        results.append((score, run_id, sample_names.get(run_id, "?"), shared, len(sample_variants)))

    results.sort(reverse=True)

    header = f"  {'Run':<25}  {'Sample':<15}  {'Match':>7}  {'Shared':>8}  {'In DB':>8}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    found = False
    for score, run_id, sample_name, shared, db_total in results:
        if score >= args.min_match:
            print(f"  {run_id:<25}  {sample_name:<15}  {score:>6.1f}%  {shared:>8,}  {db_total:>8,}")
            found = True

    if not found:
        print(f"  No matches found above {args.min_match}% threshold.")

    print(f"\nQuery contained {len(query):,} variants total.\n")


if __name__ == "__main__":
    main()
