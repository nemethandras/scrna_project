"""
Build a multi-sample donor genotype VCF from DB references for reference-guided Vireo.

When the cell lines in a pool are known in advance, this script extracts their
genotypes from the DB and writes a VCF that Vireo uses as a reference:

    vireo -c cellsnp_dir -d donor_ref.vcf -o vireo_dir

This is more accurate than genotype-free mode because Vireo uses known genotype
priors instead of inferring them from sparse single-cell data.  Small clusters
and genetically similar lines are assigned correctly even when Vireo alone would
merge or mis-separate them.

When a cell line is not yet in the DB, supply an external VCF with --extra-vcf;
its sample columns are merged with the DB-derived columns.

Usage:
    # From cell line names (most recent run per line taken from DB):
    python scripts/make_donor_vcf.py \\
        --cell-lines RKO HCT116 HT29 LIM1215 CACO2 HCT8 DLD1 \\
        --db results/variants.db \\
        --output results/demux/pool_ctr_v5/donor_ref.vcf

    # From explicit run_ids (column label = cell_line name or run_id):
    python scripts/make_donor_vcf.py \\
        --run-ids SRR5071686_hg38 SRR5071667_hg38 SRR5071672_hg38 \\
        --db results/variants.db \\
        --output donor_ref.vcf

    # Mix DB lines + external VCF for a line not yet in the DB:
    python scripts/make_donor_vcf.py \\
        --cell-lines RKO HCT116 HT29 \\
        --extra-vcf path/to/new_line.vcf \\
        --db results/variants.db \\
        --output donor_ref.vcf
"""

import argparse
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

_CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})

_GT_NORMALISE = {
    "0/0": "0/0", "0|0": "0/0",
    "0/1": "0/1", "1/0": "0/1", "0|1": "0/1", "1|0": "0/1",
    "1/1": "1/1", "1|1": "1/1",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Build multi-sample donor VCF from DB genotypes for reference-guided Vireo"
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--cell-lines", nargs="+", metavar="NAME",
        help="cell line names to fetch from DB (most recent run per line is used)"
    )
    src.add_argument(
        "--run-ids", nargs="+", metavar="ID",
        help="explicit run_ids from the DB (one VCF column per run; label = cell_line name)"
    )
    p.add_argument("--db",        required=True, metavar="PATH", help="SQLite database path")
    p.add_argument("--extra-vcf", metavar="PATH",
                   help="VCF with additional donor columns for lines not in the DB; "
                        "sample names from this file become extra VCF columns")
    p.add_argument("--output",    required=True, metavar="PATH", help="output VCF path")
    return p.parse_args()


def resolve_run_ids(conn, cell_lines):
    """Return {cell_line: run_id}, picking the most recent run for each line."""
    resolved = {}
    missing  = []
    for cl in cell_lines:
        row = conn.execute(
            """
            SELECT r.run_id
            FROM runs r
            JOIN samples s ON r.sample_id = s.sample_id
            WHERE LOWER(s.cell_line) = LOWER(?)
            ORDER BY COALESCE(r.run_date, '') DESC, r.run_id DESC
            LIMIT 1
            """,
            (cl,),
        ).fetchone()
        if row:
            resolved[cl] = row[0]
        else:
            missing.append(cl)
    if missing:
        print(
            f"WARNING: the following cell lines were not found in the DB and will be skipped:\n"
            f"  {missing}\n"
            f"  Use --extra-vcf to supply genotypes for lines not in the DB.",
            file=sys.stderr,
        )
    return resolved


def load_db_genotypes(conn, run_id_to_label):
    """
    Fetch all genotype calls for the given run_ids.

    Two-pass approach so that 0/0 (ref-homozygous) calls are included for
    positions where another donor carries an ALT allele.  The DB stores 0/0
    calls as separate variant records with alt_allele='.'; a naive query
    filtered to ALT variants misses them and writes './' (missing) instead,
    stripping Vireo of the discriminating signal that a sample is ref-homozygous.

    Pass 1: collect every (chrom, pos) where any queried run has an ALT call.
    Pass 2: for those positions, fetch ALL calls for the queried runs —
            including 0/0 records stored under the alt_allele='.' variant.

    Returns
    -------
    pos_meta : dict  (chrom, pos) → (ref, alt)
    gts      : dict  (chrom, pos) → {label: gt_str}
    """
    run_ids      = list(run_id_to_label.keys())
    placeholders = ",".join("?" * len(run_ids))

    # Pass 1 — positions with a real ALT call in any queried run
    alt_rows = conn.execute(
        f"""
        SELECT DISTINCT v.chromosome, v.position, v.ref_allele, v.alt_allele
        FROM genotype_calls gc
        JOIN variants v ON gc.variant_id = v.variant_id
        WHERE gc.run_id IN ({placeholders})
          AND v.alt_allele NOT IN ('.', '*')
        """,
        run_ids,
    ).fetchall()

    pos_meta     = {}
    alt_positions = set()
    for chrom, pos, ref, alt in alt_rows:
        key = (chrom, pos)
        if key not in pos_meta:
            pos_meta[key] = (ref, alt)
        alt_positions.add(key)

    # Pass 2 — all calls (ALT + 0/0) for the queried runs; filter to ALT positions
    all_rows = conn.execute(
        f"""
        SELECT v.chromosome, v.position, gc.run_id, gc.genotype
        FROM genotype_calls gc
        JOIN variants v ON gc.variant_id = v.variant_id
        WHERE gc.run_id IN ({placeholders})
        """,
        run_ids,
    ).fetchall()

    gts = defaultdict(dict)
    for chrom, pos, run_id, genotype in all_rows:
        key = (chrom, pos)
        if key not in alt_positions:
            continue
        label = run_id_to_label[run_id]
        gt    = _GT_NORMALISE.get(genotype, "./.")
        # Prefer ALT calls over 0/0 (guard against duplicate rows)
        existing = gts[key].get(label)
        if existing is None or (existing == "0/0" and gt not in ("0/0", "./.")):
            gts[key][label] = gt

    return pos_meta, gts


def load_extra_vcf(vcf_path):
    """
    Parse an extra multi-sample VCF with cyvcf2.

    Returns
    -------
    extra_samples : list of sample name strings
    extra_pos_meta: dict (chrom, pos) → (ref, alt)   — only new positions not in DB
    extra_gts     : dict (chrom, pos) → {sample: gt_str}
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise SystemExit(
            "cyvcf2 is required to read --extra-vcf; install with:\n  pip install cyvcf2"
        )

    GT_STR = {(0, 0): "0/0", (0, 1): "0/1", (1, 0): "0/1", (1, 1): "1/1"}

    vcf     = VCF(vcf_path)
    samples = list(vcf.samples)
    extra_pos_meta = {}
    extra_gts      = defaultdict(dict)

    for v in vcf:
        if not v.ALT:
            continue
        key = (v.CHROM, v.POS)
        extra_pos_meta[key] = (v.REF, str(v.ALT[0]))
        for i, s in enumerate(samples):
            a1, a2 = v.genotypes[i][0], v.genotypes[i][1]
            extra_gts[key][s] = GT_STR.get((a1, a2), "./.")

    vcf.close()
    return samples, extra_pos_meta, extra_gts


def _chrom_sort_key(chrom):
    return _CHROM_ORDER.get(chrom, 99), chrom


def write_vcf(pos_meta, gts, all_labels, output_path):
    """Write a minimal multi-sample VCF with a single GT FORMAT field."""
    sorted_positions = sorted(
        pos_meta.keys(),
        key=lambda k: (_chrom_sort_key(k[0]), k[1])
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as fh:
        fh.write("##fileformat=VCFv4.1\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        cols += all_labels
        fh.write("\t".join(cols) + "\n")

        for key in sorted_positions:
            chrom, pos = key
            ref, alt   = pos_meta[key]
            row_gts    = [gts[key].get(lbl, "./.") for lbl in all_labels]
            fh.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
                + "\t".join(row_gts) + "\n"
            )

    return len(sorted_positions)


def main():
    args = parse_args()
    conn = sqlite3.connect(args.db)

    # ── Resolve sources ───────────────────────────────────────────────────────
    if args.cell_lines:
        print(f"Resolving {len(args.cell_lines)} cell line(s) from DB …")
        cl_to_run = resolve_run_ids(conn, args.cell_lines)
        if not cl_to_run:
            raise SystemExit("ERROR: no cell lines could be resolved from the DB.")
        run_id_to_label = {rid: cl for cl, rid in cl_to_run.items()}
        for cl, rid in cl_to_run.items():
            print(f"  {cl:<14} → {rid}")
    else:
        print(f"Using {len(args.run_ids)} explicit run_id(s) …")
        run_id_to_label = {}
        for rid in args.run_ids:
            row = conn.execute(
                "SELECT COALESCE(s.cell_line, s.name) "
                "FROM runs r JOIN samples s ON r.sample_id = s.sample_id "
                "WHERE r.run_id = ?",
                (rid,),
            ).fetchone()
            label = row[0] if row else rid
            run_id_to_label[rid] = label
            print(f"  {rid} → {label}")

    # ── Load DB genotypes ─────────────────────────────────────────────────────
    print(f"\nLoading genotypes from DB …")
    pos_meta, gts = load_db_genotypes(conn, run_id_to_label)
    conn.close()
    print(f"  {len(pos_meta):,} positions, {len(run_id_to_label)} DB donor(s)")

    all_labels = list(run_id_to_label.values())

    # ── Merge extra VCF ───────────────────────────────────────────────────────
    if args.extra_vcf:
        print(f"\nLoading extra donors from {args.extra_vcf} …")
        extra_samples, extra_pos_meta, extra_gts = load_extra_vcf(args.extra_vcf)
        print(f"  {len(extra_samples)} extra sample(s): {extra_samples}")

        new_positions = 0
        for key, sample_gts in extra_gts.items():
            if key not in pos_meta:
                pos_meta[key] = extra_pos_meta[key]
                new_positions += 1
            for s, gt in sample_gts.items():
                gts[key][s] = gt

        all_labels += extra_samples
        if new_positions:
            print(f"  {new_positions:,} new positions added from extra VCF")

    # ── Write output ──────────────────────────────────────────────────────────
    n = write_vcf(pos_meta, gts, all_labels, args.output)
    print(f"\nWrote {n:,} positions × {len(all_labels)} donor(s) → {args.output}")


if __name__ == "__main__":
    main()
