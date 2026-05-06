# scripts/load_to_db.py
"""
Loads pipeline outputs into:
  - Grist (run metadata + QC summary)
  - SQLite (full variant data)

Usage:
    python scripts/load_to_db.py \
        --run-id SRR1258218_chr22 \
        --sample SRR1258218 \
        --vcf results/SRR1258218_chr22/vcf/SRR1258218.filtered.vcf \
        --flagstat results/SRR1258218_chr22/bam/SRR1258218.flagstat.txt \
        --db variants.db
"""

import os
import re
import sqlite3
import argparse
import subprocess
from datetime import date
from dotenv import load_dotenv
from grist_api import GristDocAPI
from cyvcf2 import VCF

load_dotenv()


# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Load pipeline outputs to Grist + SQLite")
    p.add_argument("--run-id",   required=True)
    p.add_argument("--sample",   required=True)
    p.add_argument("--vcf",      required=True, help="path to filtered VCF")
    p.add_argument("--flagstat", required=True, help="path to samtools flagstat output")
    p.add_argument("--db",       default="results/variants.db", help="SQLite database path")
    p.add_argument("--reference", default="data/reference/genome.fa")
    p.add_argument("--sequencing", default="single")
    p.add_argument("--star-log",  default=None, help="path to STAR Log.final.out")
    return p.parse_args()


# ─────────────────────────────────────────────
# Parse flagstat
# ─────────────────────────────────────────────
def parse_flagstat(path):
    """Extract total reads, mapped reads, and mapping rate from flagstat output."""
    total, mapped, rate = 0, 0, 0.0
    with open(path) as f:
        for line in f:
            if "in total" in line:
                total = int(line.split()[0])
            if "mapped (" in line:
                mapped = int(line.split()[0])
                match = re.search(r'([\d.]+)%', line)
                if match:
                    rate = float(match.group(1))
    return total, mapped, rate


# ─────────────────────────────────────────────
# Parse STAR Log.final.out
# ─────────────────────────────────────────────
def parse_star_log(path):
    """Extract annotated splice % and pct_too_short from STAR's final log."""
    if not path or not os.path.exists(path):
        return None, None
    total_splices, annotated_splices, pct_too_short = 0, 0, 0.0
    with open(path) as f:
        for line in f:
            if "Number of splices: Total" in line:
                total_splices = int(line.split("|")[1].strip())
            elif "Number of splices: Annotated (sjdb)" in line:
                annotated_splices = int(line.split("|")[1].strip())
            elif "% of reads unmapped: too short" in line:
                pct_too_short = float(line.split("|")[1].strip().replace("%", ""))
    annotated_splice_pct = (
        round(annotated_splices / total_splices * 100, 2) if total_splices > 0 else None
    )
    return annotated_splice_pct, pct_too_short


# ─────────────────────────────────────────────
# Get tool versions
# ─────────────────────────────────────────────
def get_version(cmd):
    """Run a command and return first line of output."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True
        )
        # STAR prints to stdout, bcftools to stderr
        output = result.stdout or result.stderr
        return output.strip().splitlines()[0]
    except Exception:
        return "unknown"


# ─────────────────────────────────────────────
# SQLite setup
# ─────────────────────────────────────────────
def init_sqlite(db_path):
    """Create tables if they don't exist."""
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS samples (
            sample_id    INTEGER PRIMARY KEY AUTOINCREMENT,
            name         TEXT UNIQUE,
            date_added   TEXT
        );

        CREATE TABLE IF NOT EXISTS runs (
            run_id       TEXT PRIMARY KEY,
            sample_id    INTEGER REFERENCES samples(sample_id),
            reference    TEXT,
            sequencing   TEXT,
            run_date     TEXT,
            mapping_rate REAL,
            total_raw    INTEGER,
            total_filt   INTEGER
        );

        CREATE TABLE IF NOT EXISTS variants (
            variant_id   INTEGER PRIMARY KEY AUTOINCREMENT,
            chromosome   TEXT,
            position     INTEGER,
            ref_allele   TEXT,
            alt_allele   TEXT,
            variant_type TEXT,
            UNIQUE(chromosome, position, ref_allele, alt_allele)
        );

        CREATE TABLE IF NOT EXISTS genotype_calls (
            call_id      INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id       TEXT REFERENCES runs(run_id),
            variant_id   INTEGER REFERENCES variants(variant_id),
            genotype     TEXT,
            quality      REAL,
            depth        INTEGER,
            ref_depth    INTEGER,
            alt_depth    INTEGER,
            allele_freq  REAL,
            filter_status TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos
            ON variants(chromosome, position);
        CREATE INDEX IF NOT EXISTS idx_calls_run
            ON genotype_calls(run_id);
    """)
    conn.commit()
    return conn


# ─────────────────────────────────────────────
# Load variants into SQLite
# ─────────────────────────────────────────────
def load_variants_sqlite(conn, vcf_path, run_id):
    """Parse VCF and insert variants + genotype calls into SQLite."""
    vcf = VCF(vcf_path)
    rows_inserted = 0

    for variant in vcf:
        chrom    = variant.CHROM
        pos      = variant.POS
        ref      = variant.REF
        alt      = str(variant.ALT[0]) if variant.ALT else "."
        qual     = float(variant.QUAL) if variant.QUAL else None
        filt     = variant.FILTER or "PASS"
        dp       = variant.INFO.get("DP") or 0
        vtype    = "SNP" if len(ref) == len(alt) == 1 else "INDEL"

        # Parse genotype fields
        gt_str   = "."
        rd, ad   = 0, 0
        af       = None

        if variant.genotypes:
            gt_raw = variant.genotypes[0]
            # gt_raw is a list like [0, 1, True] meaning allele1, allele2, phased
            alleles = [str(a) for a in gt_raw[:-1]]
            gt_str  = "/".join(alleles)

        # AD field: ref depth, alt depth
        try:
            ad_field = variant.format("AD")
            if ad_field is not None:
                rd = int(ad_field[0][0])
                ad = int(ad_field[0][1])
                af = round(ad / dp, 4) if dp > 0 else None
        except Exception:
            pass

        # Insert or ignore the variant position
        conn.execute("""
            INSERT OR IGNORE INTO variants
                (chromosome, position, ref_allele, alt_allele, variant_type)
            VALUES (?, ?, ?, ?, ?)
        """, (chrom, pos, ref, alt, vtype))

        variant_id = conn.execute("""
            SELECT variant_id FROM variants
            WHERE chromosome=? AND position=? AND ref_allele=? AND alt_allele=?
        """, (chrom, pos, ref, alt)).fetchone()[0]

        conn.execute("""
            INSERT INTO genotype_calls
                (run_id, variant_id, genotype, quality, depth,
                 ref_depth, alt_depth, allele_freq, filter_status)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (run_id, variant_id, gt_str, qual, dp, rd, ad, af, filt))

        rows_inserted += 1

    conn.commit()
    return rows_inserted


# ─────────────────────────────────────────────
# Push to Grist
# ─────────────────────────────────────────────
def grist_delete_existing(api, table, run_id):
    """Delete any rows for this run_id so re-runs don't create duplicates."""
    records = api.fetch_table(table)
    ids = [r.id for r in records if getattr(r, "run_id", None) == run_id]
    if ids:
        api.delete_records(table, ids)
        print(f"  Removed {len(ids)} existing row(s) from {table}")


def push_to_grist(args, total_reads, mapped_reads,
                  mapping_rate, raw_variants, filt_variants,
                  star_ver, bcftools_ver, pct_too_short,
                  mean_depth, annotated_splice_pct):
    """Push run summary, QC metrics, and sample record to Grist."""
    api = GristDocAPI(
        os.environ["GRIST_DOC_ID"],
        server=os.environ["GRIST_SERVER"],
        api_key=os.environ["GRIST_API_KEY"]
    )

    grist_delete_existing(api, "Runs", args.run_id)
    api.add_records("Runs", [{
        "sample_name":             args.sample,
        "run_id":                  args.run_id,
        "pipeline_version":        "0.1.0",
        "reference_genome":        args.reference,
        "sequencing_mode":         args.sequencing,
        "star_version":            star_ver,
        "bcftools_version":        bcftools_ver,
        "run_date":                str(date.today()),
        "mapping_rate":            mapping_rate,
        "total_raw_variants":      raw_variants,
        "total_filtered_variants": filt_variants,
        "vcf_file_path":           args.vcf,
    }])
    print(f"  ✓ Run summary pushed to Grist Runs table")

    grist_delete_existing(api, "QC_Summary", args.run_id)
    api.add_records("QC_Summary", [{
        "run_id":               args.run_id,
        "mapped_reads":         mapped_reads,
        "total_reads":          total_reads,
        "mean_depth":           mean_depth,
        "pct_too_short":        pct_too_short,
        "annotated_splice_pct": annotated_splice_pct,
        "notes":                f"{filt_variants} filtered variants",
    }])
    print(f"  ✓ QC summary pushed to Grist QC_Summary table")

    # Add sample record if not already present
    existing = api.fetch_table("Samples")
    if args.sample not in [r.name for r in existing]:
        api.add_records("Samples", [{
            "name":       args.sample,
            "date_added": str(date.today()),
        }])
        print(f"  ✓ Sample record added to Grist Samples table")
    else:
        print(f"  Sample {args.sample} already in Grist Samples table")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    args = parse_args()

    print(f"\nLoading results for run: {args.run_id}")

    # 1. Parse flagstat
    print("  Parsing flagstat...")
    total_reads, mapped_reads, mapping_rate = parse_flagstat(args.flagstat)
    print(f"  Total reads: {total_reads:,}  Mapped: {mapped_reads:,}  Rate: {mapping_rate:.2f}%")

    # 2. Count variants in raw + filtered VCF
    print("  Counting variants...")
    raw_vcf_path = args.vcf.replace(".filtered.vcf", ".raw.vcf")
    raw_variants  = sum(1 for v in VCF(raw_vcf_path))
    filt_variants = sum(1 for v in VCF(args.vcf))
    print(f"  Raw variants: {raw_variants:,}  Filtered: {filt_variants:,}")

    # 3. Get tool versions
    star_ver     = get_version(["STAR", "--version"])
    bcftools_ver = get_version(["bcftools", "--version"])
    print(f"  STAR: {star_ver}  bcftools: {bcftools_ver}")

    # 4. Load into SQLite
    print(f"  Loading variants into SQLite ({args.db})...")
    conn = init_sqlite(args.db)

    # Insert sample if not exists
    conn.execute(
        "INSERT OR IGNORE INTO samples (name, date_added) VALUES (?, ?)",
        (args.sample, str(date.today()))
    )
    sample_id = conn.execute(
        "SELECT sample_id FROM samples WHERE name=?", (args.sample,)
    ).fetchone()[0]

    # Insert run
    conn.execute("""
        INSERT OR REPLACE INTO runs
            (run_id, sample_id, reference, sequencing,
             run_date, mapping_rate, total_raw, total_filt)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, (args.run_id, sample_id, args.reference, args.sequencing,
          str(date.today()), mapping_rate, raw_variants, filt_variants))
    conn.commit()

    rows = load_variants_sqlite(conn, args.vcf, args.run_id)
    print(f"  ✓ {rows:,} variant calls loaded into SQLite")
    conn.close()

    # 5. Parse STAR log for additional QC metrics
    print("  Parsing STAR log...")
    annotated_splice_pct, pct_too_short = parse_star_log(args.star_log)
    print(f"  Annotated splice %: {annotated_splice_pct}  Too short %: {pct_too_short}")

    # 6. Calculate mean depth from SQLite
    conn = sqlite3.connect(args.db)
    mean_depth = conn.execute(
        "SELECT ROUND(AVG(depth), 1) FROM genotype_calls WHERE run_id=?", (args.run_id,)
    ).fetchone()[0]
    conn.close()
    print(f"  Mean depth: {mean_depth}x")

    # 7. Push to Grist
    print("  Pushing to Grist...")
    push_to_grist(
        args, total_reads, mapped_reads, mapping_rate,
        raw_variants, filt_variants, star_ver, bcftools_ver,
        pct_too_short, mean_depth, annotated_splice_pct,
    )

    print(f"\nDone. Run {args.run_id} loaded successfully.\n")


if __name__ == "__main__":
    main()