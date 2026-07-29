"""
Re-call genotypes at common SNP panel positions for runs that were originally
processed without the -T panel filter (genome-wide mpileup, variants-only).

For each affected run:
  1. bcftools call -m -T <panel> on the existing mpileup.bcf  → panel_raw.vcf
  2. bcftools filter (DP / QUAL thresholds)                   → panel_genotyped.vcf
  3. Delete old genotype_calls from SQLite, reload from new VCF

Run metadata (mapping rate, etc.) is unchanged. Grist is not updated.

Usage:
    python scripts/backfill_panel_genotypes.py \\
        --db results/variants.db \\
        --panel data/reference/genome1K.hg38.common_snps.vcf.gz \\
        [--run-ids SRR5071667_hg38 SRR5071672_hg38 ...]  # default: all missing 0/0

    # dry-run to preview without writing:
    python scripts/backfill_panel_genotypes.py --dry-run
"""

import argparse
import sqlite3
import subprocess
import sys
from pathlib import Path

from cyvcf2 import VCF


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Backfill 0/0 panel genotypes for runs missing ref calls"
    )
    p.add_argument("--db",        default="results/variants.db")
    p.add_argument("--panel",     default="data/reference/genome1K.hg38.common_snps.vcf.gz",
                   help="common SNP panel VCF (must be tabix-indexed)")
    p.add_argument("--results",   default="results",
                   help="results root directory")
    p.add_argument("--run-ids",   nargs="+", default=None,
                   help="specific run_ids to process (default: all missing 0/0 in DB)")
    p.add_argument("--min-depth", type=int,   default=10)
    p.add_argument("--min-qual",  type=float, default=30.0)
    p.add_argument("--dry-run",   action="store_true",
                   help="print commands without running them")
    return p.parse_args()


# ── DB helpers ────────────────────────────────────────────────────────────────

def find_runs_missing_ref(db_path):
    """Return run_ids that have zero 0/0 genotype calls."""
    conn = sqlite3.connect(db_path)
    all_runs = [r[0] for r in conn.execute("SELECT run_id FROM runs").fetchall()]
    missing = []
    for rid in all_runs:
        n = conn.execute(
            "SELECT COUNT(*) FROM genotype_calls WHERE run_id=? AND genotype='0/0'",
            (rid,)
        ).fetchone()[0]
        if n == 0:
            missing.append(rid)
    conn.close()
    return missing


def load_variants_sqlite(conn, vcf_path, run_id):
    """Delete old calls for run_id and insert from vcf_path. Returns row count."""
    conn.execute("DELETE FROM genotype_calls WHERE run_id=?", (run_id,))

    vcf = VCF(str(vcf_path))
    GT_MAP = {"0/0": "0/0", "0/1": "0/1", "1/0": "0/1", "1/1": "1/1"}
    rows = 0

    for v in vcf:
        chrom = v.CHROM
        pos   = v.POS
        ref   = v.REF
        alt   = str(v.ALT[0]) if v.ALT else "."
        qual  = float(v.QUAL) if v.QUAL else None
        filt  = v.FILTER or "PASS"
        dp    = int(v.INFO.get("DP") or 0)
        vtype = "SNP" if len(ref) == len(alt) == 1 else "INDEL"

        gt_str = "."
        if v.genotypes:
            raw = v.genotypes[0]
            alleles = [str(a) for a in raw[:-1]]
            gt_str  = GT_MAP.get("/".join(alleles), "/".join(alleles))

        rd, ad, af = 0, 0, None
        try:
            ad_field = v.format("AD")
            if ad_field is not None:
                rd = int(ad_field[0][0])
                ad = int(ad_field[0][1])
                af = round(ad / dp, 4) if dp > 0 else None
        except Exception:
            pass

        conn.execute("""
            INSERT OR IGNORE INTO variants
                (chromosome, position, ref_allele, alt_allele, variant_type)
            VALUES (?, ?, ?, ?, ?)
        """, (chrom, pos, ref, alt, vtype))

        vid = conn.execute("""
            SELECT variant_id FROM variants
            WHERE chromosome=? AND position=? AND ref_allele=? AND alt_allele=?
        """, (chrom, pos, ref, alt)).fetchone()[0]

        conn.execute("""
            INSERT INTO genotype_calls
                (run_id, variant_id, genotype, quality, depth,
                 ref_depth, alt_depth, allele_freq, filter_status)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (run_id, vid, gt_str, qual, dp, rd, ad, af, filt))
        rows += 1

    conn.commit()
    return rows


# ── Variant-calling helpers ───────────────────────────────────────────────────

def run_cmd(cmd, log_path=None, dry_run=False):
    cmd_str = " ".join(str(c) for c in cmd)
    print(f"    $ {cmd_str}")
    if dry_run:
        return
    result = subprocess.run(cmd, capture_output=True, text=True)
    if log_path:
        Path(log_path).write_text(result.stderr)
    if result.returncode != 0:
        print(result.stderr[-2000:], file=sys.stderr)
        raise RuntimeError(f"Command failed (exit {result.returncode})")


def bcf_to_run_id(run_id):
    """Derive sample name from run_id (strip _hg38 suffix)."""
    return run_id.replace("_hg38", "")


def panel_vcf_paths(results_root, run_id):
    sample   = bcf_to_run_id(run_id)
    vcf_dir  = Path(results_root) / run_id / "vcf"
    raw      = vcf_dir / f"{sample}.panel_raw.vcf"
    filtered = vcf_dir / f"{sample}.panel_genotyped.vcf"
    bcf      = vcf_dir / f"{sample}.mpileup.bcf"
    return bcf, raw, filtered


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if args.run_ids:
        run_ids = args.run_ids
    else:
        print(f"Scanning {args.db} for runs without 0/0 calls...")
        run_ids = find_runs_missing_ref(args.db)
        print(f"  Found {len(run_ids)} affected run(s): {run_ids}")

    if not run_ids:
        print("Nothing to do.")
        return

    conn = None if args.dry_run else sqlite3.connect(args.db)

    for run_id in run_ids:
        print(f"\n{'='*60}")
        print(f"Run: {run_id}")
        bcf, raw_vcf, filt_vcf = panel_vcf_paths(args.results, run_id)

        if not bcf.exists():
            print(f"  SKIP: mpileup BCF not found at {bcf}")
            continue

        # Step 1: bcftools call restricted to SNP panel
        print(f"  Step 1: bcftools call -> {raw_vcf}")
        run_cmd([
            "bcftools", "call",
            "-m",                        # multiallelic caller
            "-T", args.panel,            # restrict to panel positions (includes 0/0)
            "--output-type", "v",
            "-o", raw_vcf,
            bcf,
        ], dry_run=args.dry_run)

        # Step 2: depth / quality filter
        print(f"  Step 2: bcftools filter -> {filt_vcf}")
        run_cmd([
            "bcftools", "filter",
            "-e", f'DP < {args.min_depth} || (GT!="0/0" && QUAL < {args.min_qual})',
            "-o", filt_vcf,
            raw_vcf,
        ], dry_run=args.dry_run)

        if args.dry_run:
            print(f"  Step 3: would reload DB ({args.db}) for {run_id}")
            continue

        # Step 3: reload genotype calls in SQLite
        print(f"  Step 3: reloading {run_id} in {args.db}...")
        n = load_variants_sqlite(conn, filt_vcf, run_id)
        n_ref = conn.execute(
            "SELECT COUNT(*) FROM genotype_calls WHERE run_id=? AND genotype='0/0'",
            (run_id,)
        ).fetchone()[0]
        n_alt = conn.execute(
            "SELECT COUNT(*) FROM genotype_calls WHERE run_id=? AND genotype!='0/0'",
            (run_id,)
        ).fetchone()[0]
        print(f"  Done: {n:,} total rows  ({n_ref:,} 0/0  {n_alt:,} ALT)")

    if conn:
        conn.close()
    print("\nAll done.")


if __name__ == "__main__":
    main()
