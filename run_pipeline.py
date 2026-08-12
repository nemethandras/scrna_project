#!/usr/bin/env python3
"""Command-line wrapper for the scRNA-seq variant calling pipeline."""

import argparse
import csv
import os
import subprocess
import sys

# Default paths to reference data bundled with the workflow
DEFAULT_REFERENCE  = "data/reference/genome.fa"
DEFAULT_GTF        = "data/reference/genes.gtf"
DEFAULT_STAR_INDEX = "data/reference/star_index_hg38_oh99"
DEFAULT_SNPS_VCF   = "data/reference/genome1K.hg38.common_snps.vcf.gz"
DEFAULT_DB         = "results/variants.db"
DEFAULT_CELL_LINE_MAP = "data/cell_line_map.csv"


def lookup_cell_line(map_path, sample_id):
    """Return cell line name for a sample_id from a CSV/TSV map, or None."""
    if not map_path or not os.path.exists(map_path):
        return None
    try:
        with open(map_path) as f:
            dialect = csv.Sniffer().sniff(f.read(2048))
            f.seek(0)
            for row in csv.DictReader(f, dialect=dialect):
                if row.get("Run") == sample_id or row.get("sample_id") == sample_id:
                    return row.get("cell_line") or row.get("cell_line_name") or None
    except Exception:
        pass
    return None


# ─────────────────────────────────────────────
# Bulk mode helpers
# ─────────────────────────────────────────────

def build_snakemake_cmd(run_id, sample, args, fastq_dir=None):
    fastq_dir = fastq_dir or f"data/{sample}"
    cell_line = lookup_cell_line(
        getattr(args, "cell_line_map", None) or DEFAULT_CELL_LINE_MAP, sample
    )
    config_overrides = [
        f"run_id={run_id}",
        f"samples=[{sample}]",
        f"fastq_dir={fastq_dir}",
        f"sequencing={args.sequencing}",
        f"reference_genome={args.reference}",
        f"annotation_gtf={args.gtf}",
        f"star_index_dir={args.star_index}",
        f"sjdb_overhang={args.sjdb_overhang}",
        f"genome_sa_index_nbases={args.genome_sa_index_nbases}",
    ]
    if cell_line:
        config_overrides.append(f"cell_line={cell_line}")
    if args.db:
        config_overrides.append(f"db_path={args.db}")
    if args.no_db:
        config_overrides.append("no_db=True")

    cmd = [
        "snakemake",
        "--snakefile", "workflow/Snakefile",
        "--cores", str(args.cores),
        "--config", *config_overrides,
    ]
    if args.dry_run:
        cmd.append("-n")
    if args.rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.force:
        cmd.append("--forceall")
    return cmd


def run_single(run_id, sample, args, fastq_dir=None):
    cmd = build_snakemake_cmd(run_id, sample, args, fastq_dir)

    if args.foreground or args.dry_run:
        result = subprocess.run(cmd)
        sys.exit(result.returncode)

    log_path = f"logs/nohup_{sample}.log"
    os.makedirs("logs", exist_ok=True)
    with open(log_path, "w") as log_file:
        proc = subprocess.Popen(
            cmd,
            stdout=log_file,
            stderr=log_file,
            start_new_session=True,
        )
    print(f"Pipeline running in background (PID {proc.pid})")
    print(f"Follow progress: tail -f {log_path}")


def _notify_cmd(email, subject, body):
    """Return a shell snippet that sends a one-line email notification."""
    safe_subject = subject.replace('"', '\\"')
    safe_body    = body.replace('"', '\\"')
    return f'echo "{safe_body}" | mail -s "{safe_subject}" {email}'


def run_batch(samples, suffix, args):
    """Run samples sequentially in a single detached process."""
    os.makedirs("logs", exist_ok=True)
    log_path = "logs/batch.log"
    run_ids   = ", ".join(f"{s}{suffix}" for s in samples)

    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f"exec >> {log_path} 2>&1",
        f'echo "Batch started: $(date)"',
    ]

    if args.notify_email:
        script_lines += [
            f'_on_exit() {{',
            f'  local code=$?',
            f'  if [ $code -ne 0 ]; then',
            f'    ' + _notify_cmd(args.notify_email,
                                  f"Pipeline FAILED [{run_ids}]",
                                  f"A sample failed. Check {log_path} for details."),
            f'  fi',
            f'}}',
            f'trap _on_exit EXIT',
        ]

    for sample in samples:
        run_id = f"{sample}{suffix}"
        cmd = build_snakemake_cmd(run_id, sample, args)
        cmd_str = " ".join(cmd)
        script_lines += [
            f'echo ""',
            f'echo "--- {sample} ({run_id}) started: $(date) ---"',
            cmd_str,
            f'echo "--- {sample} ({run_id}) finished: $(date) ---"',
        ]

    script_lines.append('echo "Batch finished: $(date)"')
    if args.notify_email:
        script_lines.append(
            _notify_cmd(args.notify_email,
                        f"Pipeline done [{run_ids}]",
                        f"All {len(samples)} samples finished. Check {log_path}.")
        )

    script = "\n".join(script_lines) + "\n"

    with open(log_path, "w") as log_file:
        proc = subprocess.Popen(
            ["bash", "-c", script],
            stdout=log_file,
            stderr=log_file,
            start_new_session=True,
        )
    print(f"Batch of {len(samples)} samples running in background (PID {proc.pid})")
    print(f"Samples: {', '.join(samples)}")
    print(f"Follow progress: tail -f {log_path}")


# ─────────────────────────────────────────────
# scRNA mode
# ─────────────────────────────────────────────

def run_scrna(args):
    """Run cellsnp-lite (if needed) then demux.py, detached."""
    demux_run_id = args.demux_run_id
    cellsnp_dir  = args.cellsnp_dir or f"results/demux/{demux_run_id}/cellsnp"
    output       = f"results/demux/{demux_run_id}/assignments.tsv"
    db           = args.db or DEFAULT_DB
    log_path     = f"logs/demux/{demux_run_id}.log"

    os.makedirs(f"results/demux/{demux_run_id}", exist_ok=True)
    os.makedirs("logs/demux", exist_ok=True)

    run_ids_str      = " ".join(args.run_ids)
    load_db_flag     = "--load-db" if args.load_db else ""

    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f'echo "scRNA demux {demux_run_id} started: $(date)"',
    ]

    # cellsnp-lite — skip if the user already has a cellsnp dir
    if args.bam:
        script_lines += [
            f'echo "--- cellsnp-lite ---"',
            f'mkdir -p {cellsnp_dir}',
            (
                f'cellsnp-lite'
                f' -s {args.bam}'
                f' -b {args.barcodes}'
                f' -O {cellsnp_dir}'
                f' -R {args.snps_vcf}'
                f' -p {args.cores}'
                f' --minMAF 0 --minCOUNT 20 --gzip'
            ),
        ]
    else:
        script_lines.append(
            f'echo "Using existing cellsnp dir: {cellsnp_dir}"'
        )

    vireo_dir     = f"results/demux/{demux_run_id}/vireo"
    donor_matches = f"results/demux/{demux_run_id}/donor_matches.tsv"
    final_out     = f"results/demux/{demux_run_id}/final_assignments.tsv"

    # ── Determine Vireo mode ──────────────────────────────────────────────────
    guided = bool(args.known_lines or args.donor_vcf)
    donor_ref_vcf = args.donor_vcf or f"results/demux/{demux_run_id}/donor_ref.vcf"

    if guided:
        vireo_d_flag = f"-d {donor_ref_vcf}"
        n_flag       = ""   # donor count inferred from VCF columns
    else:
        vireo_d_flag = ""
        n_flag       = f"-N {args.n_donors}" if args.n_donors else ""

    script_lines += [
        f'echo "--- demux (binomial scorer) ---"',
        (
            f'python scripts/demux.py'
            f' --cellsnp-dir {cellsnp_dir}'
            f' --db {db}'
            f' --run-ids {run_ids_str}'
            f' --output {output}'
            f' --demux-run-id {demux_run_id}'
            f' --min-depth {args.min_depth}'
            f' --min-positions {args.min_positions}'
            f' --doublet-gap {args.doublet_gap}'
            f' --no-match-threshold {args.no_match_threshold}'
            f' {load_db_flag}'
        ),
    ]

    # ── Build donor reference VCF when --known-lines is given ─────────────────
    if args.known_lines and not args.donor_vcf:
        known_str = " ".join(args.known_lines)
        script_lines += [
            f'echo "--- building donor reference VCF from DB ---"',
            (
                f'python scripts/make_donor_vcf.py'
                f' --cell-lines {known_str}'
                f' --db {db}'
                f' --output {donor_ref_vcf}'
            ),
        ]
    elif args.known_lines and args.donor_vcf:
        # --known-lines + --donor-vcf: DB lines merged with extra lines from file.
        # The --donor-vcf is the *input* extra-lines VCF; output always goes to the
        # auto-generated path.  Guard against self-overwrite when the user points
        # --donor-vcf at the auto-generated path (would duplicate every donor).
        auto_vcf = f"results/demux/{demux_run_id}/donor_ref.vcf"
        if os.path.abspath(args.donor_vcf) == os.path.abspath(auto_vcf):
            raise SystemExit(
                f"ERROR: --donor-vcf points at '{auto_vcf}', which is the auto-generated\n"
                f"  output path for this demux run.  This would read and write the same file.\n"
                f"  To reuse an existing VCF, omit --known-lines and pass --donor-vcf only.\n"
                f"  To rebuild from the DB, omit --donor-vcf and pass --known-lines only."
            )
        known_str = " ".join(args.known_lines)
        donor_ref_vcf = auto_vcf
        vireo_d_flag  = f"-d {donor_ref_vcf}"
        script_lines += [
            f'echo "--- building donor reference VCF (DB + extra VCF) ---"',
            (
                f'python scripts/make_donor_vcf.py'
                f' --cell-lines {known_str}'
                f' --extra-vcf {args.donor_vcf}'
                f' --db {db}'
                f' --output {donor_ref_vcf}'
            ),
        ]

    # ── Vireo ─────────────────────────────────────────────────────────────────
    mode_label = "reference-guided" if guided else "genotype-free"
    script_lines += [
        f'echo "--- vireo ({mode_label} clustering) ---"',
        f'mkdir -p {vireo_dir}',
        f'vireo -c {cellsnp_dir} {n_flag} {vireo_d_flag}'
        + (' --genoTag GT' if guided else '')
        + f' -o {vireo_dir} --randSeed 42',
    ]

    # ── Post-Vireo matching (genotype-free only) ──────────────────────────────
    if guided:
        script_lines += [
            f'echo "--- merging results (reference-guided) ---"',
            (
                f'python scripts/merge_demux.py'
                f' --vireo-dir {vireo_dir}'
                f' --scorer {output}'
                f' --output {final_out}'
                f' --guided'
            ),
        ]
    else:
        script_lines += [
            f'echo "--- matching Vireo donors to DB ---"',
            (
                f'python scripts/match_vireo.py'
                f' --vireo-dir {vireo_dir}'
                f' --db {db}'
                f' --run-ids {run_ids_str}'
                f' --output {donor_matches}'
                f' --min-concordance {args.min_concordance}'
                f' --min-gap {args.min_gap}'
                + (' --unique' if args.unique else '')
            ),
            f'echo "--- merging results ---"',
            (
                f'python scripts/merge_demux.py'
                f' --vireo-dir {vireo_dir}'
                f' --matches {donor_matches}'
                f' --scorer {output}'
                f' --output {final_out}'
                f' --min-concordance {args.min_concordance}'
            ),
        ]

    if args.add_unmatched and not guided:
        script_lines += [
            f'echo "--- adding unmatched donors to DB ---"',
            (
                f'python scripts/add_unmatched_donors.py'
                f' --vireo-dir {vireo_dir}'
                f' --matches {donor_matches}'
                f' --db {db}'
                f' --demux-run-id {demux_run_id}'
            ),
        ]

    script_lines += [
        f'echo "Done: $(date)"',
        f'echo "Binomial assignments : {output}"',
        f'echo "Vireo assignments    : {vireo_dir}/donor_ids.tsv"',
        f'echo "Final merged output  : {final_out}"',
    ]
    if not guided:
        script_lines.insert(-1, f'echo "Vireo→DB matches     : {donor_matches}"')

    script = "\n".join(script_lines) + "\n"

    if args.dry_run:
        print(script)
        return

    if args.foreground:
        result = subprocess.run(["bash", "-c", script])
        sys.exit(result.returncode)

    with open(log_path, "w") as log_file:
        proc = subprocess.Popen(
            ["bash", "-c", script],
            stdout=log_file,
            stderr=log_file,
            start_new_session=True,
        )
    print(f"scRNA demux running in background (PID {proc.pid})")
    print(f"Follow progress: tail -f {log_path}")
    print(f"Output: {output}")


# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="Variant calling and scRNA demultiplexing pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # bulk — single sample
  source .env && python run_pipeline.py --mode bulk --run-id SRR5071686_hg38 --sample SRR5071686

  # bulk — batch
  source .env && python run_pipeline.py --mode bulk --samples SRR5071662 SRR5071667 --run-id-suffix _hg38

  # scrna — from CellRanger BAM (runs cellsnp-lite then demux)
  source .env && python run_pipeline.py --mode scrna \\
      --bam data/Pool_ctr/possorted_genome_bam.bam \\
      --barcodes data/Pool_ctr/barcodes.tsv.gz \\
      --demux-run-id pool_ctr_001

  # scrna — from existing cellsnp output (skip cellsnp-lite)
  source .env && python run_pipeline.py --mode scrna \\
      --cellsnp-dir data/Pool_ctr/cellsnp \\
      --demux-run-id pool_ctr_001

  # dry run
  python run_pipeline.py --mode bulk --run-id SRR5071686_hg38 --sample SRR5071686 --dry-run
        """,
    )

    parser.add_argument(
        "--mode", choices=["bulk", "scrna"], default="bulk",
        help="pipeline mode: bulk (FASTQ→genotype→DB) or scrna (BAM→cellsnp→demux) (default: bulk)",
    )

    # ── Bulk args ──────────────────────────────────────────────────────────
    bulk = parser.add_argument_group("bulk mode")
    sample_input = bulk.add_mutually_exclusive_group()
    sample_input.add_argument(
        "--sample", metavar="ID",
        help="single sample ID — use with --run-id",
    )
    sample_input.add_argument(
        "--samples", nargs="+", metavar="ID",
        help="one or more sample IDs for a batch run — use with --run-id-suffix",
    )
    bulk.add_argument(
        "--run-id", metavar="ID",
        help="run label for a single sample (required with --sample)",
    )
    bulk.add_argument(
        "--run-id-suffix", metavar="SUFFIX",
        help="suffix appended to each sample ID to form the run ID (required with --samples)",
    )
    bulk.add_argument(
        "--fastq-dir", default=None, metavar="PATH",
        help="FASTQ directory for single-sample runs (default: data/<sample>)",
    )
    bulk.add_argument(
        "--sequencing", choices=["paired", "single"], default="single",
        help="sequencing layout (default: single)",
    )
    bulk.add_argument(
        "--no-db", action="store_true",
        help="skip the load_to_database step",
    )
    bulk.add_argument(
        "--cell-line-map", metavar="PATH", default=None,
        help=f"CSV/TSV mapping sample_id → cell_line for automatic DB labelling "
             f"(default: {DEFAULT_CELL_LINE_MAP} if it exists)",
    )

    # ── scRNA args ─────────────────────────────────────────────────────────
    scrna = parser.add_argument_group("scrna mode")
    scrna.add_argument(
        "--bam", metavar="PATH",
        help="BAM file from CellRanger or STARsolo (triggers cellsnp-lite)",
    )
    scrna.add_argument(
        "--barcodes", metavar="PATH",
        help="barcodes.tsv.gz from CellRanger or STARsolo (required with --bam)",
    )
    scrna.add_argument(
        "--cellsnp-dir", metavar="PATH",
        help="existing cellsnp-lite output directory (skips cellsnp-lite step)",
    )
    scrna.add_argument(
        "--demux-run-id", metavar="ID",
        help="unique label for this demux run (required for --mode scrna)",
    )
    scrna.add_argument(
        "--run-ids", nargs="+", default=["all"], metavar="ID",
        help="run_ids in the DB to use as references, or 'all' (default: all)",
    )
    scrna.add_argument(
        "--load-db", action="store_true",
        help="write cell assignments back into the SQLite database",
    )
    scrna.add_argument(
        "--min-depth", type=int, default=1, metavar="N",
        help="min read depth at a position in a cell (default: 1 for scRNA)",
    )
    scrna.add_argument(
        "--min-positions", type=int, default=200, metavar="N",
        help="min covered positions to attempt assignment (default: 200)",
    )
    scrna.add_argument(
        "--doublet-gap", type=float, default=2.0, metavar="F",
        help="min log-likelihood gap between top two donors for a singlet call (default: 2.0)",
    )
    scrna.add_argument(
        "--no-match-threshold", type=float, default=-0.5, metavar="F",
        help="mean LL per position below which the cell gets no_match (default: -0.5)",
    )
    scrna.add_argument(
        "--snps-vcf", metavar="PATH", default=DEFAULT_SNPS_VCF,
        help=f"common SNP panel VCF for cellsnp-lite (default: {DEFAULT_SNPS_VCF})",
    )
    scrna.add_argument(
        "--n-donors", type=int, default=None, metavar="N",
        help="number of donors for Vireo (default: auto-detect)",
    )
    scrna.add_argument(
        "--min-concordance", type=float, default=0.80, metavar="F",
        help="min genotype concordance to accept a Vireo donor→DB match (default: 0.80)",
    )
    scrna.add_argument(
        "--min-gap", type=float, default=0.10, metavar="F",
        help="best match must beat second-best by at least this (default: 0.10)",
    )
    scrna.add_argument(
        "--unique", action="store_true",
        help="enforce 1:1 donor↔reference matching via Hungarian algorithm (recommended when --run-ids lists exactly the pool members)",
    )
    scrna.add_argument(
        "--add-unmatched", action="store_true",
        help="add Vireo donors with no DB match to the DB as new unlabelled samples (run add_unmatched_donors.py as final step)",
    )
    scrna.add_argument(
        "--known-lines", nargs="+", metavar="NAME", default=None,
        help="cell line names known to be in the pool (enables reference-guided Vireo); "
             "genotypes are fetched from the DB automatically (most recent run per line)",
    )
    scrna.add_argument(
        "--donor-vcf", metavar="PATH", default=None,
        help="pre-built donor genotype VCF for reference-guided Vireo; "
             "use instead of --known-lines when you already have the file, or combine with "
             "--known-lines to supply extra lines not yet in the DB",
    )

    # ── Common options ─────────────────────────────────────────────────────
    parser.add_argument(
        "--cores", type=int, default=16,
        help="number of CPU cores (default: 16)",
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true",
        help="show what would run without executing anything",
    )
    parser.add_argument(
        "--rerun-incomplete", action="store_true",
        help="re-run jobs with incomplete output files from a previous failed run",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="force re-run of all steps even if outputs already exist",
    )
    parser.add_argument(
        "--foreground", action="store_true",
        help="run in the foreground instead of detaching",
    )
    parser.add_argument(
        "--notify-email", metavar="ADDR",
        help="send an email when the batch finishes or fails (bulk mode only)",
    )
    parser.add_argument(
        "--db", metavar="PATH", default=None,
        help=f"SQLite database file (default: {DEFAULT_DB})",
    )

    # ── Reference overrides (bulk) ─────────────────────────────────────────
    ref = parser.add_argument_group(
        "reference overrides (bulk mode)",
        "optional — defaults to the hg38 reference bundled in data/reference/",
    )
    ref.add_argument(
        "--reference", metavar="PATH", default=DEFAULT_REFERENCE,
        help=f"reference genome FASTA (default: {DEFAULT_REFERENCE})",
    )
    ref.add_argument(
        "--gtf", metavar="PATH", default=DEFAULT_GTF,
        help=f"annotation GTF (default: {DEFAULT_GTF})",
    )
    ref.add_argument(
        "--star-index", metavar="PATH", default=DEFAULT_STAR_INDEX,
        help=f"STAR index directory (default: {DEFAULT_STAR_INDEX})",
    )
    ref.add_argument(
        "--sjdb-overhang", type=int, default=74, metavar="N",
        help="STAR splice junction overhang = read_length - 1 (default: 74 for 75bp reads)",
    )
    ref.add_argument(
        "--genome-sa-index-nbases", type=int, default=14, metavar="N",
        help="STAR genome index SA size; 14 for full genome, 11 for small references (default: 14)",
    )

    args = parser.parse_args()

    # ── Validate ───────────────────────────────────────────────────────────
    if args.mode == "bulk":
        if not args.sample and not args.samples:
            parser.error("--mode bulk requires --sample or --samples")
        if args.sample and not args.run_id:
            parser.error("--run-id is required with --sample")
        if args.samples and not args.run_id_suffix:
            parser.error("--run-id-suffix is required with --samples")
        if args.samples and args.run_id:
            parser.error("--run-id cannot be used with --samples; use --run-id-suffix instead")

    elif args.mode == "scrna":
        if not args.demux_run_id:
            parser.error("--demux-run-id is required for --mode scrna")
        if not args.bam and not args.cellsnp_dir:
            parser.error("--mode scrna requires either --bam (+ --barcodes) or --cellsnp-dir")
        if args.bam and not args.barcodes:
            parser.error("--barcodes is required when --bam is provided")

    # ── Dispatch ───────────────────────────────────────────────────────────
    if args.mode == "bulk":
        if args.sample:
            run_single(args.run_id, args.sample, args, args.fastq_dir)
        else:
            run_batch(args.samples, args.run_id_suffix, args)
    else:
        run_scrna(args)


if __name__ == "__main__":
    main()
