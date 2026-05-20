#!/usr/bin/env python3
"""Command-line wrapper for the scRNA-seq variant calling pipeline."""

import argparse
import os
import subprocess
import sys

# Paths to the reference data bundled with the workflow
DEFAULT_REFERENCE  = "data/reference/genome.fa"
DEFAULT_GTF        = "data/reference/genes.gtf"
DEFAULT_STAR_INDEX = "data/reference/star_index_hg38_oh99"


def build_snakemake_cmd(run_id, sample, args, fastq_dir=None):
    fastq_dir = fastq_dir or f"data/{sample}"
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


def main():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="scRNA-seq variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # single sample
  source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686

  # batch — run-id is auto-generated as {sample}{suffix}
  source .env && python run_pipeline.py --samples SRR5071662 SRR5071667 SRR5071672 --run-id-suffix _hg38

  # batch with --no-db for test samples
  source .env && python run_pipeline.py --samples SRR5071692 --run-id-suffix _hg38 --no-db

  # dry run
  python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --dry-run

  # force re-run of all steps
  source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --force

  # run in foreground to watch output live
  source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --foreground
        """,
    )

    # ── Sample / run identity ──────────────────────────────────────────────
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--sample", metavar="ID",
        help="single sample ID — use with --run-id",
    )
    mode.add_argument(
        "--samples", nargs="+", metavar="ID",
        help="one or more sample IDs for a batch run — use with --run-id-suffix",
    )

    parser.add_argument(
        "--run-id", metavar="ID",
        help="run label for a single sample (required with --sample)",
    )
    parser.add_argument(
        "--run-id-suffix", metavar="SUFFIX",
        help="suffix appended to each sample ID to form the run ID (required with --samples)",
    )

    # ── Common options ─────────────────────────────────────────────────────
    parser.add_argument(
        "--fastq-dir", default=None, metavar="PATH",
        help="FASTQ directory for single-sample runs (default: data/<sample>)",
    )
    parser.add_argument(
        "--sequencing", choices=["paired", "single"], default="single",
        help="sequencing layout (default: single)",
    )
    parser.add_argument(
        "--cores", type=int, default=16,
        help="number of CPU cores to use (default: 16)",
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
        help="run in the foreground instead of detaching (single sample only)",
    )
    parser.add_argument(
        "--notify-email", metavar="ADDR",
        help="send an email when the batch finishes or fails (requires server mail)",
    )
    parser.add_argument(
        "--no-db", action="store_true",
        help="skip the load_to_database step (results not added to SQLite or Grist)",
    )
    parser.add_argument(
        "--db", metavar="PATH", default=None,
        help="SQLite database file (default: results/variants.db)",
    )

    # ── Reference overrides ────────────────────────────────────────────────
    ref = parser.add_argument_group(
        "reference overrides",
        "optional — defaults to the hg38 reference bundled in workflow/",
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

    # ── Validate argument combinations ─────────────────────────────────────
    if args.sample and not args.run_id:
        parser.error("--run-id is required with --sample")
    if args.samples and not args.run_id_suffix:
        parser.error("--run-id-suffix is required with --samples")
    if args.samples and args.run_id:
        parser.error("--run-id cannot be used with --samples; use --run-id-suffix instead")

    # ── Dispatch ───────────────────────────────────────────────────────────
    if args.sample:
        run_single(args.run_id, args.sample, args, args.fastq_dir)
    else:
        run_batch(args.samples, args.run_id_suffix, args)


if __name__ == "__main__":
    main()
