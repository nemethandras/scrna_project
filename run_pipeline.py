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


def main():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="scRNA-seq variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # minimal — single-end, fastq-dir derived automatically
  python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686

  # paired-end
  python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --sequencing paired

  # override fastq location
  python run_pipeline.py --run-id my_run --sample SAMPLE1 --fastq-dir /data/elsewhere

  # dry run (shows what would execute without running anything)
  python run_pipeline.py --run-id test --sample SRR5071686 --dry-run

  # custom db path
  python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --db /path/to/variants.db

  # custom reference
  python run_pipeline.py --run-id hg38_run --sample SAMPLE1 \\
      --reference data/reference/hg38.fa \\
      --gtf data/reference/gencode.v47.gtf \\
      --star-index data/reference/star_index_hg38_oh99
        """,
    )

    parser.add_argument(
        "--run-id", required=True,
        help="unique label for this run — outputs go to results/<run-id>/",
    )
    parser.add_argument(
        "--sample", required=True,
        help="sample ID matching the FASTQ filename prefix",
    )
    parser.add_argument(
        "--fastq-dir", default=None, metavar="PATH",
        help="directory containing the FASTQ files (default: data/<sample>)",
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
        help="run in the foreground instead of detaching from the terminal",
    )
    parser.add_argument(
        "--no-db", action="store_true",
        help="skip the load_to_database step (results are not added to SQLite or Grist)",
    )
    parser.add_argument(
        "--db", metavar="PATH", default=None,
        help="SQLite database file for variant storage (default: results/variants.db)",
    )

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

    fastq_dir = args.fastq_dir if args.fastq_dir else f"data/{args.sample}"

    samples_yaml = f"[{args.sample}]"
    config_overrides = [
        f"run_id={args.run_id}",
        f"samples={samples_yaml}",
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

    if args.foreground or args.dry_run:
        result = subprocess.run(cmd)
        sys.exit(result.returncode)

    log_path = f"logs/nohup_{args.sample}.log"
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


if __name__ == "__main__":
    main()
