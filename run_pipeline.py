#!/usr/bin/env python3
"""Command-line wrapper for the scRNA-seq variant calling pipeline."""

import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="scRNA-seq variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # paired-end run
  python run_pipeline.py --run-id SRR1258218_chr22 --samples SRR1258218 --fastq-dir data/SRR1258218

  # single-end run
  python run_pipeline.py --run-id my_run --samples SAMPLE1 --fastq-dir data/SAMPLE1 --sequencing single

  # multiple samples
  python run_pipeline.py --run-id batch1 --samples S1 S2 S3 --fastq-dir data/batch1

  # dry run (shows what would execute without running anything)
  python run_pipeline.py --run-id test --samples SRR1258218 --fastq-dir data/SRR1258218 --dry-run

  # custom reference
  python run_pipeline.py --run-id hg38_run --samples SAMPLE1 --fastq-dir data/SAMPLE1 \\
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
        "--samples", required=True, nargs="+",
        help="one or more sample IDs matching FASTQ filename prefixes",
    )
    parser.add_argument(
        "--fastq-dir", required=True,
        help="directory containing the FASTQ files",
    )
    parser.add_argument(
        "--sequencing", choices=["paired", "single"], default="paired",
        help="sequencing layout (default: paired)",
    )
    parser.add_argument(
        "--cores", type=int, default=16,
        help="number of CPU cores to use (default: 16)",
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true",
        help="show what would run without executing anything",
    )

    ref = parser.add_argument_group(
        "reference overrides",
        "optional — falls back to values in config/config.yaml",
    )
    ref.add_argument("--reference", metavar="PATH", help="reference genome FASTA")
    ref.add_argument("--gtf",       metavar="PATH", help="annotation GTF")
    ref.add_argument("--star-index", metavar="PATH", help="STAR index directory")

    args = parser.parse_args()

    # Snakemake parses --config values as YAML, so a bracketed list works for samples
    samples_yaml = "[" + ",".join(args.samples) + "]"
    config_overrides = [
        f"run_id={args.run_id}",
        f"fastq_dir={args.fastq_dir}",
        f"sequencing={args.sequencing}",
        f"samples={samples_yaml}",
    ]
    if args.reference:
        config_overrides.append(f"reference_genome={args.reference}")
    if args.gtf:
        config_overrides.append(f"annotation_gtf={args.gtf}")
    if args.star_index:
        config_overrides.append(f"star_index_dir={args.star_index}")

    cmd = [
        "snakemake",
        "--snakefile", "workflow/Snakefile",
        "--cores", str(args.cores),
        "--config", *config_overrides,
    ]
    if args.dry_run:
        cmd.append("-n")

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
