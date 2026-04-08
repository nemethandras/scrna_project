# scRNA-seq Variant Calling Pipeline

A Snakemake pipeline for quality control, alignment, and variant calling from RNA-seq paired-end reads.

## Pipeline overview

```
FASTQ (R1 + R2)
  └─ FastQC          → quality report per read
  └─ STAR align      → BAM (unsorted)
       └─ samtools sort + index  → sorted BAM
            └─ samtools flagstat → mapping rate QC (fails if below threshold)
            └─ bcftools mpileup
                 └─ bcftools call   → raw VCF
                      └─ bcftools filter → filtered VCF  ✓ final output
```

## Dependencies

- [Snakemake](https://snakemake.readthedocs.io) >= 7
- [STAR](https://github.com/alexdobin/STAR)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [samtools](http://www.htslib.org/)
- [bcftools](http://www.htslib.org/)

## Project structure

```
scrna_project/
├── config/
│   └── config.yaml          # all run parameters — edit this before each run
├── workflow/
│   └── Snakefile
├── data/
│   ├── reference/           # genome FASTA, GTF, and STAR indices
│   │   ├── chr22.fa
│   │   ├── gencode.v47.chr22.gtf
│   │   └── star_index_chr22_oh74/   # index name encodes genome + overhang
│   └── SRR1258218/          # one folder per sample dataset
│       ├── SRR1258218_1.fastq
│       └── SRR1258218_2.fastq
├── results/
│   └── <run_id>/            # all outputs isolated per run
│       ├── fastqc/
│       ├── star/<sample>/
│       ├── bam/
│       └── vcf/
└── logs/
    └── <run_id>/            # logs mirror the results structure
```

## Setup

All dependencies are managed via conda. With [conda](https://docs.conda.io/en/latest/miniconda.html) installed, create and activate the environment once:

```bash
conda env create -f environment.yml
conda activate scrna
```

Activate the environment at the start of every session before running the pipeline.

## Running the pipeline

Place your FASTQ files in a folder under `data/` and run:

```bash
python run_pipeline.py --run-id MY_SAMPLE_chr22 --samples MY_SAMPLE --fastq-dir data/MY_SAMPLE
```

Results land in `results/<run-id>/` and logs in `logs/<run-id>/`. Previous runs are never touched.

**All options:**

```
--run-id       unique label for this run (required)
--samples      one or more sample IDs matching FASTQ prefixes (required)
--fastq-dir    directory containing the FASTQ files (required)
--sequencing   paired (default) or single
--cores        number of CPU cores (default: 16)
-n/--dry-run   show what would run without executing anything

reference overrides (optional — fall back to config/config.yaml):
--reference    path to reference genome FASTA
--gtf          path to annotation GTF
--star-index   path to STAR index directory
```

**Examples:**

```bash
# single-end
python run_pipeline.py --run-id my_run --samples SAMPLE1 --fastq-dir data/SAMPLE1 --sequencing single

# multiple samples
python run_pipeline.py --run-id batch1 --samples S1 S2 S3 --fastq-dir data/batch1

# custom reference
python run_pipeline.py --run-id hg38_run --samples SAMPLE1 --fastq-dir data/SAMPLE1 \
    --reference data/reference/hg38.fa \
    --gtf data/reference/gencode.v47.gtf \
    --star-index data/reference/star_index_hg38_oh99

# dry run
python run_pipeline.py --run-id test --samples SRR1258218 --fastq-dir data/SRR1258218 --dry-run
```

If `--star-index` points to a directory that does not exist yet, Snakemake will build the index automatically before alignment.

## Preparing input data

### FASTQ files

Place reads in a dedicated folder. Filenames must follow the pattern `<sample_id>_1.fastq` (and `<sample_id>_2.fastq` for paired-end):

```
data/
└── MY_SAMPLE/
    ├── MY_SAMPLE_1.fastq
    └── MY_SAMPLE_2.fastq
```

### Reference files

Put the genome FASTA and GTF in `data/reference/`. Name the STAR index directory to reflect the genome and read length:

```
data/reference/
├── hg38.fa
├── gencode.v47.gtf
└── star_index_hg38_oh99/    # built automatically on first run
```

Naming convention: `star_index_<genome>_oh<overhang>`, where overhang = read length − 1.

### Changing read length or genome

Edit these values in `config/config.yaml` before running:

```yaml
star:
  sjdb_overhang: 99           # read_length - 1
  genome_sa_index_nbases: 14  # min(14, floor(log2(genome_size) / 2 - 1))
                              # STAR will print the correct value if wrong
```

## Outputs

| File | Description |
|---|---|
| `results/<run_id>/fastqc/<sample>_[1,2]_fastqc.html` | Per-read quality report |
| `results/<run_id>/bam/<sample>.sorted.bam` | Sorted, indexed alignment |
| `results/<run_id>/bam/<sample>.flagstat.txt` | Mapping rate summary |
| `results/<run_id>/vcf/<sample>.raw.vcf` | Unfiltered variants |
| `results/<run_id>/vcf/<sample>.filtered.vcf` | Final filtered variants |

## Tuning parameters

The parameters below in `config.yaml` have sensible defaults and only need changing if you have a specific reason:

| Parameter | Default | Effect |
|---|---|---|
| `samtools.min_mapping_rate` | 85 | Pipeline aborts if mapping rate falls below this (%) |
| `bcftools.min_base_quality` | 20 | Minimum base quality to count a read at a position |
| `bcftools.min_mapping_quality` | 20 | Minimum mapping quality to include a read |
| `filters.min_qual` | 10 | Minimum variant quality score (QUAL in VCF) |
| `filters.min_depth` | 3 | Minimum read depth at a variant site |
