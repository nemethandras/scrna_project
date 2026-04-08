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

From the project root:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

Add `-n` for a dry run (shows what would execute without running anything):

```bash
snakemake --snakefile workflow/Snakefile --cores 4 -n
```

## Setting up a new run

### 1. Add your FASTQ files

Place paired-end reads in a dedicated folder. The filenames must follow the pattern `<sample_id>_1.fastq` / `<sample_id>_2.fastq`:

```
data/
└── MY_SAMPLE/
    ├── MY_SAMPLE_1.fastq
    └── MY_SAMPLE_2.fastq
```

### 2. Prepare reference files (skip if reusing an existing one)

Put the genome FASTA and GTF in `data/reference/`. Name the STAR index directory to reflect the genome and read length so different indices don't get confused:

```
data/reference/
├── hg38.fa
├── gencode.v47.gtf
└── star_index_hg38_oh99/    # will be built automatically on first run
```

The naming convention for index directories is `star_index_<genome>_oh<overhang>`, where overhang = read length − 1.

### 3. Edit `config/config.yaml`

**Always change these:**

```yaml
run_id:    MY_SAMPLE_hg38     # unique label — all outputs go under results/<run_id>/
samples:
  - MY_SAMPLE                 # must match the FASTQ filename prefix
fastq_dir: data/MY_SAMPLE     # folder containing the FASTQ files
sequencing: paired            # paired (needs _1 + _2) or single (needs _1 only)
```

**Change these only when switching to a different reference or read length:**

```yaml
reference_genome: data/reference/hg38.fa
annotation_gtf:   data/reference/gencode.v47.gtf
star_index_dir:   data/reference/star_index_hg38_oh99/

star:
  sjdb_overhang: 99           # read_length - 1
  genome_sa_index_nbases: 14  # min(14, floor(log2(genome_size) / 2 - 1))
                              # STAR will print the correct value if wrong
```

If `star_index_dir` points to a path that does not exist yet, Snakemake will build the index automatically before alignment.

### 4. Run

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

Results land in `results/<run_id>/` and logs in `logs/<run_id>/`. Previous runs are never touched.

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
