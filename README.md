# scRNA-seq Variant Calling Pipeline

A Snakemake pipeline for quality control, alignment, variant calling, and database loading from RNA-seq reads (single-end or paired-end).

## Pipeline overview

```mermaid
flowchart TD
    A([FASTQ]) --> B[FastQC\nquality report]
    A          --> C[STAR align\nBAM unsorted]

    C --> D[samtools sort\nsorted BAM]
    D --> E[samtools index\nBAM index]
    D --> F[samtools flagstat\nmapping rate QC\n⚠ aborts if below threshold]

    E --> G[bcftools mpileup\npileup BCF]
    G --> H[bcftools call\nraw VCF]
    H --> I[bcftools filter\nfiltered VCF]

    C -->|Log.final.out| J
    F --> J
    I --> J[load_to_database]

    J --> K[(SQLite\nvariants.db)]
    J --> L[(Grist\nRuns · QC_Summary · Samples)]

    style K fill:#2d6a4f,color:#fff
    style L fill:#1d3557,color:#fff
    style F fill:#e63946,color:#fff
```

## Database structure

```mermaid
erDiagram
    samples {
        INTEGER sample_id PK
        TEXT    name
        TEXT    date_added
    }
    runs {
        TEXT    run_id PK
        INTEGER sample_id FK
        TEXT    reference
        TEXT    sequencing
        TEXT    run_date
        REAL    mapping_rate
        INTEGER total_raw
        INTEGER total_filt
    }
    variants {
        INTEGER variant_id PK
        TEXT    chromosome
        INTEGER position
        TEXT    ref_allele
        TEXT    alt_allele
        TEXT    variant_type
    }
    genotype_calls {
        INTEGER call_id PK
        TEXT    run_id FK
        INTEGER variant_id FK
        TEXT    genotype
        REAL    quality
        INTEGER depth
        INTEGER ref_depth
        INTEGER alt_depth
        REAL    allele_freq
        TEXT    filter_status
    }

    samples     ||--o{ runs           : "one sample, many runs"
    runs        ||--o{ genotype_calls : "one run, many calls"
    variants    ||--o{ genotype_calls : "one position, called across samples"
```

## Dependencies

- [Snakemake](https://snakemake.readthedocs.io) >= 7
- [STAR](https://github.com/alexdobin/STAR)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [samtools](http://www.htslib.org/)
- [bcftools](http://www.htslib.org/)
- [cyvcf2](https://github.com/brentp/cyvcf2)
- [grist_api](https://github.com/gristlabs/grist-api)
- [python-dotenv](https://github.com/theskumar/python-dotenv)

All managed via conda — see Setup below.

## Project structure

```
scrna_project/
├── config/
│   └── config.yaml          # all run parameters — edit this before each run
├── workflow/
│   └── Snakefile
├── scripts/
│   ├── load_to_db.py        # loads results into SQLite and Grist
│   └── test_grist.py        # inspect Grist table schema
├── data/
│   ├── reference/           # genome FASTA, GTF, and STAR indices
│   └── <sample>/            # one folder per sample: <sample>_1.fastq [_2.fastq]
├── results/
│   ├── <run_id>/            # all outputs isolated per run
│   │   ├── fastqc/
│   │   ├── star/<sample>/
│   │   ├── bam/
│   │   ├── vcf/
│   │   └── db/              # touch files confirming db load completed
│   └── variants.db          # SQLite database accumulating all runs
├── logs/
│   └── <run_id>/            # logs mirror the results structure
├── .env                     # API credentials — never committed
└── .env.example             # template showing which variables are needed
```

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yml
conda activate scrna
```

Activate at the start of every session before running the pipeline.

### 2. API credentials

Copy `.env.example` to `.env` and fill in your values:

```bash
cp .env.example .env
```

```bash
# .env
export GRIST_SERVER=https://your-grist-instance.example.com
export GRIST_API_KEY=your_api_key_here
export GRIST_DOC_ID=your_doc_id_here
```

Source the file before running the pipeline:

```bash
source .env
```

The `.env` file is in `.gitignore` and will never be committed. The SQLite database (`results/variants.db`) is also ignored.

## Running the pipeline

Place your FASTQ files in a folder under `data/` and run:

```bash
source .env && python run_pipeline.py --run-id MY_RUN --sample MY_SAMPLE
```

The pipeline runs **in the background by default** — it detaches from your terminal so it keeps running even if you close your laptop or the SSH connection drops. Output is written to `logs/nohup_<sample>.log`. Results land in `results/<run-id>/`. Previous runs are never touched.

Check progress at any time:

```bash
tail -f logs/nohup_MY_SAMPLE.log
```

### Options

```
--run-id            unique label for this run (required)
--sample            sample ID matching the FASTQ filename prefix (required)
--fastq-dir         directory containing the FASTQ files (default: data/<sample>)
--sequencing        paired or single (default: single)
--cores             number of CPU cores (default: 16)
--db                SQLite database path (default: results/variants.db)
--force             re-run all steps even if outputs already exist
--rerun-incomplete  re-run jobs with incomplete outputs from a previous failed run
--foreground        run in the terminal instead of detaching (useful for debugging)
-n/--dry-run        show what would run without executing anything (always foreground)

reference overrides (optional — defaults to config/config.yaml):
--reference         path to reference genome FASTA
--gtf               path to annotation GTF
--star-index        path to STAR index directory
```

### Examples

```bash
# standard run — detaches immediately, safe to close laptop
source .env && python run_pipeline.py --run-id SRR5071697_hg38 --sample SRR5071697

# paired-end
source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --sequencing paired

# dry run — shows what would execute without running anything
python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --dry-run

# force re-run of all steps (e.g. after pipeline changes)
source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --force

# run in foreground to watch output live (e.g. for debugging)
source .env && python run_pipeline.py --run-id SRR5071686_hg38 --sample SRR5071686 --foreground
```

### Running only the database step

If upstream results already exist (e.g. adding a sample that was processed before the db step was introduced), call the loading script directly — no need to re-run the full pipeline:

```bash
source .env && conda run -n scrna python scripts/load_to_db.py \
    --run-id MY_RUN \
    --sample MY_SAMPLE \
    --vcf results/MY_RUN/vcf/MY_SAMPLE.filtered.vcf \
    --flagstat results/MY_RUN/bam/MY_SAMPLE.flagstat.txt
```

Re-running this script for an existing run is safe — it removes any previous Grist entries for that run before inserting, so no duplicates are created.

## Preparing input data

### FASTQ files

Place reads in a dedicated folder. Filenames must follow the pattern `<sample>_1.fastq` (and `<sample>_2.fastq` for paired-end):

```
data/
└── MY_SAMPLE/
    ├── MY_SAMPLE_1.fastq
    └── MY_SAMPLE_2.fastq   # paired-end only
```

### Reference files

Put the genome FASTA and GTF in `data/reference/`. The STAR index is built automatically on first run if the directory does not exist:

```
data/reference/
├── genome.fa
├── genes.gtf
└── star_index_hg38_oh99/   # built automatically; name encodes genome + overhang
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
| `results/<run_id>/fastqc/<sample>_1_fastqc.html` | Per-read quality report |
| `results/<run_id>/bam/<sample>.sorted.bam` | Sorted, indexed alignment |
| `results/<run_id>/bam/<sample>.flagstat.txt` | Mapping rate summary |
| `results/<run_id>/vcf/<sample>.raw.vcf` | Unfiltered variants |
| `results/<run_id>/vcf/<sample>.filtered.vcf` | Final filtered variants |
| `results/<run_id>/db/<sample>.loaded` | Touch file confirming db load completed |
| `results/variants.db` | SQLite database with all runs, variants, and genotype calls |

## Tuning parameters

| Parameter | Default | Effect |
|---|---|---|
| `samtools.min_mapping_rate` | 85 | Pipeline aborts if mapping rate falls below this (%) |
| `bcftools.min_base_quality` | 20 | Minimum base quality to count a read at a position |
| `bcftools.min_mapping_quality` | 20 | Minimum mapping quality to include a read |
| `filters.min_qual` | 30 | Minimum variant quality score (QUAL in VCF) |
| `filters.min_depth` | 10 | Minimum read depth at a variant site |
| `db_path` | `results/variants.db` | SQLite database file location |
