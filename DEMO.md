# Pipeline demo

All commands run from the project root. Activate the environment once at the start:

```bash
cd /groups/nils/members/andras/scrna_project
conda activate scrna
source .env
```

---

## 1. Database — what's already loaded

Quick summary of everything in the DB:

```bash
python - <<'EOF'
import sqlite3, pandas as pd
conn = sqlite3.connect("results/variants.db")

print("=== Database summary ===")
for tbl, n in conn.execute("""
    SELECT 'runs',           COUNT(*) FROM runs           UNION ALL
    SELECT 'samples',        COUNT(*) FROM samples        UNION ALL
    SELECT 'variants',       COUNT(*) FROM variants       UNION ALL
    SELECT 'genotype_calls', COUNT(*) FROM genotype_calls
"""):
    print(f"  {tbl:<20} {n:>8,} rows")

print("\n=== Cell lines ===")
df = pd.read_sql("SELECT cell_line, COUNT(*) AS runs FROM samples WHERE cell_line IS NOT NULL GROUP BY cell_line ORDER BY cell_line", conn)
print(df.to_string(index=False))
EOF
```

**Expected output:** 44 runs · 44 samples · ~272k variant sites · ~1.95M genotype calls across 7 cell lines (CACO2, DLD1, HCT116, HCT8, HT29, LIM1215, RKO) plus others.

---

## 2. Example SQL query

Which cell lines have the most called ALT variants, and what is the typical depth?

```bash
python - <<'EOF'
import sqlite3, pandas as pd
conn = sqlite3.connect("results/variants.db")
df = pd.read_sql("""
    SELECT s.cell_line,
           COUNT(*)                    AS n_positions,
           ROUND(AVG(gc.depth), 1)    AS avg_depth,
           ROUND(AVG(gc.allele_freq), 3) AS avg_af
    FROM genotype_calls gc
    JOIN runs    r ON gc.run_id      = r.run_id
    JOIN samples s ON r.sample_id   = s.sample_id
    WHERE gc.genotype != '0/0'
      AND s.cell_line IS NOT NULL
    GROUP BY s.cell_line
    ORDER BY n_positions DESC
""", conn)
print(df.to_string(index=False))
EOF
```

---

## 3. Bulk pipeline — dry run

Shows exactly what Snakemake would do for a new RNA-seq run without executing anything.
The sample FASTQ is already on disk; this is what a real first-time run looks like:

```bash
python run_pipeline.py --mode bulk \
    --run-id SRR5071686_hg38 \
    --sample SRR5071686 \
    --dry-run
```

For a **new** sample you haven't run before, drop `--dry-run` and it executes:

```bash
python run_pipeline.py --mode bulk \
    --run-id SRR5071686_hg38 \
    --sample SRR5071686
```

Check progress any time:
```bash
tail -f logs/batch.log
```

---

## 4. scRNA demultiplexing — dry run

The Pool_ctr CellRanger BAM is already on disk. This re-runs the demux from
scratch (cellsnp → Vireo → match → merge) without writing any files:

```bash
python run_pipeline.py --mode scrna \
    --bam data/Pool_ctr/possorted_genome_bam.bam \
    --barcodes data/Pool_ctr/barcodes.tsv.gz \
    --demux-run-id pool_ctr_demo \
    --n-donors 7 \
    --dry-run
```

To skip the cellsnp step and reuse the existing pileup (much faster):

```bash
python run_pipeline.py --mode scrna \
    --cellsnp-dir data/Pool_ctr/cellsnp \
    --demux-run-id pool_ctr_demo \
    --n-donors 7 \
    --dry-run
```

The most recent real result is `pool_ctr_v4`. Inspect the final assignments:

```bash
head -5 results/demux/pool_ctr_v4/final_assignments.tsv | column -t
wc -l results/demux/pool_ctr_v4/final_assignments.tsv
```

---

## 5. Genotype matching — two scenarios

### Scenario A: sample IS in the database → exact match

Feed a bulk RNA-seq VCF back into `match_vcf.py`. The tool computes how many of
the query's ALT variants are found in each DB sample.

```bash
conda run -n scrna python scripts/match_vcf.py \
    --vcf results/SRR5071686_hg38/vcf/SRR5071686.panel_genotyped.vcf \
    --metric overlap
```

Expected output — RKO is found immediately at 100%, all other cell lines cluster
well below it (≤ 44%):

```
  Run                        Sample             Match    Shared     In DB
  -----------------------------------------------------------------------
  SRR5071686_hg38            SRR5071686        100.0%    11,639    39,186
  SRR5071667_hg38            SRR5071667         44.3%     5,161    46,862
  SRR5071662_hg38            SRR5071662         43.9%     5,105    51,257
  ...
```

The 56-point gap between rank 1 and rank 2 is unambiguous — the database
identifies the sample correctly.

---

### Scenario B: sample NOT in the database → nearest match, no exact match

Pool_ctr is a mixture of 7 cell lines. Vireo clusters cells into anonymous donors
without knowing the genotypes — its inferred donor genotypes are then compared
against every reference in the DB. Because the pool contains cell lines not well
represented in the DB (under-sequenced / no reference run for some donors),
several donors cannot be confidently assigned.

```bash
conda run -n scrna python scripts/match_vireo.py \
    --vireo-dir results/demux/pool_ctr_v4/vireo \
    --db results/variants.db \
    --output /tmp/donor_matches_demo.tsv \
    --min-concordance 0.80 \
    --min-gap 0.10 \
    --unique

cat /tmp/donor_matches_demo.tsv | column -t
```

Expected output:

```
vireo_donor  assigned_line    cell_line  concordance  n_positions  second_line      second_concordance  gap      confidence
donor0       no_match                    0.709        45930        SRR5071686_hg38  0.834               -0.125   no_match
donor1       no_match                    0.746        31861        SRR5071667_hg38  0.789               -0.044   no_match
donor2       SRR5071686_hg38  RKO        0.896        33729        SRR5071663_hg38  0.700                0.196   high
donor3       SRR5071672_hg38  HT29       0.825        30373        SRR5071680_hg38  0.745                0.081   low
donor4       no_match                    0.741        33752        SRR5071693_hg38  0.743               -0.002   no_match
donor5       no_match                    0.719        39585        SRR5071686_hg38  0.827               -0.108   no_match
donor6       SRR5071667_hg38  HCT116     0.869        40079        SRR5071691_hg38  0.750                0.119   high
```

- **donor2 → RKO** and **donor6 → HCT116**: high-confidence matches (concordance >> threshold, large gap)
- **donor3 → HT29**: assigned but low-confidence — gap too small to be certain
- **donor0, 1, 4, 5**: `no_match` — best available reference shown in `second_line`, but the
  concordance gap is negative (another reference scores equally or better after 1:1 assignment).
  These donors likely represent cell lines present in the pool but with too few reference variants
  in the DB to exceed the 0.80 threshold.

---

## 7. Explore results in the notebook

```bash
jupyter lab notebooks/01_explore_fastq.ipynb
```

The notebook connects to `results/variants.db` automatically. Sections to show:

- **Cross-sample QC overview** — mapping rates and variant counts, all 44 runs
- **Variant overlap heatmap** — Jaccard / containment between cell lines + distribution curves
- **scRNA demultiplexing results** — select `pool_ctr_v4` from the dropdown
- **Validation — reference-guided Vireo** — ground-truth composition of Pool_ctr
- **Concordance profiles** — full genotype concordance matrix for `pool_ctr_v4`
- **Raw SQL explorer** — paste any query and run it live
