# Thesis Presentation Numbers — v6

Generated: 2026-11-06
Figures: `figures/presentation_v6/`  (PNG + SVG, 200 dpi, white background)

---

## Task 1 — Mapping gate audit

**Short answer: the 85 % gate cannot fail for STAR runs.**  
STAR is configured without `--outSAMunmapped Within`, so unmapped reads are
never written to the BAM.  `samtools flagstat` on that BAM always reports
100 % mapped, regardless of how many reads STAR actually aligned.
The gate (`config["samtools"]["min_mapping_rate"]` = 85 %) checks the flagstat
rate, not STAR's internal unique-mapping rate.

### STAR unique-mapping rates (44 sc/bulk runs)

| Cell line     | Run ID            | STAR unique % | flagstat % |
|---------------|-------------------|---------------|------------|
| SNUC2B        | SRR5071688_hg38   | **68.86**     | 100 %      |
| VACO4S        | SRR5071696_hg38   | **84.96**     | 100 %      |
| VACO4A        | SRR5071689_hg38   | 88.57         | 100 %      |
| DLD1          | SRR5071662_hg38   | 88.61         | 100 %      |
| *(all others)*| …                 | 89.45–92.98   | 100 %      |

All 44 runs pass the 85 % flagstat gate.

**Hypothetical failure counts (unique-mapping threshold, not the actual gate):**
- At 85 %: **2 runs** would fail (SNUC2B 68.86 %, VACO4S 84.96 %)
- At 75 %: **1 run** would fail (SNUC2B 68.86 %)
- At 70 %: **1 run** would fail (SNUC2B 68.86 %)

---

## Task 2 — Concordance threshold audit

### Original pool_ctr_v4 run (--min-concordance 0.65 --min-gap 0.10 --unique)

| Vireo donor | Assigned line | Concordance | 2nd best | Gap    | Confidence |
|-------------|---------------|-------------|----------|--------|------------|
| donor0      | no_match      | 0.7091      | 0.8341   | −0.1251 | no_match  |
| donor1      | no_match      | 0.7456      | 0.7894   | −0.0438 | no_match  |
| donor2      | RKO           | 0.8961      | 0.7001   | +0.1960 | **high**  |
| donor3      | HT29          | 0.8252      | 0.7445   | +0.0807 | low       |
| donor4      | no_match      | 0.7406      | 0.7427   | −0.0020 | no_match  |
| donor5      | no_match      | 0.7191      | 0.8271   | −0.1080 | no_match  |
| donor6      | HCT116        | 0.8692      | 0.7503   | +0.1189 | **high**  |

### Re-run with stricter thresholds (--min-concordance 0.80 --min-gap 0.15 --min-positions 200 --unique)

| Vireo donor | Assigned line | Change?                         |
|-------------|---------------|---------------------------------|
| donor0–1    | no_match      | No change                       |
| donor2      | RKO           | No change (high conf.)          |
| donor3      | HT29          | No change (stays low conf.)     |
| donor4–5    | no_match      | No change                       |
| donor6      | HCT116        | **high → low** (gap 0.1189 < 0.15) |

**Summary:** No assignment changes — both RKO and HCT116 are still identified.
Only donor6 (HCT116) drops from high to low confidence due to gap < 0.15.

### match_vcf threshold and metric
`scripts/match_vcf.py` uses **overlap** (= containment = shared / query_total)
as its metric by default; `--metric jaccard` is the only alternative.
No concordance metric is available.  Default threshold: `--min-match 40.0`
(percent containment).

### LOO analysis metric
The LOO identification analysis (notebook cell34) uses **concordance**, not
match_vcf containment.  The separation cutoff (midpoint between max diff-donor
and min same-donor ref-vs-ref concordance) is **0.886**.  This is an analysis
threshold used to characterise separability, not a parameter of match_vcf.

### PASS-only concordance check
Pairwise concordance computed with and without `filter_status = 'PASS'`
gives **identical results** because backfilled 0/0 rows always carry
`filter_status = 'PASS'`.

---

## Task 3 — Regenerated figures

All three figures in `figures/presentation_v6/`.

### fig_why_concordance.png  (horizontal bar chart)
Separation metric = margin ÷ IQR of diff-donor pairs.

| Metric       | Margin   | Norm. margin | Overlap |
|--------------|----------|--------------|---------|
| Jaccard      | +0.019   | +0.42        | 0       |
| Containment  | −0.204   | −1.74        | 1,274   |
| Concordance  | +0.114   | +4.36        | 0       |

Labels inside the bars show margin and overlap only (AUC removed per request).

### fig_thresholds.png  (4 panels)
- **(a) Containment:** diff-donor and same-donor distributions; dashed line at
  40 % (`match_vcf --min-match 40`, annotated).
- **(b) Concordance ref-vs-ref:** separating cutoff at 0.886 (midpoint between
  max_diff = 0.829 and min_same = 0.943).
- **(c) Vireo concordance:** strip plot of all 7 Vireo donors × 54 references;
  matched donors highlighted; dashed line at 0.80
  (`match_vireo --min-concordance 0.80`).
- **(d) Gap:** bar chart for pool_ctr_v4 donors; dashed line at 0.15
  (`match_vireo --min-gap 0.15`).

### fig_loo_identify.png  (SW620 LOO bar chart)
SW480 (same donor) scores 0.963; highest non-partner (LIM2099) scores 0.802.
Separating cutoff line at 0.886, labelled: *"separating cutoff (0.886) —
analysis only, not match_vcf's rule."*

---

## Task 4 — Runtime estimates

All times are wall-clock, measured from snakemake logs on the same 16-core host
(`menten.itb.biologie.hu-berlin.de`).

### Bulk RNA / scRNA reference run (STAR pipeline)
Cores: 16 (STAR), 8 (samtools/bcftools)

Measured from the 40 runs with a complete, error-free snakemake log.

| Step               | Min  | Median | Max  |
|--------------------|------|--------|------|
| STAR align         | 2.1 min | 3.4 min | 21.5 min |
| bcftools mpileup   | 2.1 min | 5.5 min | 63.4 min |
| **Total (per run)** | **6.4 min** | **13.8 min** | **92 min** |

The 92-min maximum is the one run that also built the STAR genome index
(~41 min, done once; the index is now shared). Excluding it, max = 86 min.

**Exact command (single sample):**
```bash
python run_pipeline.py --mode bulk --sample {sample} --run-id {sample}_hg38 \
    --fastq-dir data/{sample} --dry-run      # drop --dry-run to execute
```

### WES (BWA + bcftools)
Cores: 16 (BWA), 8 (bcftools/samtools)

| Run    | Wall-clock |
|--------|------------|
| B2006  | 5.4 h      |
| B2040  | 3.8 h      |
| B2057  | 3.4 h      |
| C2019  | 3.6 h      |
| C2153  | 3.1 h      |
| C2159  | 3.4 h      |
| IKPM   | 3.6 h      |
| IKPT   | 3.1 h      |

**Min 3.1 h, median 3.6 h, max 5.4 h** (8 runs).

**Exact command:**
```bash
python run_pipeline.py --mode wes --sequencing paired \
    --sample {sample} --run-id {sample}_wes_v1 --dry-run
```

### scRNA pool demultiplexing (pool_ctr_v4, 10,308 cells, 12 GB BAM)

**cellsnp-lite has never been run by this pipeline.** Every pool run used a
pre-computed pileup (`--cellsnp-dir data/Pool_ctr/cellsnp`, copied in on
2026-05-27 together with the BAM), so there is no cellsnp-lite timing.

| Step              | Wall-clock   | Evidence |
|-------------------|--------------|----------|
| cellsnp-lite      | not measured | pileup was pre-computed |
| Vireo (unguided)  | ~64 min      | assignments.tsv 15:32 → vireo/ 16:36 (2026-06-17) |
| match_vireo.py    | ~30 s        | timed directly |

**Exact command (unguided, from the existing pileup):**
```bash
python run_pipeline.py --mode scrna --cellsnp-dir data/Pool_ctr/cellsnp \
    --demux-run-id {demux_run_id} --n-donors 7 --unique --dry-run
```
`--dry-run` prints the full shell script (demux → Vireo → match_vireo → merge)
without running it. To run only the matching step:
```bash
python scripts/match_vireo.py --vireo-dir results/demux/pool_ctr_v4/vireo \
    --db /groups/nils/data/scrna_pipeline/db/variants.db \
    --min-concordance 0.80 --min-gap 0.15 --min-positions 200 --unique \
    --output /tmp/donor_matches.tsv
```
`match_vireo.py` and `match_vcf.py` have no `--dry-run` mode.

---

## Key numbers summary

| Metric              | Value                           |
|---------------------|---------------------------------|
| Runs in DB (sc/bulk)| 44                              |
| Runs in DB (WES)    | 10                              |
| Mapping gate (flagstat) | 85 % — cannot fail for STAR BAMs |
| STAR unique rate range | 68.86 %–92.98 %              |
| Concordance same-donor min | 0.943                  |
| Concordance diff-donor max | 0.829                  |
| Separating cutoff   | 0.886                           |
| Concordance margin  | +0.114 (norm. +4.36 IQR units)  |
| Jaccard margin      | +0.019 (norm. +0.42 IQR units)  |
| Containment margin  | −0.204 (norm. −1.74 IQR units, 1,274 overlaps) |
| pool_ctr_v4: donors identified | 2/7 (RKO high, HCT116 low) |
| LOO SW620 top score | 0.963 (SW480, same donor)        |
| LOO SW620 runner-up | 0.802 (LIM2099)                  |
