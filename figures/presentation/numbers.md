# Slide numbers — generated 2026-10-28

## Cohort overview
- sc/bulk runs in DB: **44** (source: runs table, dtype != WES)
- WES runs in DB: **10**
- Total runs: **54**

## QC (Fig 5)
- sc/bulk mapping rates: 100.0% – 100.0%  (source: runs.mapping_rate)
- WES mapping rates:     99.42% – 99.57%
- Pipeline abort threshold: 85% (source: config.yaml samtools.min_mapping_rate)

## Same-donor vs different-donor (Fig 3)
- Same-donor Jaccard pairs: **8** (8 biological pairs)
- Same-donor Jaccard median: **0.5289**  (source: pairwise sc_jac matrix)
- Diff-donor Jaccard median: **0.3256**
- Min same-donor Jaccard: 0.4192
- Max diff-donor Jaccard: 0.4007
- Diff-donor pairs above same-donor min: 0
- Same-donor containment median: **0.7554**
- Diff-donor containment median: **0.5029**

## Pool_ctr cell assignment (Fig 6)
- Guided run: pool_ctr_v5   total cells: **10,308** (source: final_assignments.tsv)
- Guided singlets: **80.9%**
- Guided doublets: **6.8%**
- Guided unassigned: **12.3%**
- High-confidence cells: **80.9%**
- Unguided run: pool_ctr_v4   total cells: **10,308**
- Unguided singlets: **63.1%** (only 3/7 lines identified)

## Donor matching (Fig 8)
- Donors matched with high confidence: 2/7
- Max concordance (best match): 0.8961  (RKO)
- Threshold used (inferred): min_concordance=0.80, min_gap=0.15
  source: pool_ctr_v4/donor_matches.tsv — post-backfill thresholds

## WES similarity (Fig 10)
- WES all-pairs Jaccard median: **0.2182**  (source: pairwise wes_jac matrix)
- WES all-pairs Jaccard range: 0.1881 – 0.3212
- IKPT vs IKPM Jaccard: **0.3116** (vs background median 0.2182)
- NCO vs NCO_A Jaccard: **0.3212** (vs background median 0.2182)
