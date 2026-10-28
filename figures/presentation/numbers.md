# Slide numbers — generated 2026-10-28

## Cohort overview
- sc/bulk runs in DB: **44** (source: runs table, dtype != WES)
- WES runs in DB: **10**
- Total runs: **54**

## QC (Fig 5)
- Mapping metric: STAR uniquely mapped reads % (from Log.final.out)
- sc/bulk range: 68.9% – 93.0%
- Runs below 85%: 2 (SNUC2B 68.9%, VACO4S 85.0%)
- WES BWA mapping rates: 99.42% – 99.57%
- 85% Snakemake gate: **never fires for sc/bulk** (flagstat on STAR BAM always 100%)

## Same-donor vs different-donor (Fig 3)
- Same-donor pairs: **8** (8 biological pairs across 4 donor groups)
- Jaccard  same-donor median: **0.5289**
- Jaccard  diff-donor median: **0.3256**
- **Clean separation: min same-donor 0.4192 > max diff-donor 0.4007**
- Containment same-donor median: **0.7554**
- Containment diff-donor median: **0.5029**

## Pool_ctr cell assignment (Fig 6)
- Guided run: pool_ctr_v5   total cells: **10,308**
- Assigned singlets: **80.9%** (8,338 cells)
- Doublets: **6.8%**
- Unassigned: **12.3%**
- Per-line counts (guided, as % of all cells):
    - RKO: 4,502 (43.7%)
    - HCT116: 2,968 (28.8%)
    - HT29: 722 (7.0%)
    - CACO2: 80 (0.8%)
    - DLD1: 33 (0.3%)
    - LIM1215: 30 (0.3%)
    - HCT8: 3 (0.0%)
- **Pool composition is highly unequal: RKO ~44%, HCT8 <1%**
- Unguided run: pool_ctr_v4   total cells: **10,308**
- Unguided assigned singlets: **63.1%** (3/7 lines)

## Donor matching (Fig 8)
- High-confidence matches: 2/7  (RKO, HCT116)
- Low-confidence matches: 1/7  (HT29)
- Best concordance: 0.8961 (RKO)
- **min_gap used: ≈0.10** (inferred from data; no run log; config default 0.15 contradicts data)
- min_concordance used: ≈0.80 (post-backfill default)

## Sample identification (Fig 9)
- Query: SRR5071686_hg38 (RKO sc/bulk run)
- Self-match overlap: **100.0%**
- Next best sc/bulk (HCT116): **44.3%**
- Gap: **56 pp**

## WES similarity (Fig 10)
- All-pairs Jaccard median: **0.2182**
- All-pairs range: 0.1881 – 0.3212
- IKPT vs IKPM Jaccard: **0.3116** (bg median 0.2182)
- NCO vs NCO_A Jaccard: **0.3212** (bg median 0.2182)
