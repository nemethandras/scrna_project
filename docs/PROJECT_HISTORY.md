# scrna_project — project history and handoff

Written 2026-11-17 at the end of a long Claude Code session, so that a new
Claude session (new account) can continue without losing context. Read this
together with `README.md` (user-facing documentation) and `CLAUDE.md`.

---

## 1. What the project is

A Snakemake pipeline (Andras Nemeth, Nils lab, HU Berlin; server `menten`) that

1. builds a **genotype fingerprint** for each cell line from bulk RNA-seq
   (STAR) or WES (BWA) at ~7M common SNPs (1000 Genomes, AF > 5 %, hg38), and
   stores genotypes (including 0/0) in a **shared SQLite database**;
2. **demultiplexes pooled scRNA-seq** (10x) by comparing per-cell allele counts
   (cellsnp-lite output) with those fingerprints: a binomial scorer
   (`scripts/demux.py`) plus Vireo, either **guided** (known lines → donor VCF
   from the DB) or **unguided** (N donors, then `match_vireo.py` matches each
   donor to the DB by concordance).

It is the basis of Andras's thesis and is now being rolled out to lab
colleagues on the same server.

Entry point: `run_pipeline.py --mode bulk|wes|scrna` (or the shared launcher
`scrna-pipeline`). Snakemake handles bulk/wes; scrna mode builds and runs a
bash script (demux.py → [make_donor_vcf] → vireo → [match_vireo] → merge_demux).

---

## 2. Where everything lives

| What | Path |
|---|---|
| Development repo (this one) | `/groups/nils/members/andras/scrna_project` |
| GitHub remote | `https://github.com/nemethandras/scrna_project` |
| Shared root | `/groups/nils/data/scrna_pipeline/` (group `nils`, setgid) |
| Shared reference (read-only, 30 GB) | `…/reference/hg38/` — genome.fa(+fai), genes.gtf, star_index_hg38_oh99, bwa_index_hg38, raw + normalised 1000G SNP panel |
| Shared database | `…/db/variants.db` (+ `variants.db.lock`, `backups/`) |
| Shared stable code copy | `…/code/` — a git clone **of the local dev repo**, tracking branch **v4** |
| Shared conda env | `…/env/` (built from `environment.yml`, plus `sqlite`) |
| Launcher | `…/bin/scrna-pipeline` → symlink to `code/bin/scrna-pipeline` |
| Workshop workspace | `…/demo/` (group-writable `results/demux`, `logs/demux`; read-only `data/Pool_ctr/cellsnp`; `fastq/` symlinks to dev-repo FASTQs) |
| Pooled scRNA input | `data/Pool_ctr/` in dev repo: 12 GB CellRanger BAM, barcodes, pre-computed `cellsnp/` (copied in 2026-05-27) |
| Andras's own conda env | `/home/nemeth/miniforge3/envs/scrna` (home dir is private — colleagues cannot use it) |
| Old local DB (stale since 2026-11-11) | `results/variants.db` in dev repo |
| Old local references (duplicate, deletable) | `data/reference/` in dev repo (~30 GB) |
| Tutorial page source | `docs/workshop.html` (published as a claude.ai artifact under the old account) |

Server facts: 4 CPU cores, ~1 TB RAM but heavily shared (often only 50–200 GB
free), `/groups` is **NFS4**, ~30 users logged in. Docker is installed but the
user is **not** in the `docker` group; no Apptainer.

---

## 3. Branches and versions

Convention: one branch per version, merged into `main` at the end; old version
branches are never deleted.

| Branch | Content |
|---|---|
| `v1` | Genome-wide calling, ALT-only, batch runs, Grist push |
| `v2` | Genotyping at common SNPs, 0/0 stored, duplicate run_id guard, external VCF validation, binomial demux, Vireo (unguided) + DB matching, merge_demux, 1:1 Hungarian matching (`--unique`), unmatched-donor insertion, 0/0 backfill |
| `v3` | Reference-guided Vireo, WES mode (BWA), .fastq.gz support, temp intermediates, split analysis notebooks 01–04, thesis figures notebook 05 |
| `v4` (current) | Packaging / lab-wide rollout (see §4). Branched from main after main was fast-forwarded to v3. |

State on 2026-11-17: `main` = `v3` = `5c12a68`; `v4` = `923a39b` (6 commits
ahead). **Nothing since `27c3e17` has been pushed** — the session had no GitHub
credentials. Andras must run:

```
git push origin main v3 && git push -u origin v4
```

v4 has **not** been merged into main yet. The shared code copy tracks v4
directly; once v4 is merged, switch it: `git -C /groups/nils/data/scrna_pipeline/code checkout main`.

---

## 4. Timeline of this session (Nov 2026)

### 4.1 Thesis presentation work (notebook 05, figures)
- Fixed a `TypeError` in notebook 05 (fig22: `correct` vs `correct_arr`).
- Added same-donor box plots (orange #E76F51) to fig3, fig3b, fig13, fig14;
  panel titles on fig3b; moved fig21 labels; fig22 gap legend.
- Checked concordance with a PASS-only filter: **identical** results (backfilled
  0/0 rows are always PASS).
- Produced `figures/presentation_v6/` (fig_why_concordance, fig_thresholds,
  fig_loo_identify; PNG 200 dpi, white, no in-image titles) and
  `presentation_numbers_v6.md`.
- `figures/` is now **git-ignored** (personal presentation output; files kept on
  disk). Caveat: v1/v2 still track figures, so checking them out and back
  deletes the figure files from disk — back up `figures/` first.

Key thesis numbers (44 sc/bulk runs, same-donor groups SW480/SW620,
IS1/IS2/IS3, COLO201/COLO205, HCT8/HCT15/DLD1):
- Concordance (binary ALT-presence, key chrom:pos): same-donor min 0.943,
  different-donor max 0.829, margin +0.114, separation 4.36 IQR units, 0 overlaps.
  Separating midpoint cutoff 0.886 (an analysis threshold, **not** a tool rule).
- Jaccard (PASS, chrom:pos:ref:alt): margin +0.019 (0.42 IQR), 0 overlaps.
- Containment: margin −0.204 (−1.74 IQR), 1,274 overlapping pairs — fails.
- LOO (concordance-based, not match_vcf): SW620 → SW480 0.963, runner-up LIM2099 0.802.
- `match_vcf.py` supports only overlap(=containment, default, `--min-match 40`) and jaccard; **no concordance**.
- pool_ctr_v4 re-run with `--min-concordance 0.80 --min-gap 0.15 --min-positions 200 --unique`:
  same assignments as the original 0.65/0.10 run; only HCT116 drops high → low (gap 0.119).
- **Mapping gate cannot fail for STAR runs**: STAR runs without
  `--outSAMunmapped Within`, so flagstat always reports 100 %. STAR unique rates:
  68.86–92.98 %; at 85 % SNUC2B (68.86) and VACO4S (84.96) would fail. Not fixed
  (pipeline defaults are only changed with Andras's agreement).
- Runtimes (measured): bulk RNA run median 13.8 min (6.4–92, n=40; the 92-min
  run included building the STAR index ~41 min); WES 3.1–5.4 h (median 3.6 h);
  Vireo unguided on Pool_ctr ~64–79 min on 1 core; match_vireo ~30 s.
- **Corrections made later** to `presentation_numbers_v6.md`: cellsnp-lite was
  never run by this pipeline (pileup was pre-computed), Pool_ctr has
  **10,308 cells** (not 299), bulk runtimes recomputed, commands fixed.

### 4.2 Repo housekeeping
- main fast-forwarded to v3; v4 created from main.

### 4.3 v4 step 1 — shared reference + shared database (`cdf2b39`)
- `scripts/setup_shared.sh` creates the shared layout (references read-only via
  `chmod a-w`; DB dir 2775; DB and lock 664). Ran it: 40 GB copied.
- `scripts/shared_paths.py`: defaults for all paths; overridable with
  `SCRNA_SHARED_ROOT` (whole tree) or `SCRNA_DB` (DB only).
- `scripts/db_lock.py`: `write_lock(db)` — POSIX `lockf` on `<db>.lock`
  (SQLite locking is unreliable on NFS). Used by every writer: load_to_db,
  demux `--load-db`, add_unmatched_donors, backfill_panel_genotypes, backup_db.
  Tested on NFS (second writer waits, prints holder).
- DB journal mode is `delete` (never switch to WAL on NFS).
- `runs` gained `added_by` and `pipeline_version` (git describe, or
  `SCRNA_PIPELINE_VERSION`). Existing 54 runs backfilled as `nemeth`/`pre-v4`.
- `scripts/backup_db.py`: SQLite online backup under the lock, keeps 14.
  First backup taken. **Cron job not installed** (Andras to add).
- Grist push now skipped when `GRIST_*` env vars are absent (colleagues have no `.env`).
- Fixed: fresh DB crashed on the `samples` migration (ran before CREATE TABLE).
- `run_pipeline.py` passes `--rerun-triggers mtime` (code changes must not
  trigger rebuilds of read-only shared indexes).
- Notebooks 01–04 and DEMO.md point at the shared DB; notebook 05 still reads the
  old local `results/variants.db` (intentionally, it's the finished thesis notebook).

### 4.4 v4 step 2 — packaging (`68c60e8`)
- `environment.yml` rewritten: previously listed only 5 tools; now pins
  everything (snakemake 9.17.3, STAR 2.7.10b, bwa 0.7.19, samtools/bcftools/htslib 1.23.1,
  cellsnp-lite 1.2.3, numpy/scipy/pandas, pip: cyvcf2 0.32.1, vireoSNP 0.5.9,
  grist-api, python-dotenv, jupyterlab, plotly, kaleido). Channels conda-forge +
  bioconda only (dropped `defaults` for licensing). Verified by building a fresh env.
- Pipeline runs from **any working directory**: code paths resolve from the repo
  (`workflow.basedir` / `REPO_DIR`), `data/ results/ logs/` stay relative to cwd.
  Removed the hardcoded `/home/nemeth/...` conda PATH hack from the Snakefile.
- `run_pipeline.py` prepends the interpreter's bin dir to PATH and fails early
  with a clear message if tools are missing; batch mode honours `--foreground`.
- `Dockerfile` + `.dockerignore` written (micromamba base, mounts `/shared` and
  `/work`, entrypoint adds `--foreground`). **Never built or tested** (no Docker access).
- README Setup rewritten (shared install / own clone / container / outside the lab,
  shared database section).

### 4.5 v4 step 3 — lab-wide install + workshop (`bf33bd4`, `db90f5a`, `04c1075`, `923a39b`)
- `bin/scrna-pipeline` launcher; shared `code/`, `env/`, `bin/` set up.
- **demux.py scoring rewrite**: `score_cells` used dense matrices (26 GB, 10+ min
  for 10k cells × 54 refs). Now computes on non-zero entries only: identical
  scores (max diff 6e-12) and identical assignments, ~20× faster, 3.4 GB peak.
- `sqlite` added to environment.yml (CLI for DB demo).
- `ancient()` on inputs of normalize_snp_panel, star_index, bwa_index: a git
  clone/pull reset timestamps and Snakemake tried to rebuild the shared SNP panel.
- Workshop: notebook `notebooks/06_demo_compare.ipynb`, tutorial page
  `docs/workshop.html`, four pre-computed runs in `demo/results/demux/precomputed_*`.

---

## 5. The workshop (ready to run)

Plan: explain pipeline → dry-run (bulk DLD1 `SRR5071662`, WES `B2006`) → DB
queries (`sqlite3 -readonly`) → four groups demultiplex Pool_ctr → compare in
notebook 06 → discussion. Tutorial page has all copy-paste commands.

| Group | Run id | Arguments |
|---|---|---|
| 1 | `group1_guided7` | `--known-lines RKO HCT116 HT29 HCT8 CACO2 DLD1 LIM1215` |
| 2 | `group2_unguided` | `--n-donors 7 --unique` (≈80 min — launch first; compare with pre-computed) |
| 3 | `group3_no_hct8` | `--known-lines RKO HCT116 HT29 CACO2 DLD1 LIM1215` |
| 4 | `group4_no_dld1` | `--known-lines RKO HCT116 HT29 HCT8 CACO2 LIM1215` |

All from `/groups/nils/data/scrna_pipeline/demo` with
`scrna-pipeline --mode scrna --cellsnp-dir data/Pool_ctr/cellsnp --cores 1 --demux-run-id <id> …`.
One person per group launches (same id twice collides on file ownership).

Pre-computed results (DLD1 and HCT8 are from the same donor):

| Run | Assigned to a line | DLD1 | HCT8 | HCT8+DLD1 doublet | Unassigned | Time |
|---|---|---|---|---|---|---|
| Guided 7 | 8,338 | 33 | 3 | 198 | 1,269 | 66 s |
| Unguided | 6,508 | 0 | 0 | 0 | 662 (+1,803 unmatched donors) | 79 min, 21 GB |
| No HCT8 | 8,649 | 342 | – | 0 | 1,119 | 67 s |
| No DLD1 | 8,636 | – | 340 | 0 | 1,123 | 66 s |

Unguided matches RKO (0.90, high), HCT116 (0.87, high), HT29 (0.83, low);
four donors unmatched. Guided-7 reproduces the earlier `pool_ctr_v5` exactly.
Vireo `-p 2` cuts unguided to 43 min with identical output but peaks at 35 GB —
decided **not** to use it for the workshop.

DB demo queries (all < 1 s): pool lines in `runs`; TP53 rs1042522
(chr17:7676154) genotypes across the 7 lines (DLD1 = HCT8 = 1/1, CACO2 no call).
Cross-run joins on `genotype_calls` are **slow (> 2 min)** because there is no
index on `variant_id` — avoid them live, or discuss adding an index first.

---

## 6. Decisions and preferences to respect

- **Never change pipeline defaults silently** — Andras asked to be told first.
- Commit per logical step on the version branch; only commit/push when asked.
  Commit messages explain *why*.
- One branch per version; never delete old branches.
- Shared references stay read-only; shared DB is written only under the lock.
- Lab colleagues use the shared install (conda route); Docker only for other
  machines. If containers are needed on the server, ask IT for **Apptainer**
  rather than docker-group membership (root-equivalent).
- Presentation figures: PNG, white background, 200 dpi, no titles in the image.
- Andras likes copy-paste-ready terminal commands and concise recaps.
- Be careful with numbers: several earlier figures were wrong because they were
  estimated instead of measured (see §4.1 corrections). Measure, cite the source.

---

## 7. Open items / known issues

1. **Push** main, v3, v4 to GitHub (credentials needed on Andras's side).
2. **Merge v4 into main** when the rollout is accepted; then point the shared code copy at main.
3. **Install the backup cron job**:
   `0 3 * * * cd /groups/nils/data/scrna_pipeline/code && /groups/nils/data/scrna_pipeline/env/bin/python scripts/backup_db.py`
4. **Mapping-rate gate is ineffective for STAR** (flagstat always 100 %). Options:
   gate on STAR's unique-mapping % from `Log.final.out`, or add `--outSAMunmapped Within`.
   Needs Andras's decision (changes pipeline behaviour).
5. **Threshold inconsistencies**: `config.yaml` demux `min_concordance 0.55`, `min_gap 0.15`
   vs `run_pipeline.py` defaults `0.80` / `0.10`; Snakemake `match_vireo_donors` rule has its own defaults.
6. **No index on `genotype_calls(variant_id)`** → cross-run SQL is slow. Adding one is a
   shared-DB schema change; ask first.
7. **Dockerfile untested.**
8. `load_to_db.py` inserts the `runs` row and commits before loading genotypes — a crash
   mid-load leaves a partial run (now at least serialised by the lock).
9. scrna `--dry-run` still creates empty `results/demux/<id>` and `logs/demux` dirs.
10. Old local copies: `results/variants.db` (stale) and `data/reference/` (30 GB) in the dev
    repo can be renamed/deleted once Andras is happy with the shared copies.
11. Earlier idea (May 2026), not started: STARsolo rule so pooled scRNA FASTQs can be
    processed without CellRanger. Also: cellsnp-lite has never actually been run by the
    pipeline end-to-end (Pool_ctr pileup was pre-computed).
12. Grist credentials live in Andras's private `.env` (mode 600); never commit.
13. The tutorial page was published as a claude.ai artifact under the old account, which the
    new account cannot open. Republish `docs/workshop.html` from the new account. Sources of
    **all** old-account artifacts (this project and others) are in
    `/home/nemeth/claude_account_export/` (see its README), including a rendered transcript of
    an earlier scrna_project session (May–June 2026) in `artifacts/transcript.md.html`.

---

## 8. Useful commands

```bash
# Use the shared install (any project folder)
export PATH=/groups/nils/data/scrna_pipeline/bin:/groups/nils/data/scrna_pipeline/env/bin:$PATH
scrna-pipeline --mode bulk --sample <ID> --run-id <ID>_hg38 --fastq-dir <dir> --dry-run

# Deploy committed changes to the shared copy
git -C /groups/nils/data/scrna_pipeline/code pull

# Update shared env after environment.yml changes
mamba env update -p /groups/nils/data/scrna_pipeline/env -f environment.yml

# Query the DB safely
sqlite3 -readonly -header -column /groups/nils/data/scrna_pipeline/db/variants.db "SELECT COUNT(*) FROM runs;"

# Re-run a pre-computed workshop configuration
cd /groups/nils/data/scrna_pipeline/demo
scrna-pipeline --mode scrna --cellsnp-dir data/Pool_ctr/cellsnp --cores 1 --foreground \
    --demux-run-id precomputed_guided7 --known-lines RKO HCT116 HT29 HCT8 CACO2 DLD1 LIM1215
```
