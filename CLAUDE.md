# CLAUDE.md

Snakemake pipeline that fingerprints cell lines at ~7M common SNPs (bulk RNA via
STAR, WES via BWA) into a shared SQLite DB, and demultiplexes pooled scRNA-seq
against it (binomial scorer + Vireo, guided or unguided). Owner: Andras (Nils lab,
server `menten`). Full history, decisions and open items: `docs/PROJECT_HISTORY.md`
— read it before making changes.

## Current state
- Active branch **v4** (lab-wide rollout). `main` = `v3`. One branch per version; never delete old ones.
- Unpushed commits on main/v3/v4 may exist — check `git status -sb` / `git log origin/v4..v4`.
- Shared install: `/groups/nils/data/scrna_pipeline/` — `reference/hg38` (read-only),
  `db/variants.db` (shared, lock-protected), `code/` (clone of this repo tracking v4),
  `env/` (conda env), `bin/scrna-pipeline` (launcher), `demo/` (workshop workspace).
- After committing, deploy with `git -C /groups/nils/data/scrna_pipeline/code pull`.

## Rules
- Do not change pipeline defaults or thresholds without asking Andras first.
- Commit only when asked; commit per logical step; messages explain why.
- Any DB write must go through `scripts/db_lock.write_lock`. Keep SQLite journal mode `delete` (NFS).
- Never write into the shared reference directory; reference-building rules use `ancient()` inputs.
- Paths: code resolves from the repo; `data/ results/ logs/` are relative to the working directory.
  Defaults come from `scripts/shared_paths.py` (`SCRNA_SHARED_ROOT`, `SCRNA_DB`).
- `.env` holds private Grist credentials — never read into commits or print.
- Report measured numbers with their source; don't estimate runtimes or counts.

## Environment
- Conda env: `/groups/nils/data/scrna_pipeline/env` (shared) or `/home/nemeth/miniforge3/envs/scrna` (Andras's).
- Server: 4 cores, shared RAM, `/groups` is NFS. Docker not usable (no group membership).
