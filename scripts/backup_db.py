#!/usr/bin/env python3
"""Snapshot the shared genotype database and prune old snapshots.

Intended for a nightly cron job, e.g.
    0 3 * * *  cd /path/to/scrna_project && python scripts/backup_db.py
"""
import argparse
import sqlite3
import time
from pathlib import Path

from db_lock import write_lock
from shared_paths import DEFAULT_DB


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", default=DEFAULT_DB)
    p.add_argument("--backup-dir", default=None,
                   help="default: <db dir>/backups")
    p.add_argument("--keep", type=int, default=14, help="snapshots to keep (default: 14)")
    args = p.parse_args()

    db = Path(args.db)
    backup_dir = Path(args.backup_dir) if args.backup_dir else db.parent / "backups"
    backup_dir.mkdir(parents=True, exist_ok=True)
    target = backup_dir / f"{db.stem}_{time.strftime('%Y%m%d_%H%M%S')}.db"

    with write_lock(str(db)):
        src = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        dst = sqlite3.connect(target)
        src.backup(dst)
        dst.close()
        src.close()
    target.chmod(0o444)
    print(f"Backed up {db} -> {target}")

    snapshots = sorted(backup_dir.glob(f"{db.stem}_*.db"))
    for old in snapshots[:-args.keep]:
        old.unlink()
        print(f"Removed old snapshot {old.name}")


if __name__ == "__main__":
    main()
