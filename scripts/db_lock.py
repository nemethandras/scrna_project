"""Exclusive write lock for the shared SQLite database.

SQLite's own file locking is unreliable on NFS, so every writer also takes a
POSIX lock on a sidecar file (<db>.lock), which NFSv4 enforces server-side.
"""
import contextlib
import fcntl
import getpass
import os
import socket
import sys
import time


@contextlib.contextmanager
def write_lock(db_path):
    fd = os.open(db_path + ".lock", os.O_RDWR | os.O_CREAT, 0o664)
    try:
        try:
            fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            holder = os.pread(fd, 200, 0).decode(errors="replace").strip() or "another process"
            print(f"  Database locked by {holder} — waiting...", file=sys.stderr, flush=True)
            fcntl.lockf(fd, fcntl.LOCK_EX)
        with contextlib.suppress(PermissionError):
            os.fchmod(fd, 0o664)
        os.ftruncate(fd, 0)
        os.pwrite(fd, (f"{getpass.getuser()}@{socket.gethostname()} pid {os.getpid()} "
                       f"since {time.strftime('%Y-%m-%d %H:%M:%S')}\n").encode(), 0)
        yield
    finally:
        os.close(fd)
