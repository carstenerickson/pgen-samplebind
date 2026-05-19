"""Advisory output-prefix lock + filesystem warnings.

`output_lock(prefix)` takes a non-blocking `fcntl.flock` on `{prefix}.lock`.
Released on context exit. Used by the `merge` orchestrator so two concurrent
invocations writing to the same prefix can't silently corrupt each other's
outputs.

`detect_network_filesystem(path)` is a best-effort check warning the user when
the lock file lives on NFS/SMB/CIFS — `fcntl.flock` over network filesystems
is documented as unreliable (no-op on some implementations), so the lock
becomes advisory-only-on-paper.
"""

from __future__ import annotations

import contextlib
import fcntl
import os
import platform
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from .errors import IOFailure

_NETWORK_FS_TYPES = frozenset(
    {
        "nfs",
        "nfs4",
        "smbfs",
        "smb",
        "smb2",
        "smb3",
        "cifs",
        "fuse.sshfs",
        "afpfs",
    }
)


def detect_network_filesystem(path: Path) -> str | None:
    """Best-effort detection of NFS/SMB/CIFS at `path`. Returns the fs type
    string when recognized as networked, else None.

    Linux: parse `/proc/self/mountinfo` for the longest matching mount.
    macOS / other: `os.statvfs` doesn't expose fs type; we fall back to None
    (no warning) rather than emit false positives. The advisory lock is still
    taken; the warning is purely diagnostic.
    """
    try:
        resolved = path.resolve().parent
    except OSError:
        return None

    if platform.system() != "Linux":
        return None

    mountinfo = Path("/proc/self/mountinfo")
    try:
        lines = mountinfo.read_text().splitlines()
    except OSError:
        return None

    best_match = ("", "")  # (mountpoint, fstype)
    resolved_str = str(resolved)
    for line in lines:
        # mountinfo format:
        #  ID parent_id maj:min root mountpoint mount_opts - fstype source super_opts
        parts = line.split(" ")
        try:
            sep = parts.index("-")
        except ValueError:
            continue
        if sep + 1 >= len(parts):
            continue
        mountpoint = parts[4]
        fstype = parts[sep + 1]
        if (
            resolved_str == mountpoint or resolved_str.startswith(mountpoint.rstrip("/") + "/")
        ) and len(mountpoint) > len(best_match[0]):
            best_match = (mountpoint, fstype)

    if best_match[1].lower() in _NETWORK_FS_TYPES:
        return best_match[1]
    return None


@contextmanager
def output_lock(prefix: Path) -> Iterator[None]:
    """Take a non-blocking exclusive `fcntl.flock` on `{prefix}.lock`.

    Released on context exit. The lock file is created if absent and left in
    place after release (next acquirer just re-locks the existing file). A
    networked-filesystem warning is emitted to stderr at acquire time when
    the lock file lives on NFS/SMB/CIFS, since flock semantics on those are
    implementation-defined and effectively no-op.

    Raises:
        IOFailure: another process already holds the lock for this prefix.
    """
    lock_path = Path(str(prefix) + ".lock")

    fs_type = detect_network_filesystem(lock_path)
    if fs_type is not None:
        sys.stderr.write(
            f"warning: output prefix is on {fs_type} — fcntl.flock may be a no-op; "
            "concurrent writes to the same prefix are not safely detectable.\n"
        )

    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise IOFailure(
                f"output prefix {prefix} is locked by another pgen-samplebind "
                f"process (lock file: {lock_path}). Wait for it to finish or "
                "remove the lock file if the holding process has died."
            ) from e
        try:
            yield
        finally:
            with contextlib.suppress(OSError):
                fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)
