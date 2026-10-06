"""Provenance helpers for the figure scripts.

Published figures record their input files by name and a short SHA-256,
never by local path.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


def short_sha256(path: Path | str | None, length: int = 12) -> str | None:
    """The first ``length`` hex digits of a file's SHA-256, or None if it is missing."""
    if path is None:
        return None
    path = Path(path).expanduser()
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()[:length]


def hash_provenance_paths(scene: dict) -> dict:
    """Replace the local paths older scenes recorded in ``provenance`` with hashes.

    Each ``path`` / ``*_path`` entry becomes ``sha256`` / ``*_sha256`` in place
    (None when the file is not on this machine), next to the file name the
    entry already records.
    """

    def is_path(key: str) -> bool:
        return key == "path" or key.endswith("_path")

    provenance = scene.get("provenance")
    if not isinstance(provenance, dict):
        return scene
    for name, entry in provenance.items():
        if isinstance(entry, dict) and any(is_path(key) for key in entry):
            provenance[name] = {
                (key[:-4] + "sha256" if is_path(key) else key): (short_sha256(value) if is_path(key) else value)
                for key, value in entry.items()
            }
    return scene
