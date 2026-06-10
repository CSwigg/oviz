#!/usr/bin/env python3
"""Compatibility wrapper for the test-owned young-stars figure runner."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "tests" / "young_stars.py"),
        run_name="__main__",
    )
