#!/usr/bin/env python3
"""Compatibility wrapper for the test-owned mobile main figure runner."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


if __name__ == "__main__":
    tests_dir = Path(__file__).resolve().parents[1] / "tests"
    sys.path.insert(0, str(tests_dir))
    runpy.run_path(
        str(tests_dir / "main_figure_mobile.py"),
        run_name="__main__",
    )
