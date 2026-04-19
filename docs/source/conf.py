from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


project = "oviz"
copyright = "2024, Cameren Swiggum"
author = "Cameren Swiggum"
release = "0.1.0"


extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = []


html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
