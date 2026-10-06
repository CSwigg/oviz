"""The Oviz viewer: a fast, self-contained WebGL2 figure runtime.

``compile_scene_spec`` turns a scene spec into a columnar binary bundle and
``OvizFigure`` writes it into one standalone HTML file.
"""

from .bundle import Bundle, BundleBuilder
from .compile import CompileError, compile_scene_spec
from .figure import OvizFigure, render_bundle_html
from .upgrade import read_legacy_scene_spec, upgrade_html

__all__ = [
    "Bundle",
    "BundleBuilder",
    "CompileError",
    "OvizFigure",
    "compile_scene_spec",
    "read_legacy_scene_spec",
    "render_bundle_html",
    "upgrade_html",
]
