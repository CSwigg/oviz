"""The Oviz viewer: a fast, self-contained WebGL2 figure runtime.

``compile_scene_spec`` turns a scene spec into a columnar binary bundle and
``OvizFigure`` writes it into one standalone HTML file.
"""

from .bundle import Bundle, BundleBuilder
from .compile import CompileError, compile_scene_spec

__all__ = ["Bundle", "BundleBuilder", "CompileError", "compile_scene_spec"]
