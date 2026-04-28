# __init__.py
import importlib.util

from .scene import Scene3D
from .threejs_profiles import (
	build_threejs_profile,
	galactic_lite_profile,
	galactic_simple_profile,
	lite_profile,
	merge_threejs_profile,
	normalize_threejs_initial_state,
	threejs_profile,
	website_background_profile,
)
from .traces import Layer, LayerCollection, Trace, TraceCollection
from .viz import Animate3D

if importlib.util.find_spec("dash") is None:  # pragma: no cover - optional dash dependency path
	create_dash_app = None
	run_dash_app = None
	run_dash_app_in_notebook = None
	launch_from_animate3d = None
else:
	from .app import (
		create_dash_app,
		run_dash_app,
		run_dash_app_in_notebook,
		launch_from_animate3d,
	)

from . import orbit_maker
from . import point_sizes

__all__ = [
	"Animate3D",
	"Scene3D",
	"Trace",
	"TraceCollection",
	"Layer",
	"LayerCollection",
	"threejs_profile",
	"lite_profile",
	"galactic_lite_profile",
	"galactic_simple_profile",
	"website_background_profile",
	"build_threejs_profile",
	"merge_threejs_profile",
	"normalize_threejs_initial_state",
	"create_dash_app",
	"run_dash_app",
	"run_dash_app_in_notebook",
	"launch_from_animate3d",
]
