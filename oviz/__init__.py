# __init__.py
import importlib.util

from .traces import Trace, TraceCollection
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

