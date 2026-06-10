# oviz

`oviz` is a Python package for 3D astronomical spatial visualization, with orbit tracing as a core workflow.

It combines:
- orbital integration (`galpy`)
- a standalone Three.js viewer for interactive 3D exploration
- on-sky Aladin Lite backgrounds and sky controls inside the Three.js viewer

## Installation

```bash
python -m pip install -e .
```

or

```bash
python -m pip install .
```

Main dependencies:
- `numpy`
- `pandas`
- `astropy`
- `galpy`

The supported renderer is the standalone Three.js export path.

## Core Data Model

`Trace` requires a dataframe with these columns:
- `x`, `y`, `z` (pc)
- `U`, `V`, `W` (km/s)
- `name`
- `age_myr`

Optional:
- `n_stars` (required only if `size_by_n_stars=True`)

For a more general astronomy-facing API, `oviz` also exposes:
- `Layer`
- `LayerCollection`
- `Scene3D`

`Layer` preserves the same rendering/orbit behavior as `Trace`, but it can also represent static spatial layers when you pass `assume_stationary=True`.

## Quick Start

```python
import numpy as np
import pandas as pd

from oviz import Layer, LayerCollection, Scene3D, build_threejs_profile

# Example minimal cluster dataframe
cluster_df = pd.DataFrame(
    {
        "x": [10.0, 12.0],
        "y": [20.0, 19.0],
        "z": [5.0, 6.0],
        "U": [-12.0, -11.5],
        "V": [220.0, 221.0],
        "W": [7.0, 7.2],
        "name": ["Cluster A", "Cluster A"],
        "age_myr": [8.0, 8.0],
    }
)

layer = Layer.from_dataframe(
    cluster_df,
    layer_name="Cluster A",
    color="cyan",
    marker_style="circle",
)
collection = LayerCollection([layer])

viz = Scene3D(
    data_collection=collection,
    xyz_widths=(1000, 1000, 400),
    figure_theme="dark",
)

time = np.round(np.arange(0, -20.5, -0.5), 1)
fig = viz.make_plot(
    time=time,
    show=False,
    galactic_mode=True,
    show_galactic_guides=False,
    threejs_initial_state=build_threejs_profile("full"),
)
fig.write_html("oviz_scene.html")
```

## Sky View

Enable sky features directly on the Three.js viewer:

```python
fig = viz.make_plot(
    time=time,
    show=False,
    enable_sky_panel=True,
    sky_radius_deg=7.0,
    sky_frame="galactic",
    sky_survey="P/DSS2/color",
    cluster_members_file="/absolute/path/to/members.csv",  # optional
)
```

Sky behavior:
- click a cluster member in 3D at `t = 0` to set the footprint and sky target
- cone footprint is shown only at `t = 0`
- sky layers, Aladin backgrounds, and spectrum aperture controls are built into the viewer

## Optional Member File Format

If `cluster_members_file` is passed, it should be a CSV with:
- `name` (cluster label)
- `l` (galactic longitude, deg)
- `b` (galactic latitude, deg)

When a cluster is selected, only that cluster's member stars are rendered in Aladin.

## Public Entry Points

From `oviz`:
- `Trace`
- `TraceCollection`
- `Animate3D`
- `Layer`
- `LayerCollection`
- `Scene3D`
- `build_threejs_profile`
- `threejs_profile`
- `lite_profile`
- `galactic_lite_profile`
- `website_background_profile`

## Testing

Run the test suite with:

```bash
pytest -q tests
```

## Notes

- Live Aladin sky backgrounds load remote Aladin/HiPS assets, so network access is required at runtime for those layers.
