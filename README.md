# oviz

`oviz` is a Python package for 3D orbit visualization of stellar clusters and associations.

It combines:
- orbital integration (`galpy`)
- animated 3D Plotly figures
- an optional Dash app with an on-sky Aladin Lite panel synced to cluster selection

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
- `plotly`
- `dash`
- `astropy`
- `galpy`

## Core Data Model

`Trace` requires a dataframe with these columns:
- `x`, `y`, `z` (pc)
- `U`, `V`, `W` (km/s)
- `name`
- `age_myr`

Optional:
- `n_stars` (required only if `size_by_n_stars=True`)

## Quick Start

```python
import numpy as np
import pandas as pd

from oviz import Trace, TraceCollection, Animate3D
from oviz.app import run_dash_app_in_notebook

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

trace = Trace(
    cluster_df,
    data_name="Cluster A",
    color="cyan",
    marker_style="circle",
)
collection = TraceCollection([trace])

viz = Animate3D(
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
)

app = run_dash_app_in_notebook(
    figure=fig,
    mode="external",
    port=8061,
    enable_age_filter=False,
)
```

## Dash Sky Panel (Aladin Lite)

Enable the sky panel when launching Dash:

```python
app = run_dash_app_in_notebook(
    figure=fig,
    mode="external",
    port=8061,
    title="oviz sky panel",
    enable_age_filter=False,
    enable_sky_panel=True,
    sky_radius_deg=7.0,
    sky_frame="galactic",
    sky_survey="P/DSS2/color",
    cluster_members_file="/absolute/path/to/members.csv",  # optional
)
```

Sky panel behavior:
- click a cluster member in 3D at `t = 0` to set the footprint and sky target
- cone footprint is shown only at `t = 0`
- fullscreen, restore, and hide/show controls are built in
- panel can be dragged in normal mode using the top drag strip

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
- `create_dash_app`
- `run_dash_app`
- `run_dash_app_in_notebook`
- `launch_from_animate3d`

## Testing

Run the sky-panel tests with:

```bash
pytest -q tests/test_sky_panel.py
```

## Notes

- `sky_frame="icrs"` is currently not supported in the Dash sky panel path.
- The Dash sky panel uses Aladin Lite loaded from CDN, so network access is required at runtime.
