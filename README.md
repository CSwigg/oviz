# Oviz

Oviz turns three-dimensional stellar data into interactive web figures that can
be explored in space and time. Its current focus is Gaia data for young star
clusters and associations, together with emerging three-dimensional maps of the
interstellar medium (ISM).

Oviz combines:

- [Three.js](https://threejs.org/) for the interactive 3D scene;
- [Aladin Lite](https://aladin.cds.unistra.fr/AladinLite/doc/) for registered
  all-sky HiPS imagery;
- [galpy](https://docs.galpy.org/en/latest/) for Galactic orbit integration;
- Astropy, NumPy, and pandas for astronomical coordinates and tabular data.

Together these tools let a viewer move between a Galactic 3D scene and the sky,
follow clusters through time, compare stars with dust or emission maps, and save
a scientific narrative as an ordered sequence of viewer States.

## Scientific scope

Oviz is designed around data products such as:

- [Gaia Data Release 3](https://doi.org/10.1051/0004-6361/202243940);
- the [Gaia DR3 all-sky open-cluster catalogue](https://doi.org/10.1051/0004-6361/202346285);
- the [SigMA view of the Sco-Cen association](https://doi.org/10.1051/0004-6361/202243690);
- the [Edenhofer et al. parsec-scale 3D dust map](https://doi.org/10.1051/0004-6361/202347628);
- the [Vergely et al. Gaia-2MASS 3D dust maps](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/664/A174?format=html&tex=true).

The data model can also represent other point, line, image, and volume layers.
Broader astronomical views and analysis tools are planned.

## Install

```bash
python -m pip install -e .
```

Oviz requires Python 3.8 or newer. Its core Python dependencies are NumPy,
pandas, Astropy, and galpy.

## Minimal example

A time-varying `Trace` needs Cartesian Galactic positions in parsecs, velocities
in km/s, a name, and an age:

```python
import numpy as np
import pandas as pd

from oviz import Animate3D, Trace, TraceCollection

clusters = pd.DataFrame(
    {
        "x": [35.0, -80.0],
        "y": [120.0, 60.0],
        "z": [15.0, -25.0],
        "U": [-11.0, -9.5],
        "V": [-18.0, -21.0],
        "W": [-7.0, -5.0],
        "name": ["Cluster A", "Cluster B"],
        "age_myr": [12.0, 28.0],
        "n_stars": [180, 75],
    }
)

trace = Trace(
    clusters,
    data_name="Young clusters",
    color="#e34a4a",
    size_by_n_stars=True,
)
scene = Animate3D(
    TraceCollection([trace]),
    xyz_widths=(2000, 2000, 600),
    figure_theme="dark",
)

figure = scene.make_plot(
    time=np.arange(-30.0, 0.5, 0.5),
    galactic_mode=True,
    enable_sky_panel=True,
    renderer="threejs",
    show=False,
    compress_scene_spec=True,
)
figure.write_html("young_clusters.html")
```

Use `Layer(..., assume_stationary=True)` for a static XYZ catalogue. Volume
layers and optional cluster-member catalogues can be passed to
`Animate3D.make_plot()`; see the [Python API](docs/source/python_api.rst).

## States and sharing

A State is a complete viewer snapshot: camera, time, 3D/Sky mode, trace and
volume settings, Aladin layers, selections, panels, and presentation settings.
Open **States**, arrange the views you want to explain, and export either an
editable States HTML or a present-only HTML with forward/back controls.

The result is one HTML file containing the scientific scene and saved States.
It can be opened locally, served from any static web host, or sent directly to a
collaborator. This makes it practical to link someone to a guided view of a Gaia
sample without requiring them to install Python or reproduce the analysis.

Three.js and Aladin Lite are loaded from public CDNs, and HiPS backgrounds load
from their survey servers, so those features require an internet connection.

States can also be controlled from JavaScript:

```javascript
const root = document.querySelector('[id^="oviz-three-"]');
const viewer = window.Oviz.get(root.id);
await viewer.states.goTo(1);       // one-based State index
await viewer.states.next();
const html = await viewer.states.exportHtml({ download: false });
```

See the [browser API reference](docs/source/browser_api.rst) for navigation,
authoring, events, and `postMessage` commands.

## Input conventions

- Positions `x`, `y`, and `z` are parsecs in Galactic Cartesian coordinates.
- Velocities `U`, `V`, and `W` are km/s.
- `age_myr` is in Myr; `n_stars` is optional unless size scaling uses it.
- Every time array must include `t = 0`.
- A member-star table must include a cluster-name column and either Galactic
  (`l`, `b`) or ICRS (`ra`, `dec`) sky coordinates.

## Documentation and tests

- [Overview](docs/source/overview.rst)
- [Python API](docs/source/python_api.rst)
- [Browser API](docs/source/browser_api.rst)
- [Agent guide](AGENTS.md)

Run the maintained test suite with:

```bash
pytest -q tests
```
