---
name: oviz
description: Build, extend and publish interactive 3D HTML figures of Milky Way data with Oviz. Covers star clusters, associations or stars traced back in time with galpy orbits (optionally in a spiral-arm potential), 3D dust and emission volumes, member stars on a sky view registered to Aladin Lite, published spiral arms, and saved views that play as a presentation. Use it whenever someone wants to visualize Gaia clusters, stellar orbits or tracebacks, 3D dust maps, or Galactic structure interactively, or wants to make, update, debug or share an Oviz figure or presentation, even if they don't say "Oviz".
---

# Oviz

Oviz turns tables of Galactic objects into one self-contained HTML file. Each
object's orbit is integrated with galpy over a time grid, and the result is
written as a WebGL2 viewer. In the figure, readers scrub time, orbit in 3D,
switch to the sky, inspect objects, lasso, and present saved views. Nobody
needs Python to open it.

The worked reference is `tests/main_figure_oct1.py` (the "October 1 figure",
live at https://cswigg.github.io/cam_website/oviz_figures/oviz_oct1.html). It
shows:
- young clusters integrated in a spiral potential;
- 3D dust and Hα volumes;
- Sky member stars;
- two published sets of spiral arms and two pulsar traces, all hidden until
  switched on.

Read it when you need a full-scale example. This file distils how it works.

## Setup

```bash
pip install git+https://github.com/CSwigg/oviz.git   # or: pip install -e . in a clone
```

Oviz needs Python ≥ 3.10; pip brings in NumPy, pandas, Astropy, galpy and
SciPy. Checking a figure needs a browser (Chrome works headless). Sky imagery
needs network access; the 3D scene does not.

## 1. Prepare the data

Use a pandas DataFrame per group of objects, with these units and frames.
Getting them wrong silently produces wrong orbits.

| column | meaning |
| --- | --- |
| `x, y, z` | heliocentric Galactic Cartesian position, pc. x points to the Galactic centre, y to l = 90°, z to the north Galactic pole |
| `U, V, W` | heliocentric velocity, km/s, same axes. Not corrected for solar motion: Oviz applies Schönrich et al. (2010) itself |
| `name` | unique object name (used for search, links and member matching) |
| `age_myr` | age in Myr. Objects fade in at birth when time runs back |
| `n_stars` | optional: member count, for size scaling (`size_by_n_stars=True`) |

Objects without velocities (a fixed catalogue, a cloud outline) go in
`Layer(df, layer_name, assume_stationary=True)` instead of `Trace`.

## 2. Build the scene

```python
import numpy as np
from oviz import Scene3D, Trace, TraceCollection
from oviz.spiral_models import CASTRO_GINARD2021_ARMS, KHALIL2025_ARMS, khalil2025_potential

young = Trace(young_df, data_name="Clusters (< 15 Myr)", color="#ff0000", size_by_n_stars=True)
older = Trace(older_df, data_name="Clusters (< 60 Myr)", color="#2f80ff", size_by_n_stars=True)
scene = Scene3D(TraceCollection([older, young]), figure_theme="dark",
                potential=khalil2025_potential())   # omit for MWPotential2014

figure = scene.make_plot(
    time=np.arange(0, -61, -1.0),        # Myr; must contain 0; 1 Myr steps are plenty
    galactic_mode=True,                  # Galaxy-scale guides (GC, radius rings)
    enable_sky_panel=True,               # the registered sky view
    volumes=[dust_volume],               # see section 3
    cluster_members_file="members.csv",  # see section 4
    show_cluster_members_in_sky=True,
    spiral_arm_models=(KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS),  # hidden until toggled
)
figure.write_html("figure.html")
```

`make_plot` returns an `OvizFigure` (the WebGL2 viewer). Pass
`viewer="classic"` only if you need the older Three.js runtime's Slides or
Paper exports.

Other options worth knowing:
- `viewer_mode="detailed"`: open with every panel visible (default `"focus"`).
- `camera_anchor="lsr" | "sun" | "free"`: what the camera orbits and follows
  through time (default: the LSR).
- `show_sun=True`: adds the Sun as a trace.
- `fade_in_time`: Myr over which objects fade in at birth.

Potentials:
- The default is galpy's MWPotential2014 at R0 = 8.122 kpc, v0 = 236 km/s.
- `khalil2025_potential()` adds the Khalil et al. (2025) m = 2 and m = 3
  spirals. In the plane it reproduces their released SPIBACK model.
- Any galpy potential list works. The default reference frame (the LSR)
  always uses only its axisymmetric part, so the frame never wobbles.
- If you build galpy `SpiralArmsPotential` terms yourself, galpy negates N and
  alpha internally. Pass alpha = +pitch for trailing arms, and a density
  amplitude of K/(G·H·R0) for a Cox & Gomez depth K. `oviz.spiral_models`
  already does this.

## 3. Volumes (3D dust, emission, densities)

Each volume is a dict. Either point at a FITS cube:

```python
dust_volume = {"name": "Edenhofer+2024 Dust", "path": "mean_and_std_xyz.fits", "hdu": "MEAN",
               "max_resolution": 512, "opacity": 1.0, "samples": 200, "alpha_coef": 105,
               "colormap": "Greys", "vmin": 0, "vmax": 0.0385}
```

Or pass an array: `{"name": ..., "data": cube_zyx, "bounds": {"x": [-1250, 1250], "y": [...], "z": [...]}}`.
The cube is ordered (z, y, x) and the bounds are heliocentric pc.

`max_resolution` downsamples the cube. Keep cubes at or under 512³ for desktop
and expect phones to stride further. Hidden or present-day-only volumes still
ship in the file, so mind the size.

## 4. Sky view members

`cluster_members_file` is a CSV (or anything pandas reads) with one row per
star:
- a cluster-name column (`name`, `cluster_name`, `cluster` or `group_name`)
  matching the traces' `name`;
- sky coordinates: `l`/`b` or `ra`/`dec`;
- optionally `source_id`, `pmra`/`pmdec`, and `parallax` or a distance.

In Sky view each cluster marker crossfades into its members. Members never
appear in 3D, by design.

## 5. Saved views and presentation

Saved views ("States") capture the whole viewer: camera, time, layers,
styles, volumes, selection and Sky layers.
- **By hand:** open the figure, press **N** to save views and **Y** to open
  the story. **P** presents. **⌘S** saves a copy of the figure with its
  views, editable or present-only.
- **From a script:** run the browser API on an open figure.
  - `Oviz.viewer.states.add()`, `.goTo(i)` and `.present()`.
  - `Oviz.viewer.getState()` and `applyState(s)`.
  - `Oviz.viewer.states.exportHtml({download: false})` returns the HTML of a
    figure that includes the views.
- **To share:** `Oviz.viewer.viewLink()` gives a link to one view;
  `presentationLink()` gives one that opens presenting all saved views.

## 6. Large builds: science once, viewer many times

The October 1 build runs a slow science pipeline (catalog joins, orbit
integration, volume resampling) once. It writes a classic-format source figure
(`write_classic` in `tests/main_figure_oct1.py`). The viewer is then
regenerated from that source in seconds:

```python
from oviz import upgrade_html
upgrade_html("science_source.html", "figure.html", mode="focus", camera_anchor="lsr")
```

Use the same pattern whenever data processing dominates and you are iterating
on presentation. `python -m oviz.viewer.upgrade old.html new.html` also
converts figures made by older Oviz versions, saved views included.

You may need objects whose tracks come from elsewhere, such as the Oct 1
pulsars traced back along their proper motions. Add them to every frame of
the scene spec before compiling; `add_pulsar_traces` in
`tests/main_figure_oct1.py` shows the shape. Keep one stable `key` per trace,
give it a legend item, and set it to `"legendonly"` in every
`group_visibility` entry so it starts hidden.

## 7. Verify before you share

1. If you changed Oviz itself, run its tests: `pytest -q tests` and
   `node --test tests/viewer_js/*.test.mjs`.
2. Open the figure in a real browser with a real window size. A zero-size
   hidden iframe breaks Aladin. Check that:
   - the console is clean;
   - layers toggle from the key;
   - scrubbing time moves the objects;
   - Sky view loads its survey;
   - every saved view restores.
3. Check the HTML size. Under about 30 MB stays pleasant on phones.
4. If you changed the potential, compare a few tracks with the baseline. The
   present day must be identical; the past should differ smoothly.

## Pitfalls

- A time grid without 0, or one that isn't evenly spaced, fails or misaligns
  frames. galpy needs evenly spaced times.
- Velocities already corrected to the LSR or Galactic rest frame give wrong
  orbits. Oviz expects heliocentric UVW.
- Duplicate object names break member matching and search.
- Figures embed whatever provenance you attach. Record file names and
  versions, not absolute local paths.
- Large HTML files are slow on phones. Keep volumes compact (`max_resolution`)
  and leave rarely used layers out rather than hidden.
