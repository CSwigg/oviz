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

The worked reference is `tests/main_figure_oct1.py` in the repository
(https://github.com/CSwigg/oviz; `tests/` is not part of a pip install). It
builds the "October 1 figure", live at
https://cswigg.github.io/cam_website/oviz_figures/oviz_oct1.html:
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

Oviz needs Python ≥ 3.10; pip installs the rest (NumPy, pandas, Astropy,
galpy, SciPy, Matplotlib and a few small packages). Checking a figure needs a
browser (Chrome works headless). Sky imagery needs network access; the 3D
scene does not.

## 1. Prepare the data

Use a pandas DataFrame per group of objects, with these units and frames.
Getting them wrong silently produces wrong orbits.

| column | meaning |
| --- | --- |
| `x, y, z` | heliocentric Galactic Cartesian position, pc. x points to the Galactic centre, y to l = 90°, z to the north Galactic pole |
| `U, V, W` | heliocentric velocity, km/s, same axes. Not corrected for solar motion: Oviz applies Schönrich et al. (2010) itself |
| `name` | unique object name (used for search, links and member matching) |
| `age_myr` | age in Myr. Going back in time, each object fades out over the `fade_in_time` (default 5 Myr) before its birth |
| `n_stars` | optional: member count, for size scaling (`size_by_n_stars=True`) |

Objects without velocities (a fixed catalogue, a cloud outline) go in
`Layer(df, layer_name, assume_stationary=True)` (`from oviz import Layer,
LayerCollection`) instead of `Trace`.

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
    time=np.arange(0, -61, -1.0),        # Myr; evenly spaced, must contain 0; 1 Myr steps are plenty
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
- `figure_theme`: `"dark"`, `"gray"`, `"light"` or `"solarized_light"`.
- `viewer_mode="detailed"`: open with every panel visible (default `"focus"`).
- `camera_anchor="lsr" | "sun" | "free"`: what the camera orbits and follows
  through time (default: the LSR).
- `show_sun=False`: leave out the Sun, which is added by default.
- `fade_in_time`: Myr over which objects fade in before their birth.

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

- An array volume needs `bounds`. Without them each voxel counts as 1 pc
  around the Sun, so a 64³ cube spans only ±32 pc.
- A FITS cube takes its extent from the header (`CRPIX`, `CRVAL`, `CDELT` in
  pc on each axis) and ignores `bounds`. Without those keys it falls back to
  1 pc voxels. `clip_bounds` (same shape as `bounds`) crops it.
- `max_resolution` downsamples the cube. Keep cubes at or under 512³ for
  desktop and expect phones to stride further. Hidden or present-day-only
  volumes still ship in the file, so mind the size.

### KT maps and flow lines

A kinetic tomography (KT) map, gas density plus its line-of-sight velocity on
a heliocentric grid, has its own helper:

```python
from oviz.kt import read_kt_map

kt = read_kt_map("cubes_full_mean_std.h5")   # DIB density, v_los and the rotation residual
figure = scene.make_plot(
    time=time_myr,
    volumes=[kt.volume()],                   # opens coloured by v_los (RdBu_r)
    flows=[kt.flows(), kt.flows(field="residual", visible=False)],
)
```

- `read_kt_map` reads `(x, y, z)`-ordered HDF5 cubes on a Sun-centred grid
  with voxel sizes in `distances` (the DIB KT maps' layout). Other layouts go
  straight into `KTMap(density=..., velocity=..., residual=..., axes=(x, y, z))`
  with `(z, y, x)` arrays.
- `kt.volume()` is a volume dict whose opacity follows density and whose
  colour readers switch between the velocity, the residual and density
  (`color_by=` picks the opening one). It is block-averaged to
  `max_resolution` (velocities density-weighted).
- `kt.flows(field=...)` traces animated flow lines. A KT map measures only
  line-of-sight motion, so the lines run along sight lines: pulses leave the
  Sun where gas recedes and approach it where it approaches. Use
  `field="residual"` for motion relative to Galactic rotation. Both fade
  away from t = 0 and stay out of Sky view.
- Readers set each flow's colour, range, width, pulse speed, spacing and
  trail in the layers panel; States capture them.
- For a full 3D velocity field (for example a NIFTy reconstruction), use
  `oviz.viewer.flow.trace_streamlines` and `flow_layer` and pass the result
  in `flows=`.

## 4. Sky view members

`cluster_members_file` is a CSV (read with `pandas.read_csv`) with one row per
star:
- a cluster-name column (`name`, `cluster_name`, `cluster` or `group_name`)
  matching the traces' `name`, ignoring case and punctuation;
- sky coordinates in degrees: `l`/`b` or `ra`/`dec`;
- optionally `source_id`, `pmra`/`pmdec` (mas/yr), and `parallax` (mas) or a
  distance in pc (`distance_pc`, `r_med_geo`, `dist_pc` or `distance`; other
  names are ignored).

In Sky view each cluster marker crossfades into its members. Members never
appear in 3D, by design. A `ValueError` asking for "a readable
cluster_members_file" means the file could not be read, lacks the name or
coordinate columns, or no row matched a trace's object names.

## 5. Annotations

Labels, curves, arrows, bubbles and shells go in with `make_plot(annotations=[...])`.
Positions are pc in the figure's frame (heliocentric x, y, z at the present day):

```python
annotations = [
    {"kind": "shell", "center": [0, 0, 0], "radius": 150, "label": "Local Bubble", "present": True},
    {"kind": "arrow", "points": [[100, -300, 0], [50, -150, 20]], "color": "#ffd27a"},
    {"kind": "curve", "points": [[0, 0, 0], [200, 100, 50], [400, 0, 0]], "dash": "dash"},
    {"kind": "bubble", "center": [10, 0, 0], "radii": [30, 20, 10], "rot": [0, 0, 45]},
    {"kind": "text", "at": [300, 200, 0], "text": "Sco-Cen", "size": 16, "group": "Regions"},
    {"kind": "box", "center": [-200, 50, 0], "size": [300, 150, 80], "rot": [0, 0, 25], "select": "isolate"},
    {"kind": "curve", "points": [[0, -400, 0], [0, 400, 0]], "select": "highlight", "reach": 60},
]
```

- Kinds: `text`, `curve`, `arrow`, the ellipsoids `bubble`, `shell` and
  `wire` (`radius` or `radii` in pc), and `box` (`size`: edge lengths in pc).
  Shapes take `rot` in degrees about x, y, z.
- `select="highlight"` colours the objects a shape holds (a curve: within
  `reach` pc of it) like the annotation; `select="isolate"` shows only them.
  Both follow the objects through time.
- A shape's `label` adds a label at its centre, in the same key entry.
- Items with the same `group` share one key entry.
- `present: True` shows an item only around t = 0.
- Readers draw and edit annotations in the figure too: **K**, or the pen
  button. A label clicked onto a cluster follows that cluster's orbit, which
  Python input cannot do yet.

## 6. Saved views and presentation

Saved views ("States") capture the whole viewer: camera, time, layers,
styles, volumes, selection and Sky layers.
- **By hand:** open the figure, press **N** to save views and **Y** to open
  the story. **P** presents. **⌘S** saves a copy of the figure with its
  views, editable or present-only.
- **From a script:** run the browser API on an open figure.
  - `Oviz.viewer.states.add()`, `.goTo(target)` and `.present()`. `goTo`
    takes a 1-based index, a view's id or name, or `"original"`.
  - `Oviz.viewer.getState()` and `applyState(s)`.
  - `Oviz.viewer.states.exportHtml({download: false})` returns the HTML of a
    figure that includes the views.
- **To share:** `Oviz.viewer.viewLink()` gives a link to one view;
  `presentationLink()` gives one that opens presenting all saved views.

## 7. Large builds: science once, viewer many times

When data processing dominates (catalogue joins, orbit integration, volume
resampling), run it once and keep the result as a classic-format source
figure. The viewer is then regenerated from that source in seconds:

```python
from oviz import OvizFigure
from oviz.viewer import compile_scene_spec, read_legacy_scene_spec

scene.make_plot(..., viewer="classic").write_html("science_source.html")   # slow, once

spec = read_legacy_scene_spec("science_source.html")                        # fast, every time
# ...edit spec here: add traces, provenance, visibility...
OvizFigure(bundle=compile_scene_spec(spec), mode="focus", camera_anchor="lsr").write_html("figure.html")
```

Without edits, `oviz.upgrade_html("science_source.html", "figure.html")` does
the same. `python -m oviz.viewer.upgrade old.html new.html` also converts
figures made by older Oviz versions, saved views included. For a small edit to
a fresh figure, change `figure.scene_spec` before `write_html`.

To add objects whose tracks come from elsewhere, such as the October 1
pulsars traced back along their proper motions, add them to every frame of the
spec before compiling. `add_pulsar_traces` in `tests/main_figure_oct1.py`
shows the shape. Keep one stable `key` per trace, give it a legend item, and
set it to `"legendonly"` in every `group_visibility` entry so it starts
hidden. The figure embeds the spec's `provenance` verbatim, so record input
files by name and hash, not by local path.

## 8. Verify before you share

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

- A time grid without 0 raises an error; an unevenly spaced one fails in
  galpy.
- Velocities already corrected to the LSR or Galactic rest frame give wrong
  orbits. Oviz expects heliocentric UVW.
- Duplicate object names pass silently but confuse member matching, search
  and links.
- A member distance in kpc, or under an unrecognised column name, is read
  wrongly or ignored without a warning.
- Large HTML files are slow on phones. Keep volumes compact (`max_resolution`)
  and leave rarely used layers out rather than hidden.
