# Oviz agent guide

Use this file when an AI coding agent creates, edits, tests, exports, or
publishes an Oviz figure.

## Purpose

Oviz builds interactive HTML figures for Galactic 3D data. The normal scientific
workflow combines Gaia cluster or association data, galpy orbit integration,
the Oviz WebGL2 viewer, Aladin Lite sky backgrounds, optional member stars, and
optional 3D ISM volumes. Preserve the established visual style and the registered
3D/Sky behavior unless the user explicitly asks to change them.

Two runtimes exist:

- **Oviz viewer** (`oviz/viewer/`, the default). A custom WebGL2 engine and UI
  that read a compact binary bundle compiled from the scene spec.
- **Classic viewer** (`oviz/threejs_*.py`, `viewer="classic"`). The previous
  Three.js runtime. Keep it byte-for-byte stable: tests pin its exact source
  strings and function-body hashes. It alone provides Slides and Paper.
  AR (Apple Quick Look, USDZ) lives in the Oviz viewer: `web/src/ar/`.

## Start with the public API

Prefer `Trace` and `TraceCollection` for phase-space samples. Use `Layer` and
`LayerCollection` when a scene mixes orbiting and stationary spatial data.
`Scene3D` is the astronomy-facing alias of `Animate3D`.

```python
from oviz import Scene3D, Trace, TraceCollection, build_threejs_profile

trace = Trace(table, data_name="Young clusters", color="#e34a4a")
scene = Scene3D(
    TraceCollection([trace]),
    xyz_widths=(2000, 2000, 600),
    figure_theme="dark",
)
figure = scene.make_plot(
    time=time_myr,                 # must include 0
    galactic_mode=True,
    enable_sky_panel=True,
    renderer="threejs",
    threejs_initial_state=build_threejs_profile("full"),
    compress_scene_spec=True,
    show=False,
)
figure.write_html(output_path)
```

Required time-varying columns are `x`, `y`, `z` (pc), `U`, `V`, `W` (km/s),
`name`, and `age_myr`. `n_stars` is required only when
`size_by_n_stars=True`. For a static layer, provide XYZ and pass
`assume_stationary=True`.

## Scientific defaults

- Keep `t = 0` in every time grid and make time order explicit.
- Use a physically stated Galactic potential and record non-default `ro`, `vo`,
  and `zo` values.
- Do not silently change coordinate frames, sign conventions, units, sample
  cuts, velocities, ages, or source membership.
- Keep member stars out of 3D. When enabled, they replace or crossfade with the
  parent cluster marker only in Sky mode and inherit its visibility and styling.
- Treat volume layers as scientific data. Preserve their coordinate bounds,
  units, stretch, colormap, opacity transfer function, and time visibility.
- Keep Aladin backgrounds registered to Three.js traces. Do not fake spherical
  motion with CSS translation or scaling.

## Scene and runtime rules

- The Oviz viewer is the maintained runtime. Plotly paths exist for
  compatibility.
- Change the Python builders, the viewer sources in `oviz/viewer/web/`, or
  `oviz/threejs_runtime_*.py` for the classic runtime. Do not hand-edit a
  generated HTML artifact as the source of truth.

## Oviz viewer architecture

- `oviz/viewer/compile.py` turns the scene spec into a `Bundle`: a JSON
  manifest plus gzip (+ byte-shuffle) binary blobs.
  - Positions over time are stored static, rigid (base + per-frame offset), or
    per-frame.
  - Per-object attributes are stored once.
  - Volumes keep raw uint8 voxels plus a precomputed occupancy grid.
  - Nothing scientific is recomputed.
- `oviz/viewer/web/src` holds real ES modules (`core/`, `engine/`, `layers/`,
  `sky/`, `app/`, `ui/`). `oviz/viewer/build.py` bundles them into one script:
  each module gets its own scope, only named relative imports are allowed, and
  exported names must be unique. Styles live in `oviz/viewer/web/styles`
  (tokens in `00-tokens.css`).
- The viewer has one design (tokens in `00-tokens.css`: see-through glass,
  system type, calm motion) and two modes, keyed on `<html data-oviz-mode>`:
  `focus` (default: figure, key and one bar; selection shows a callout) and
  `detailed` (layers, inspector, toolbar and full transport in view). Mode
  CSS lives in `styles/80-mode-focus.css`, `81-mode-detailed.css` and
  `85-modes-mobile.css`; behaviour flags in `src/ui/modes.js`; arrangements
  in `src/ui/layout.js`, which moves (never rebuilds) the shared components
  and undoes every move on a mode or breakpoint change; the Python list is
  `figure.VIEWER_MODES`. On phones Focus keeps its single bar (the toolbar's
  actions live in More); phone-specific mode rules belong in
  `85-modes-mobile.css`. Never make a feature reachable in only one mode.
- `src/app/export.js` mirrors `web/template.html` for in-browser self-export.
  Keep the two in sync; `tests/test_viewer.py` checks this.
- Time is continuous: the vertex shader interpolates frame textures, so
  scrubbing must never re-upload geometry. Rendering is demand-driven: call
  `renderer.invalidate()` after state changes.
- Volume optics reproduce the classic model exactly: texture-space steps and
  the legacy composite. Do not "fix" the blend without the user's agreement;
  it changes tuned figures.
- Sky registration relies on a TAN projection with the same centre,
  horizontal FOV, and Galactic-north-up orientation. The camera sits at the
  origin in Sky view.
- The camera has an anchor (`src/engine/anchor.js`). It defaults to the LSR,
  which is the origin of Oviz's LSR-centred frame. The other anchors are the
  Sun, a followed object, and free.
  - While anchored, zoom keeps the anchor centred.
  - Time moves the camera by the anchor's displacement; it never snaps.
  - Pans and flights elsewhere free the camera; Home re-anchors it.
  - States and view links carry the anchor. Older ones infer it from the pose.
- Per-object GPU state is a bitfield: 1 = dimmed, 2 = hidden, 4 = replaced by
  member stars.
- A lost WebGL context must be recoverable. New GPU resources must be
  recreated in `Viewer._restoreGPU` or on the `gpu-restored` event.

Fast UI iteration: upgrade an existing figure and serve it over HTTP. There is
no need to rebuild the science.

```bash
python -m oviz.viewer.upgrade tests/main_figure_july25.html /tmp/july25.html
python -m http.server 8812 --bind 127.0.0.1 --directory /tmp
```
- Preserve States, Actions, Sky controls, presentation mode, mobile controls,
  lasso/selection behavior, and exact final-state restoration.
- A State captures the whole viewer. A State whose camera behavior is `keep`
  applies its other properties without replacing or stopping the live camera.
- Exported State-enabled HTML must remain editable unless a present-only export
  was explicitly requested.
- Avoid blocking first render on every remote HiPS layer. Remote tiles can be
  slow or unavailable.
- Desktop behavior is the baseline. Mobile controls should activate
  automatically on iPhone/mobile browsers without changing desktop layout.

## Canonical figure workflow

- The main figure runner belongs in `tests/main_figure.py`.
- The generated main figure belongs in `tests/main_figure.html`.
- Keep `scripts/main_figure.py` as a compatibility wrapper that delegates to
  `tests/main_figure.py`.
- Generate stored figures with compact payloads and keep
  `tests/main_figure.html` below 100 MB.
- Preserve unrelated dirty and untracked files. Stage only the source, focused
  tests, and canonical artifact required by the task.

## Verification

Run focused tests while iterating (`tests/test_viewer.py` and
`node --test tests/viewer_js/*.test.mjs` for the viewer), then the maintained
tracked suite:

```bash
git ls-files -z 'tests/test*.py' \
  | xargs -0 python -m pytest -q
```

Use the project environment when the default interpreter lacks astronomy
packages; on this workstation that is normally:

```bash
conda run -n p311 python -m pytest -q <test paths>
```

For viewer changes, regenerate the canonical HTML and inspect the relevant 3D,
Sky, time, State, presentation, and mobile workflows in a browser. Check the
artifact size and run `git diff --check` before committing.

## Publishing

When the user asks to upload an Oviz figure, prefer:

```bash
python scripts/upload_oviz_figure.py <path-to-html> --dry-run --no-push
python scripts/upload_oviz_figure.py <path-to-html>
```

The helper copies only the requested file to
`/Users/swiggumc/Desktop/astro_research/cam_website/oviz_figures`, checks its
size, commits that path, pushes the website repository, and prints the expected
GitHub Pages URL. Its default upload limit is 25 MiB, so use a compact or
mobile-safe export when needed. Verify the live URL or hash before reporting a
publication complete.
