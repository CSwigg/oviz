# Oviz agent guide

Use this file when an AI coding agent creates, edits, tests, exports, or
publishes an Oviz figure.

## Purpose

Oviz builds interactive HTML figures for Galactic 3D data. The normal scientific
workflow combines Gaia cluster or association data, galpy orbit integration,
Three.js rendering, Aladin Lite sky backgrounds, optional member stars, and
optional 3D ISM volumes. Preserve the established visual style and the registered
3D/Sky behavior unless the user explicitly asks to change them.

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

- Three.js is the maintained renderer. Plotly paths exist for compatibility.
- Change Python builders or `oviz/threejs_runtime_*.py`; do not hand-edit a
  generated HTML artifact as the source of truth.
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

Run focused tests while iterating, then the maintained tracked suite:

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
