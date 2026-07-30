# Claude guide for Oviz

Read [AGENTS.md](AGENTS.md) before changing this repository. It contains the
authoritative build, verification, artifact, and publishing rules.

## Build a figure

1. Put Gaia-derived phase-space rows in a `Trace` (`x`, `y`, `z` in pc;
   `U`, `V`, `W` in km/s; plus `name` and `age_myr`).
2. Combine traces in `TraceCollection` and pass them to `Scene3D` or
   `Animate3D`.
3. Call `make_plot()` with a time array containing zero, the Three.js renderer,
   and an explicit initial-state profile.
4. Write the returned `ThreeJSFigure` to one compact HTML file.

```python
from oviz import Scene3D, Trace, TraceCollection, build_threejs_profile

trace = Trace(cluster_table, "Clusters", color="#e34a4a")
scene = Scene3D(TraceCollection([trace]), figure_theme="dark")
figure = scene.make_plot(
    time=time_myr,
    renderer="threejs",
    galactic_mode=True,
    enable_sky_panel=True,
    threejs_initial_state=build_threejs_profile("full"),
    compress_scene_spec=True,
)
figure.write_html(output_path)
```

Use `Layer(..., assume_stationary=True)` for XYZ data that should not be orbit
integrated. Member-star tables need a cluster-name field plus `l`/`b` or
`ra`/`dec`. Volume maps need explicit bounds and rendering controls.

## Do not break these invariants

- Member stars render only in Sky mode; bulk cluster markers render in 3D.
- Aladin imagery and Three.js traces remain registered during pan, zoom, and
  3D/Sky transitions.
- States restore camera, time, visibility, styles, volumes, selections, Sky
  layers, panels, and presentation settings exactly.
- Present-only exports begin at the first saved State and navigate States in
  order.
- Desktop appearance stays unchanged when mobile behavior is added.
- Generated HTML is an artifact, not the primary implementation.

Prefer narrow changes, focused regression tests, a regenerated canonical
artifact when runtime output changes, and the full maintained test suite before
commit. Preserve unrelated working-tree files.
