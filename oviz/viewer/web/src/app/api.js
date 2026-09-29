// Public JS API: window.Oviz.viewer / window.Oviz.get(id).

export function installApi(root, viewer, ui) {
  const api = {
    root,
    viewer,
    ui,
    get time() { return viewer.timeline.time; },
    setTime: (t) => viewer.timeline.setTime(Number(t)),
    play: () => viewer.timeline.play(),
    pause: () => viewer.timeline.pause(),
    setViewMode: (m) => ui.setViewMode(m === "sky" ? "sky" : "3d"),
    // Whole-viewer States, the same snapshots saved views use.
    getState: () => ui.plugins.find((p) => p.name === "states")?.api().capture() ?? null,
    applyState: (s, opts) => ui.plugins.find((p) => p.name === "states")?.applyState(s, opts) ?? Promise.resolve(),
    screenshot: (opts) => viewer.renderer.capture(opts),
    // The current view as a USDZ for AR Quick Look: {blob, summary}.
    // {moment: true} for a still of the current time instead of the time-lapse.
    arModel: (opts) => ui.arModel(opts),
    viewInAr: (opts) => ui.viewInAr(opts),
    get states() { return ui.plugins.find((p) => p.name === "states")?.api(); },
    stats: () => ({
      frames: viewer.renderer.frameCount,
      avgFrameMs: viewer.renderer.fps(),
      lastFrameMs: viewer.renderer.lastFrameMs,
      volumeMarchMs: viewer.volumes.lastMarchMs,
      timing: viewer.timing,
    }),
  };
  const registry = (window.Oviz ||= { __viewers: new Map() });
  registry.__viewers ||= new Map();
  registry.__viewers.set(root.id, api);
  registry.viewer = api;
  registry.get = (id) => registry.__viewers.get(id || root.id);
  return api;
}
