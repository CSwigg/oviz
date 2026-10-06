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
    // What the camera orbits, zooms toward and moves with through time:
    // {kind: "lsr" | "sun" | "free"} or {kind: "object", trace, index}.
    get cameraAnchor() { return { ...viewer.anchor }; },
    setCameraAnchor: (a) => {
      const anchor = a && typeof a === "object" && a.kind === undefined ? { kind: "object", ...a } : a;
      if (anchor === "free" || anchor?.kind === "free") {
        viewer.releaseAnchor();
        return Promise.resolve({ done: true });
      }
      return viewer.setAnchor(anchor);
    },
    // Whole-viewer States, the same snapshots saved views use.
    getState: () => ui.plugins.find((p) => p.name === "states")?.api().capture() ?? null,
    applyState: (s, opts) => ui.plugins.find((p) => p.name === "states")?.applyState(s, opts) ?? Promise.resolve(),
    screenshot: (opts) => viewer.renderer.capture(opts),
    // A link to the current view (everything in it), as "Copy link to this view" makes.
    viewLink: (opts) => ui.viewLinkUrl(opts),
    // A link that opens presenting the saved views (they travel in the link).
    presentationLink: () => ui.viewLinkUrl({ present: true }),
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
  installMessageBridge(api);
  const registry = (window.Oviz ||= { __viewers: new Map() });
  registry.__viewers ||= new Map();
  registry.__viewers.set(root.id, api);
  registry.viewer = api;
  registry.get = (id) => registry.__viewers.get(id || root.id);
  return api;
}

/**
 * Commands from the page embedding this figure in an iframe (the classic
 * Oviz postMessage bridge). The parent posts `{ type: "oviz", command,
 * args, id }` to the frame and gets `{ type: "oviz:done", command, id }`
 * back when it has run. Only the parent window is heard, and only these
 * named commands exist (nothing is evaluated).
 */
function installMessageBridge(api) {
  if (window.parent === window) return;
  const commands = {
    goTo: (target) => api.states?.goTo(target),
    next: () => api.states?.next(),
    previous: () => api.states?.previous(),
    original: () => api.states?.original?.(),
    present: (on = true) => api.states?.present(on !== false),
    setTime: (t) => api.setTime(t),
    play: () => api.play(),
    pause: () => api.pause(),
    setViewMode: (m) => api.setViewMode(m),
  };
  window.addEventListener("message", (e) => {
    if (e.source !== window.parent) return;
    const d = e.data;
    if (!d || typeof d !== "object" || d.type !== "oviz" || !Object.prototype.hasOwnProperty.call(commands, d.command)) return;
    const args = Array.isArray(d.args) ? d.args : d.args === undefined ? [] : [d.args];
    Promise.resolve()
      .then(() => commands[d.command](...args))
      .then(() => window.parent.postMessage({ type: "oviz:done", command: d.command, id: d.id ?? null }, "*"))
      .catch(() => {});
  });
}
