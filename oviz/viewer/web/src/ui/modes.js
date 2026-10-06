// Modes: one design, two amounts of interface.
//
//   focus     the figure, its key and one quiet bar that fades while you
//             look; clicking an object shows a small label beside it
//   detailed  everything in view: the layers panel, a details panel for
//             the selection, the toolbar and every time control
//
// Each mode names a layout (ui/layout.js) and a few behaviour switches.
// Figures open in their default mode (focus unless the author chose
// otherwise); `?mode=detailed` in the URL overrides it, and U switches.

export const MODES = [
  {
    id: "focus",
    name: "Focus",
    desc: "The figure, its key and one quiet bar.",
    layout: "focus", selection: "callout", autoHide: true, layersOpen: false,
  },
  {
    id: "detailed",
    name: "Detailed",
    desc: "Layers, details and every control in view.",
    layout: "detailed", selection: "panel", autoHide: false, layersOpen: true,
  },
];

export function modeById(id) {
  return MODES.find((m) => m.id === id) || MODES[0];
}

/** The mode to start in: `?mode=` in the URL, else the figure's default. */
export function initialMode() {
  let fromUrl = null;
  try {
    const q = new URLSearchParams(location.search);
    fromUrl = q.get("mode") || q.get("ui");
  } catch (_) { /* file:// quirks */ }
  const id = fromUrl || document.documentElement.dataset.ovizMode || "focus";
  return modeById(id).id;
}

export function applyMode(id) {
  const mode = modeById(id);
  document.documentElement.dataset.ovizMode = mode.id;
  return mode.id;
}
