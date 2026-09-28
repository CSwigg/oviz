// Interface styles ("skins"). Each style is a stylesheet layer keyed on
// <html data-oviz-skin> plus an optional layout (ui/layout.js) that
// arranges the shared components, and a few behaviour switches. Every
// style wears the same components, so each one keeps every feature. The
// figure chooses a default; `?ui=<id>` or the Display menu switches live.
//
//   layout     which arrangement from ui/layout.js to build (null: default)
//   selection  "panel" opens the inspector on click; "callout" shows a
//              small card beside the object, with details on demand
//   legend     show the compact legend of data layers
//   autoHide   fade the controls while the pointer rests (desktop)
//   layersOpen whether the layers panel starts open on desktop

export const SKINS = [
  {
    id: "observatory",
    name: "Observatory",
    desc: "The original: floating glass panels over the night sky.",
    swatch: ["#f4c46a", "#1c2230", "#0e1118"],
    layout: null, selection: "panel", legend: false, autoHide: false, layersOpen: true,
  },
  {
    id: "focus",
    name: "Focus",
    desc: "Just the figure and one quiet bar. Controls fade while you look; clicking shows a small label beside the object.",
    swatch: ["#ffd27a", "#f2f4f7", "#0c0e13"],
    layout: "focus", selection: "callout", legend: true, autoHide: true, layersOpen: false,
  },
  {
    id: "maps",
    name: "Maps",
    desc: "Works like a maps app: search and details at the top left, zoom and view buttons on the right, time along the bottom.",
    swatch: ["#8ab4f8", "#26282d", "#0c0e13"],
    layout: "maps", selection: "panel", legend: true, autoHide: false, layersOpen: false,
  },
  {
    id: "studio",
    name: "Studio",
    desc: "One tidy sidebar and a timeline along the bottom. The figure resizes into the space left, so nothing covers it.",
    swatch: ["#a5b0ff", "#16181e", "#0c0e13"],
    layout: "studio", selection: "panel", legend: false, autoHide: false, layersOpen: true,
  },
];

export function skinById(id) {
  return SKINS.find((s) => s.id === id) || SKINS[0];
}

/** The skin to start with: `?ui=` in the URL, else the figure's default. */
export function initialSkin() {
  let fromUrl = null;
  try { fromUrl = new URLSearchParams(location.search).get("ui"); } catch (_) { /* file:// quirks */ }
  const id = fromUrl || document.documentElement.dataset.ovizSkin || "observatory";
  return skinById(id).id;
}

export function applySkin(id) {
  const skin = skinById(id);
  document.documentElement.dataset.ovizSkin = skin.id;
  return skin.id;
}
