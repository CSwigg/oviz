// Interface styles ("skins"). A skin is a stylesheet layer keyed on
// <html data-oviz-skin>; every skin wears the same components, so each one
// keeps every feature. The figure chooses a default; `?ui=<id>` or the
// Display menu switches at any time.

export const SKINS = [
  {
    id: "observatory",
    name: "Observatory",
    desc: "Quiet glass instruments floating over the night sky.",
    swatch: ["#f4c46a", "#1c2230", "#0e1118"],
  },
  {
    id: "spatial",
    name: "Spatial",
    desc: "Glass panels float at different depths and drift with your pointer; menus grow out of what you tap.",
    swatch: ["#ffffff", "#9fb8ff", "#141a28"],
  },
  {
    id: "island",
    name: "Island",
    desc: "A living capsule shows what is happening and flows open into search; iOS-style grouped lists.",
    swatch: ["#0a84ff", "#000000", "#101624"],
  },
  {
    id: "cards",
    name: "Cards",
    desc: "Panels are physical cards: past picks stack up behind, and you flick a card away to dismiss it.",
    swatch: ["#ff5a36", "#ffffff", "#101624"],
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
