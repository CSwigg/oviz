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
    id: "orbit",
    name: "Orbit",
    desc: "Soft and playful: rounded pills, a breathing play orb and springy panels.",
    swatch: ["#a78bff", "#ff8fd0", "#16123a"],
  },
  {
    id: "instrument",
    name: "Instrument",
    desc: "Mission-control HUD: edge-docked panels, live target readouts, phosphor type.",
    swatch: ["#6effc0", "#ffb547", "#030a06"],
  },
  {
    id: "atlas",
    name: "Atlas",
    desc: "An editorial page framing the sky, with serif captions and a ruler timeline.",
    swatch: ["#d9472b", "#1f1b14", "#f5f0e6"],
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
