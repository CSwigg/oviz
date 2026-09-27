// Colour parsing (CSS strings → [r, g, b, a] in 0..1) with a small cache.

const cache = new Map();
let probe = null;

export function parseColor(value) {
  if (Array.isArray(value)) return value.length === 4 ? value : [...value, 1];
  const key = String(value ?? "").trim().toLowerCase();
  const hit = cache.get(key);
  if (hit) return hit;
  let out = null;
  let m;
  if ((m = /^#([0-9a-f]{3,8})$/.exec(key))) {
    let h = m[1];
    if (h.length === 3 || h.length === 4) h = [...h].map((c) => c + c).join("");
    const n = (i) => parseInt(h.slice(i, i + 2), 16) / 255;
    out = [n(0), n(2), n(4), h.length === 8 ? n(6) : 1];
  } else if ((m = /^rgba?\(([^)]+)\)$/.exec(key))) {
    const parts = m[1].split(/[\s,/]+/).filter(Boolean).map((p) => (p.endsWith("%") ? parseFloat(p) * 2.55 : parseFloat(p)));
    out = [parts[0] / 255, parts[1] / 255, parts[2] / 255, parts.length > 3 ? (m[1].includes("%") && parts[3] > 1 ? parts[3] / 255 : parts[3]) : 1];
  } else if (typeof document !== "undefined") {
    probe = probe || document.createElement("canvas").getContext("2d");
    probe.fillStyle = "#000";
    probe.fillStyle = key;
    const norm = probe.fillStyle;
    if (norm !== key) out = parseColor(norm);
  }
  if (!out || out.some((v) => !(v === v))) out = [1, 1, 1, 1];
  cache.set(key, out);
  return out;
}

export function toHex(c) {
  const [r, g, b] = parseColor(c);
  const h = (v) => Math.round(Math.min(Math.max(v, 0), 1) * 255).toString(16).padStart(2, "0");
  return `#${h(r)}${h(g)}${h(b)}`;
}

export function withAlpha(c, a) {
  const [r, g, b] = parseColor(c);
  return `rgba(${Math.round(r * 255)}, ${Math.round(g * 255)}, ${Math.round(b * 255)}, ${a})`;
}

/** Relative luminance (sRGB, WCAG) for choosing readable text colours. */
export function luminance(c) {
  const lin = (v) => (v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4));
  const [r, g, b] = parseColor(c);
  return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b);
}

/** Sample an RGBA LUT (Uint8Array width*4) at u ∈ [0, 1] → CSS rgb(). */
export function sampleLut(lut, u) {
  const w = lut.length / 4;
  const x = Math.min(Math.max(u, 0), 1) * (w - 1);
  const i = Math.floor(x), j = Math.min(i + 1, w - 1), t = x - i;
  const ch = (k) => Math.round(lut[i * 4 + k] * (1 - t) + lut[j * 4 + k] * t);
  return `rgb(${ch(0)}, ${ch(1)}, ${ch(2)})`;
}

/** CSS linear-gradient for a LUT (for legends and colorbars). */
export function lutGradient(lut, direction = "to right", stops = 12) {
  const parts = [];
  for (let i = 0; i <= stops; i++) {
    const u = i / stops;
    parts.push(`${sampleLut(lut, u)} ${(u * 100).toFixed(1)}%`);
  }
  return `linear-gradient(${direction}, ${parts.join(", ")})`;
}
