// Pure edits of the Sky layer stack (state.sky.layers, top first). The
// bottom visible layer is the background; visible layers above it blend
// over it with their opacity. Hidden layers stay in the list, so a figure
// remembers every survey it has shown.

import { normSurveyId, surveyInfo, SKY_SURVEYS } from "./catalog.js";

/** Index of the layer showing survey `id` (any spelling), or -1. */
export function findSurveyLayer(layers, id) {
  const n = normSurveyId(id);
  const info = surveyInfo(id);
  const alt = normSurveyId(info.id);
  return layers.findIndex((l) => {
    const m = normSurveyId(l.survey || l.key);
    return m === n || m === alt;
  });
}

/** A new layer entry for survey `id` (hidden until shown). */
export function newSurveyLayer(id, { label, opacity = 1 } = {}) {
  const info = surveyInfo(id);
  return { key: String(id), survey: String(id), label: label || info.label, opacity, visible: false };
}

export function isShown(l, backgroundVisible = true) {
  return backgroundVisible && l.visible !== false && (l.opacity ?? 1) > 0.001;
}

/** The layers in view, bottom (the background) first. */
export function layersInView(layers, backgroundVisible = true) {
  return layers.filter((l) => isShown(l, backgroundVisible)).reverse();
}

/**
 * Show survey `id` alone as the background: every other layer is hidden
 * (their opacities kept for later) and it is shown fully opaque. Returns
 * the new list; the layer keeps its place, so nothing re-attaches.
 */
export function showOnlySurvey(layers, id, meta = {}) {
  const out = layers.map((l) => ({ ...l }));
  let i = findSurveyLayer(out, id);
  if (i < 0) {
    out.push(newSurveyLayer(id, meta));
    i = out.length - 1;
  }
  out.forEach((l, j) => { l.visible = j === i; });
  out[i].opacity = 1;
  return out;
}

/**
 * Blend survey `id` over what is in view: it moves to the top of the stack
 * (only that layer re-attaches in Aladin) and shows at half strength, or at
 * the blend it had before. If nothing else is in view it simply becomes
 * the background.
 */
export function blendSurvey(layers, id, { opacity = 0.5, ...meta } = {}) {
  const out = layers.map((l) => ({ ...l }));
  const i = findSurveyLayer(out, id);
  const layer = i < 0 ? newSurveyLayer(id, meta) : out.splice(i, 1)[0];
  const others = out.filter((l) => isShown(l));
  layer.visible = true;
  layer.opacity = others.length ? (layer.opacity > 0.001 && layer.opacity < 0.999 ? layer.opacity : opacity) : 1;
  out.unshift(layer);
  return out;
}

/** Hide one layer (the next one down becomes the background). */
export function hideSurvey(layers, key) {
  return layers.map((l) => (l.key === key ? { ...l, visible: false } : { ...l }));
}

// ---------------------------------------------------------------- wavelength

/**
 * Stops of the wavelength slider: every all-sky survey with a known
 * wavelength among the figure's layers and the curated list, shortest
 * wavelength first (one per wavelength).
 */
export function wavelengthStops(layers) {
  const seen = new Set();
  const stops = [];
  const add = (id, label) => {
    const info = surveyInfo(id);
    const k = normSurveyId(info.id);
    if (seen.has(k) || !Number.isFinite(info.lambda) || info.spectrum === false) return;
    seen.add(k);
    stops.push({ id: info.known ? info.id : String(id), label: label || info.label, lambda: info.lambda, band: info.band });
  };
  for (const l of layers) add(l.survey || l.key, l.label);
  for (const s of SKY_SURVEYS) add(s.id);
  stops.sort((a, b) => a.lambda - b.lambda);
  // Two surveys at the same wavelength would make a dead stretch of slider.
  return stops.filter((s, i) => i === 0 || s.lambda > stops[i - 1].lambda * 1.0001);
}

/**
 * Crossfade the stack to position `pos` on the stops (0 … n−1): the two
 * neighbouring surveys are shown and all else hidden. Whichever of the two
 * sits higher in the stack takes its own weight as opacity over the other
 * at full opacity, so the result is (1 − f)·A + f·B in either direction
 * without reordering the stack (which would reload tiles).
 */
export function crossfadeToWavelength(layers, stops, pos) {
  if (stops.length < 1) return layers.map((l) => ({ ...l }));
  const n = stops.length;
  const p = Math.min(Math.max(Number(pos) || 0, 0), n - 1);
  const lo = Math.min(Math.floor(p), Math.max(0, n - 2));
  const f = n > 1 ? p - lo : 0;
  let out = layers.map((l) => ({ ...l }));
  // New surveys go on top, where Aladin adds them.
  const ensure = (stop) => {
    if (findSurveyLayer(out, stop.id) < 0) out = [newSurveyLayer(stop.id, { label: stop.label }), ...out];
  };
  const a = stops[lo], b = stops[Math.min(lo + 1, n - 1)];
  ensure(a);
  ensure(b);
  const ia = findSurveyLayer(out, a.id), ib = findSurveyLayer(out, b.id);
  out.forEach((l) => { l.visible = false; });
  if (ia === ib || f <= 0.0005) {
    out[ia].visible = true;
    out[ia].opacity = 1;
    return out;
  }
  if (f >= 0.9995) {
    out[ib].visible = true;
    out[ib].opacity = 1;
    return out;
  }
  const [top, bottom, wTop] = ia < ib ? [ia, ib, 1 - f] : [ib, ia, f];
  out[bottom].visible = true;
  out[bottom].opacity = 1;
  out[top].visible = true;
  out[top].opacity = Math.round(wTop * 1000) / 1000;
  return out;
}

/**
 * Where the stack sits on the wavelength slider, or null when what is in
 * view is not a single stop or a crossfade of two neighbouring stops.
 */
export function wavelengthPosition(layers, stops, backgroundVisible = true) {
  const view = layersInView(layers, backgroundVisible);
  const stopIndex = (l) => {
    const n = normSurveyId(surveyInfo(l.survey || l.key).id);
    return stops.findIndex((s) => normSurveyId(s.id) === n);
  };
  if (view.length === 1) {
    const i = stopIndex(view[0]);
    return i >= 0 ? i : null;
  }
  if (view.length !== 2) return null;
  const [bottom, top] = view;
  const ib = stopIndex(bottom), it = stopIndex(top);
  if (ib < 0 || it < 0 || Math.abs(ib - it) !== 1 || (bottom.opacity ?? 1) < 0.999) return null;
  const w = top.opacity ?? 1;
  return it > ib ? ib + w : ib - w;
}
