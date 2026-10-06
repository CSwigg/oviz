// Relative SFH (classic Oviz "Age KDE"): for each visible layer, a Gaussian
// kernel density of when its clusters formed (lookback time, each curve
// scaled to its own peak), restricted to the lasso selection when there is
// one, with the current time marked. Dragging across the plot filters the
// figure by age (the distribution filter's "Age today"); double-click clears.

import { h, clear } from "./dom.js";
import { slider } from "./controls.js";
import { clamp } from "../core/math.js";
import { rafThrottle } from "../core/emitter.js";

/** Gaussian KDE of `values` on `grid` with bandwidth `bw` (classic gaussianKde). */
export function gaussianKde(values, grid, bw) {
  const vals = [];
  for (const x of values) if (Number.isFinite(x)) vals.push(x);
  const out = new Float64Array(grid.length);
  if (!vals.length) return out;
  const b = Math.max(Number(bw) || 1, 1e-6);
  const norm = vals.length * b * Math.sqrt(2 * Math.PI);
  for (let g = 0; g < grid.length; g++) {
    let sum = 0;
    const x = grid[g];
    for (let i = 0; i < vals.length; i++) {
      const u = (x - vals[i]) / b;
      if (u > -8 && u < 8) sum += Math.exp(-0.5 * u * u);
    }
    out[g] = sum / norm;
  }
  return out;
}

/** Scale a curve to its own peak (the classic "relative" SFH). */
export function relativeCurve(density) {
  let max = 0;
  for (const y of density) max = Math.max(max, y);
  return max > 0 ? Array.from(density, (y) => y / max) : Array.from(density);
}

export class SfhPlugin {
  constructor() {
    this.name = "relative-sfh";
    this.bandwidth = null;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    const cfg = v.manifest.widgets?.ageKde || {};
    this.title = cfg.title || "Relative SFH";
    this.bandwidth = Number.isFinite(cfg.bandwidthMyr) ? cfg.bandwidthMyr : 2;
    this.xRange = Array.isArray(cfg.xRange) ? cfg.xRange.slice() : null;
    ui.widgets.register({
      id: "relative-sfh", title: this.title, icon: "chart",
      description: "When each layer's clusters formed: a density of birth times",
      size: { w: 480, h: 300 }, minSize: { w: 320, h: 220 },
      available: () => this.traces().length > 0,
      mount: (panel) => this.mount(panel),
    });
  }

  traces() {
    const v = this.viewer;
    return v.traces.filter((t) => t.points && t.showInLegend && t.points.ageNow?.blob && v.state.traces[t.key]?.inGroup !== false);
  }

  range() {
    if (this.xRange) return [Math.min(...this.xRange), Math.max(...this.xRange)];
    const tl = this.viewer.timeline;
    const lo = Math.min(tl.times[0], tl.times[tl.count - 1]);
    return [Math.min(lo, -1), 0];
  }

  /** One curve per visible layer: lookback times (−age) of the clusters it shows. */
  series() {
    const v = this.viewer;
    const lasso = this.ui.plugins.find((p) => p.name === "lasso");
    const sel = lasso?.selection && !lasso.filterOff ? lasso.selection : null;
    const [x0, x1] = this.range();
    const N = 300;
    const grid = Array.from({ length: N }, (_, i) => x0 + (i / (N - 1)) * (x1 - x0));
    const out = [];
    for (const t of this.traces()) {
      const ts = v.state.traces[t.key];
      if (!ts?.visible) continue;
      const d = v.data.get(t.key);
      if (!d?.ageNow) continue;
      const bits = v.points.byKey.get(t.key)?.state;
      const pick = sel ? sel.get(t.key) : null;
      if (sel && !pick) continue;
      const vals = [];
      // Only births the plot can show count (a few bandwidths of margin).
      const lo = x0 - 3 * this.bandwidth, hi = x1 + 3 * this.bandwidth;
      const add = (i) => {
        const age = d.ageNow[i];
        // Hidden or dimmed by the filter: not part of the SFH shown.
        if (!(age === age) || (bits && bits[i] & 3)) return;
        const tb = -Math.abs(age);
        if (tb >= lo && tb <= hi) vals.push(tb);
      };
      if (pick) for (const i of pick) add(i);
      else for (let i = 0; i < t.points.count; i++) add(i);
      if (!vals.length) continue;
      out.push({ trace: t, count: vals.length, x: grid, y: relativeCurve(gaussianKde(vals, grid, this.bandwidth)), color: ts.color || t.legendColor || t.color, opacity: clamp((ts.opacity ?? t.opacity) / (t.opacity || 1), 0.15, 1) });
    }
    out.sort((a, b) => a.trace.name.localeCompare(b.trace.name));
    return { grid, series: out, restricted: !!sel };
  }

  mount(panel) {
    const v = this.viewer;
    const body = panel.body;
    body.classList.add("ov-sfh");
    this.canvas = h("canvas", { class: "ov-sfh-canvas", "aria-label": "Relative star-formation history; drag to filter by age" });
    this.key = h("div", { class: "ov-sfh-key" });
    this.foot = h("div", { class: "ov-sfh-foot" });
    const bw = slider({ label: "Bandwidth", min: 0.25, max: 10, value: this.bandwidth, scale: "log", format: (x) => `${x.toFixed(x < 1 ? 2 : 1)} Myr`, onInput: (x) => { this.bandwidth = x; this.dirty = true; redraw(); } });
    bw.classList.add("ov-sfh-bw");
    body.append(h("div", { class: "ov-sfh-plot" }, this.canvas), h("div", { class: "ov-sfh-row" }, this.key, bw), this.foot);
    const redraw = rafThrottle(() => this.draw());
    this.redraw = redraw;
    const data = () => { this.dirty = true; redraw(); };
    const offs = [
      v.on("style", data), v.on("state-applied", data), v.on("theme", redraw), v.on("time", redraw),
      v.on("select", data), v.on("filter", data),
    ];
    // The lasso and the filter change the curves; poll their versions cheaply.
    this._poll = setInterval(() => {
      const sig = v.points.batches.map((b) => b.stateVersion || 0).join(",");
      if (sig !== this._bitsSig) { this._bitsSig = sig; data(); }
    }, 400);
    this.bindBrush();
    return {
      onOpen: () => { this.dirty = true; redraw(); },
      onResize: () => redraw(),
      capture: () => ({ bandwidth: this.bandwidth }),
      apply: (s) => { if (s && Number.isFinite(s.bandwidth)) { this.bandwidth = s.bandwidth; bw.set(s.bandwidth); data(); } },
      destroy: () => { offs.forEach((off) => off?.()); clearInterval(this._poll); },
    };
  }

  layout() {
    const cv = this.canvas;
    const W = cv.clientWidth || 1, H = cv.clientHeight || 1;
    const m = { left: 36, right: 12, top: 10, bottom: 22 };
    const [x0, x1] = this.range();
    const pw = Math.max(40, W - m.left - m.right), ph = Math.max(30, H - m.top - m.bottom);
    return {
      W, H, m, pw, ph, x0, x1,
      x: (t) => m.left + ((t - x0) / Math.max(x1 - x0, 1e-9)) * pw,
      t: (px) => x0 + ((px - m.left) / pw) * (x1 - x0),
      y: (val) => m.top + ph - clamp(val, 0, 1) * ph,
    };
  }

  draw() {
    const cv = this.canvas;
    if (!cv || !cv.isConnected) return;
    if (this.dirty || !this._data) { this._data = this.series(); this.dirty = false; }
    const { series, restricted } = this._data;
    const L = this.layout();
    const dpr = Math.max(window.devicePixelRatio || 1, 1);
    if (cv.width !== Math.round(L.W * dpr) || cv.height !== Math.round(L.H * dpr)) { cv.width = Math.round(L.W * dpr); cv.height = Math.round(L.H * dpr); }
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    const css = getComputedStyle(this.ui.root);
    const axis = css.getPropertyValue("--ov-text-3").trim() || "rgba(236,240,246,0.5)";
    const accent = css.getPropertyValue("--ov-accent").trim() || "#ffd27a";
    const { m, pw, ph } = L;
    // The age filter as a band (lookback time = −age).
    const f = this.ui.filter;
    if (f?.range && f.param?.key === "ageNow") {
      const a = L.x(clamp(-f.range[1], L.x0, L.x1)), b = L.x(clamp(-f.range[0], L.x0, L.x1));
      ctx.fillStyle = accent; ctx.globalAlpha = 0.1;
      ctx.fillRect(Math.min(a, b), m.top, Math.abs(b - a), ph);
      ctx.globalAlpha = 1;
    }
    ctx.strokeStyle = axis; ctx.fillStyle = axis; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(m.left, m.top); ctx.lineTo(m.left, m.top + ph); ctx.lineTo(m.left + pw, m.top + ph); ctx.stroke();
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (const fr of [0, 0.5, 1]) {
      const py = L.y(fr);
      ctx.globalAlpha = 0.14; ctx.beginPath(); ctx.moveTo(m.left, py); ctx.lineTo(m.left + pw, py); ctx.stroke(); ctx.globalAlpha = 1;
      ctx.fillText(fr.toFixed(1), m.left - 6, py);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (const t of [L.x0, (L.x0 + L.x1) / 2, L.x1]) ctx.fillText(`${Math.abs(t) < 1e-9 ? "0" : t.toFixed(0).replace("-", "−")} Myr`, L.x(t), m.top + ph + 6);
    for (const s of series) {
      ctx.beginPath();
      ctx.moveTo(L.x(s.x[0]), m.top + ph);
      for (let i = 0; i < s.x.length; i++) ctx.lineTo(L.x(s.x[i]), L.y(s.y[i]));
      ctx.lineTo(L.x(s.x[s.x.length - 1]), m.top + ph);
      ctx.closePath();
      ctx.fillStyle = s.color; ctx.globalAlpha = Math.min(0.28, 0.14 + 0.18 * s.opacity); ctx.fill();
      ctx.beginPath();
      for (let i = 0; i < s.x.length; i++) { const px = L.x(s.x[i]), py = L.y(s.y[i]); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
      ctx.strokeStyle = s.color; ctx.globalAlpha = s.opacity; ctx.lineWidth = 2; ctx.stroke();
    }
    ctx.globalAlpha = 1;
    // Now: a dashed line at the current time.
    const t = this.viewer.timeline.time;
    const tx = L.x(clamp(t, L.x0, L.x1));
    ctx.save(); ctx.setLineDash([5, 5]); ctx.strokeStyle = axis; ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.moveTo(tx, m.top); ctx.lineTo(tx, m.top + ph); ctx.stroke(); ctx.restore();
    if (this._brush) {
      const [a, b] = this._brush.map((x) => L.x(x));
      ctx.fillStyle = accent; ctx.globalAlpha = 0.2; ctx.fillRect(Math.min(a, b), m.top, Math.abs(b - a), ph); ctx.globalAlpha = 1;
    }
    clear(this.key);
    for (const s of series) this.key.append(h("span", { class: "ov-sfh-item" }, h("i", { style: { background: s.color } }), `${s.trace.name} · ${s.count.toLocaleString()}`));
    if (!series.length) this.key.append(h("span", { class: "ov-faint" }, "No visible layers with ages"));
    this.foot.textContent = `${restricted ? "Lasso selection · " : ""}each curve scaled to its peak · drag to filter by age, double-click to clear`;
  }

  bindBrush() {
    const cv = this.canvas;
    let start = null;
    const at = (e) => {
      const r = cv.getBoundingClientRect();
      const L = this.layout();
      return clamp(L.t(e.clientX - r.left), L.x0, L.x1);
    };
    cv.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      start = at(e);
      this._brush = [start, start];
      try { cv.setPointerCapture(e.pointerId); } catch (_) { /* synthetic pointer */ }
      this.redraw();
    });
    cv.addEventListener("pointermove", (e) => {
      if (start == null) return;
      this._brush = [start, at(e)];
      this.redraw();
    });
    const end = () => {
      if (start == null) return;
      const [a, b] = this._brush || [start, start];
      start = null;
      this._brush = null;
      const L = this.layout();
      if (Math.abs(b - a) < (L.x1 - L.x0) / 150) { this.redraw(); return; }
      // Lookback time t is −age: the range [t0, t1] is ages [−t1, −t0].
      this.setAgeFilter([Math.max(0, -Math.max(a, b)), Math.max(0, -Math.min(a, b))]);
    };
    cv.addEventListener("pointerup", end);
    cv.addEventListener("pointercancel", end);
    cv.addEventListener("dblclick", () => this.setAgeFilter(null));
  }

  setAgeFilter(range) {
    const f = this.ui.filter;
    if (!f) return;
    if (!range) {
      if (f.param?.key === "ageNow") f.setRange(null);
    } else {
      f.applyState({ param: "ageNow", range, mode: f.mode || "dim" });
    }
    this.dirty = true;
    this.redraw();
  }
}
