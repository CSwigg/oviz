// Distribution & filter: a live histogram of a per-object quantity for the
// visible layers, with a range brush that dims or hides everything outside.
// Replaces the classic Age KDE and Cluster Filter widgets.

import { h } from "./dom.js";
import { miniSeg } from "./controls.js";
import { parseColor } from "../core/color.js";
import { clamp } from "../core/math.js";

const PARAMS = [
  { key: "ageAt", label: "Age at t", unit: "Myr", needs: "ageNow" },
  { key: "ageNow", label: "Age today", unit: "Myr", needs: "ageNow" },
  { key: "nStars", label: "Members", unit: "stars", needs: "nStars", log: true },
  { key: "dist", label: "Distance from Sun", unit: "pc" },
  { key: "z", label: "Height z", unit: "pc" },
];

const BINS = 48;

export class FilterPlugin {
  constructor() {
    this.name = "filter";
    this.param = null;
    this.range = null; // [lo, hi] in value space, null = no filter
    this.mode = "dim";
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    // Objects only: samples of a surface or model (not pickable) have no distribution to filter.
    this.traces = v.traces.filter((t) => t.points && t.showInLegend && t.pickable !== false && t.points.count > 1);
    this.params = PARAMS.filter((p) => !p.needs || this.traces.some((t) => t.points[p.needs]));
    if (!this.traces.length || !this.params.length) return;
    this.param = this.params[0];
    const redraw = () => this.scheduleUpdate();
    v.on("time", () => { if (this.param?.key === "ageAt" || this.param?.key === "dist" || this.param?.key === "z") redraw(); });
    v.on("style", redraw);
    v.on("gpu-restored", redraw);
    v.on("theme", redraw); // axis text and range colours follow the theme
    v.stateExtensions.set("filter", { capture: () => this.captureState(), apply: (f) => this.applyState(f) });
    ui.filter = this;
  }

  commands() {
    if (!this.param) return [];
    return [
      { title: "Filter by age, members or distance", icon: "sliders", keywords: "histogram distribution range", run: () => { this.ui.setLayersOpen(true); this.el?.scrollIntoView({ block: "nearest", behavior: "smooth" }); } },
      ...(this.range ? [{ title: "Clear filter", icon: "close", run: () => this.setRange(null) }] : []),
    ];
  }

  /** Section element inserted into the layers panel. */
  section() {
    if (!this.param) return null;
    const sel = h("select", { "aria-label": "Quantity" }, this.params.map((p) => h("option", { value: p.key, selected: p === this.param }, p.label)));
    this.sel = sel;
    sel.addEventListener("change", () => {
      this.param = this.params.find((p) => p.key === sel.value);
      this.setRange(null);
    });
    this.canvas = h("canvas", { class: "ov-hist", width: 560, height: 160, "aria-label": "Distribution histogram. Drag to filter." });
    this.readout = h("div", { class: "ov-hist-readout" });
    this.clearBtn = h("button", { class: "ov-btn ov-btn--sm ov-btn--ghost", type: "button", onclick: () => this.setRange(null) }, "Clear");
    this.modeSeg = miniSeg({ value: this.mode, options: [{ value: "dim", label: "Dim" }, { value: "hide", label: "Hide" }], onChange: (m) => { this.mode = m; this.apply(); } });
    this.el = h("div", { class: "ov-section" },
      h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Distribution")),
      h("div", { class: "ov-filter" },
        h("div", { class: "ov-group", style: { margin: "0" } }, h("label", null, "By"), sel),
        this.canvas,
        h("div", { class: "ov-filter-foot" }, this.readout, h("span", { class: "ov-grow" }), this.modeSeg, this.clearBtn)));
    this.bindBrush();
    this.scheduleUpdate();
    return this.el;
  }

  // ------------------------------------------------------------ data

  values(trace) {
    const v = this.viewer;
    const d = v.data.get(trace.key);
    const N = trace.points.count;
    const out = new Float32Array(N);
    const t = v.timeline.time;
    const key = this.param.key;
    const meta = key === "dist" ? v.objectMetaSync(trace.key) : null;
    const sun = v.sunPosition();
    for (let i = 0; i < N; i++) {
      let x = NaN;
      if (key === "ageNow") x = d.ageNow ? d.ageNow[i] : NaN;
      else if (key === "ageAt") x = d.ageNow ? d.ageNow[i] + t : NaN;
      else if (key === "nStars") x = d.nStars && d.nStars[i] >= 0 ? d.nStars[i] : NaN;
      else if (key === "dist" || key === "z") {
        const p = v.objectPosition(trace.key, i);
        if (p) x = key === "z" ? p[2] : Math.hypot(p[0] - sun[0], p[1] - sun[1], p[2] - sun[2]);
        if (key === "dist" && Math.abs(t) < 1e-6 && meta?.dist?.[i] != null) x = meta.dist[i];
      }
      out[i] = x;
    }
    return out;
  }

  scheduleUpdate() {
    if (this._queued) return;
    this._queued = true;
    requestAnimationFrame(() => {
      this._queued = false;
      this.update();
    });
  }

  update() {
    if (!this.canvas || !this.canvas.isConnected) return;
    const v = this.viewer;
    const visible = this.traces.filter((t) => v.state.traces[t.key]?.visible);
    const series = visible.map((t) => ({ trace: t, values: this.values(t) }));
    let lo = Infinity, hi = -Infinity;
    const log = !!this.param.log;
    const tf = (x) => (log ? Math.log10(Math.max(x, 1e-9)) : x);
    for (const s of series) for (const x of s.values) if (x === x && (!log || x > 0) && !(this.param.key === "ageAt" && x < 0)) { lo = Math.min(lo, tf(x)); hi = Math.max(hi, tf(x)); }
    if (!(hi > lo)) { hi = lo + 1; }
    this.axis = { lo, hi, log };
    const counts = series.map((s) => {
      const c = new Float32Array(BINS);
      for (const x of s.values) {
        if (!(x === x) || (log && x <= 0) || (this.param.key === "ageAt" && x < 0)) continue;
        const b = clamp(Math.floor(((tf(x) - lo) / (hi - lo)) * BINS), 0, BINS - 1);
        c[b]++;
      }
      return c;
    });
    this.series = series;
    this.draw(counts, visible);
    this.apply();
  }

  draw(counts, traces) {
    const c = this.canvas;
    const ctx = c.getContext("2d");
    const W = c.width, H = c.height;
    ctx.clearRect(0, 0, W, H);
    const light = document.documentElement.dataset.ovizTheme === "light";
    const total = new Float32Array(BINS);
    for (const cs of counts) for (let b = 0; b < BINS; b++) total[b] += cs[b];
    const max = Math.max(1, ...total);
    const bw = W / BINS;
    const base = new Float32Array(BINS);
    counts.forEach((cs, k) => {
      const [r, g, b] = parseColor(traces[k].legendColor || traces[k].color);
      for (let i = 0; i < BINS; i++) {
        if (!cs[i]) continue;
        const inRange = this.inRangeBin(i);
        ctx.fillStyle = `rgba(${r * 255 | 0}, ${g * 255 | 0}, ${b * 255 | 0}, ${inRange ? 0.85 : 0.22})`;
        const y0 = H - 18 - (base[i] / max) * (H - 26);
        const hgt = (cs[i] / max) * (H - 26);
        ctx.fillRect(i * bw + 1, y0 - hgt, bw - 2, hgt);
        base[i] += cs[i];
      }
    });
    // Axis
    ctx.fillStyle = light ? "rgba(15,23,42,0.5)" : "rgba(255,255,255,0.45)";
    ctx.font = "600 18px system-ui, sans-serif";
    ctx.textBaseline = "bottom";
    const { lo, hi } = this.axis;
    ctx.textAlign = "left";
    ctx.fillText(this.fmtAxis(lo), 2, H);
    ctx.textAlign = "right";
    ctx.fillText(this.fmtAxis(hi), W - 2, H);
    if (this.range) {
      const [a, b] = this.range.map((x) => this.toX(x));
      ctx.fillStyle = light ? "rgba(184,121,27,0.12)" : "rgba(244,196,106,0.12)";
      ctx.fillRect(a, 0, b - a, H - 18);
      ctx.fillStyle = light ? "#b8791b" : "#f4c46a";
      ctx.fillRect(a - 1.5, 0, 3, H - 18);
      ctx.fillRect(b - 1.5, 0, 3, H - 18);
    }
    this.renderReadout();
  }

  fmtAxis(t) {
    const x = this.axis.log ? Math.pow(10, t) : t;
    const a = Math.abs(x);
    const s = a >= 1000 ? Math.round(x).toLocaleString() : a >= 10 ? x.toFixed(0) : x.toFixed(1);
    return `${s.replace("-", "−")} ${this.param.unit}`;
  }

  toX(t) {
    const { lo, hi } = this.axis;
    return ((t - lo) / (hi - lo)) * this.canvas.width;
  }

  fromX(px) {
    const { lo, hi } = this.axis;
    return lo + clamp(px / this.canvas.width, 0, 1) * (hi - lo);
  }

  inRangeBin(i) {
    if (!this.range) return true;
    const { lo, hi } = this.axis;
    const center = lo + ((i + 0.5) / BINS) * (hi - lo);
    return center >= this.range[0] && center <= this.range[1];
  }

  bindBrush() {
    const c = this.canvas;
    let start = null;
    const localX = (e) => {
      const r = c.getBoundingClientRect();
      return ((e.clientX - r.left) / r.width) * c.width;
    };
    c.addEventListener("pointerdown", (e) => {
      c.setPointerCapture(e.pointerId);
      start = this.fromX(localX(e));
      this.range = [start, start];
    });
    c.addEventListener("pointermove", (e) => {
      if (start == null) return;
      const cur = this.fromX(localX(e));
      this.range = [Math.min(start, cur), Math.max(start, cur)];
      this.scheduleUpdate();
    });
    const end = () => {
      if (start == null) return;
      start = null;
      if (this.range && Math.abs(this.range[1] - this.range[0]) < (this.axis.hi - this.axis.lo) / 200) this.range = null;
      this.scheduleUpdate();
    };
    c.addEventListener("pointerup", end);
    c.addEventListener("pointercancel", end);
    c.addEventListener("dblclick", () => this.setRange(null));
  }

  setRange(r) {
    this.range = r;
    this.update();
  }

  /** Keep the panel's quantity menu on the current parameter (other tools set it too). */
  syncSelect() {
    if (this.sel && this.param && this.sel.value !== this.param.key) this.sel.value = this.param.key;
  }

  /** Push dim/hide bits into every point batch. */
  apply() {
    const v = this.viewer;
    const log = this.axis?.log;
    const tf = (x) => (log ? Math.log10(Math.max(x, 1e-9)) : x);
    let kept = 0, total = 0;
    for (const s of this.series || []) {
      const batch = v.points.byKey.get(s.trace.key);
      if (!batch) continue;
      const bits = new Uint8Array(batch.count);
      if (this.range) {
        const flag = this.mode === "hide" ? 2 : 1;
        for (let i = 0; i < s.values.length; i++) {
          const x = s.values[i];
          const inside = x === x && tf(x) >= this.range[0] && tf(x) <= this.range[1];
          if (!inside) bits[i] = flag; else kept++;
          total++;
        }
      }
      batch.setStateBits(3, bits);
    }
    // Hidden layers keep no stale filter bits.
    for (const t of this.traces) {
      if (!(this.series || []).some((s) => s.trace === t)) v.points.byKey.get(t.key)?.setStateBits(3, null);
    }
    this.kept = kept;
    this.total = total;
    v.invalidatePick();
    v.renderer.invalidate();
    this.renderReadout();
  }

  renderReadout() {
    if (!this.readout) return;
    this.clearBtn.hidden = !this.range;
    if (!this.range) {
      const n = (this.series || []).reduce((s, x) => s + x.values.filter((y) => y === y).length, 0);
      this.readout.textContent = `${n.toLocaleString()} objects · drag to filter`;
      return;
    }
    const [a, b] = this.range;
    this.readout.textContent = `${this.fmtAxis(a)} – ${this.fmtAxis(b)} · ${this.kept?.toLocaleString() ?? 0} of ${this.total?.toLocaleString() ?? 0}`;
  }

  captureState() {
    return this.range ? { param: this.param.key, range: this.range.slice(), mode: this.mode } : null;
  }

  applyState(f) {
    if (!this.params) return;
    if (!f) { if (this.range) this.setRange(null); return; }
    this.param = this.params.find((p) => p.key === f.param) || this.param;
    this.mode = f.mode === "hide" ? "hide" : "dim";
    this.range = Array.isArray(f.range) ? f.range.slice(0, 2) : null;
    this.syncSelect();
    this.update();
  }
}
