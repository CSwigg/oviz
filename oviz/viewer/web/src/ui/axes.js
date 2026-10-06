// Axes (classic "Show axes"): X, Y and Z lines along the scene's range
// corner, ticked at round values in pc, with the figure's axis titles. For
// figures that are not Galactic maps, where the grid means nothing.

import { niceFloor } from "../core/math.js";
import { toggle } from "./controls.js";

const KEY = "__axes";

export class AxesPlugin {
  constructor() {
    this.name = "axes";
    this.built = false;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    const r = v.world?.ranges;
    this.ok = ["x", "y", "z"].every((k) => Array.isArray(r?.[k]) && r[k].length === 2 && r[k].every(Number.isFinite) && r[k][1] > r[k][0]);
    if (!this.ok) return;
    const sync = () => this.sync();
    v.on("global", sync);
    v.on("state-applied", sync);
    v.on("viewmode", sync);
    v.on("gpu-restored", () => { this.built = false; this.sync(); });
    this.sync();
  }

  settingsRows() {
    if (!this.ok) return [];
    const v = this.viewer;
    return [toggle({ label: "Axes", hint: "X, Y and Z in pc along the scene's edges", checked: !!v.state.global.axes, onChange: (x) => v.setGlobal({ axes: x }) })];
  }

  commands() {
    if (!this.ok) return [];
    const v = this.viewer;
    return [{ title: v.state.global.axes ? "Hide axes" : "Show axes", icon: "grid", keywords: "axis ticks x y z pc", run: () => v.setGlobal({ axes: !v.state.global.axes }) }];
  }

  sync() {
    const v = this.viewer;
    const on = !!v.state.global.axes && v.state.view.mode !== "sky";
    if (on && !this.built) this.build();
    const style = v.extraStyles.get(KEY);
    if (style) style.visible = on;
    for (const id of this.labels || []) {
      const el = v.labels.dynamic.get(id)?.el;
      if (el) el.hidden = !on;
    }
    v.renderer.invalidate();
  }

  build() {
    const v = this.viewer;
    const r = v.world.ranges;
    const titles = v.world.axisTitles || { x: "X (pc)", y: "Y (pc)", z: "Z (pc)" };
    const x0 = r.x[0], y0 = r.y[0], z0 = r.z[0];
    const tick = (v.world.maxSpan || 1000) * 0.02;
    const segs = [];
    const labels = [];
    const line = (a, b) => segs.push([a, b]);
    line([r.x[0], y0, z0], [r.x[1], y0, z0]);
    line([x0, r.y[0], z0], [x0, r.y[1], z0]);
    line([x0, y0, r.z[0]], [x0, y0, r.z[1]]);
    const axis = (k, at, off) => {
      const [lo, hi] = r[k];
      const step = niceFloor((hi - lo) / 5) || (hi - lo);
      for (let t = Math.ceil(lo / step - 1e-9) * step; t <= hi + 1e-9; t += step) {
        const p = at(t);
        line(p, off(p, 1));
        labels.push({ text: fmt(t), pos: off(p, 2.6) });
      }
    };
    axis("x", (t) => [t, y0, z0], (p, s) => [p[0], p[1] - tick * s, p[2]]);
    axis("y", (t) => [x0, t, z0], (p, s) => [p[0] - tick * s, p[1], p[2]]);
    axis("z", (t) => [x0, y0, t], (p, s) => [p[0] - tick * s, p[1], p[2]]);
    labels.push({ text: titles.x, pos: [r.x[1], y0 - tick * 5, z0], title: true });
    labels.push({ text: titles.y, pos: [x0 - tick * 5, r.y[1], z0], title: true });
    labels.push({ text: titles.z, pos: [x0 - tick * 5, y0, r.z[1]], title: true });
    const pos = new Float32Array(segs.length * 6);
    segs.forEach(([a, b], i) => { pos.set(a, i * 6); pos.set(b, i * 6 + 3); });
    const trace = { key: KEY, name: "Axes", color: "#8a93a6", showInLegend: false, lines: { count: segs.length, position: { frames: 1 }, color: "#8a93a6", width: 1.4, dash: "solid", opacity: 0.85 } };
    v.lines.removeTrace?.(KEY);
    v.lines.addTrace(trace, { position: pos, arc: new Float32Array(segs.length * 2) });
    v.extraStyles.set(KEY, { visible: true, presence: 1, opacityScale: 1, sizeScale: 1 });
    this.labels = labels.map((l, i) => {
      const id = `axis:${i}`;
      v.labels.setDynamic(id, l.text, () => l.pos, l.title ? "ov-label--axis ov-label--axis-title" : "ov-label--axis");
      return id;
    });
    this.built = true;
  }
}

function fmt(t) {
  const a = Math.abs(t);
  const s = a >= 1000 && a % 1000 === 0 ? `${t / 1000}k` : String(+t.toPrecision(6));
  return s.replace("-", "−");
}
