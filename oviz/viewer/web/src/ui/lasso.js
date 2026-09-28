// Lasso selection: draw around objects to select them; everything else dims.
// The selection bar offers framing, isolating and a CSV export.

import { h, icon, iconButton, downloadBlob } from "./dom.js";
import { clonePose } from "../engine/camera.js";

const LASSO_BIT = 8;

export class LassoPlugin {
  constructor() {
    this.name = "lasso";
    this.armed = false;
    this.selection = null; // Map traceKey → Uint32Array of indices
    this.isolate = false;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    this.svg.setAttribute("class", "ov-lasso");
    this.path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    this.svg.append(this.path);
    v.container.append(this.svg);
    this.btn = iconButton("wand", "Lasso select", () => this.arm(!this.armed), { cls: "ov-glass", shortcut: "X", pressed: false });
    ui.skyExtras.parentElement.querySelector(".ov-corner-btns:last-child")?.prepend(this.btn);
    this.bar = h("div", { class: "ov-selbar ov-glass", hidden: true, role: "region", "aria-label": "Lasso selection" });
    ui.ui.append(this.bar);
    v.stateExtensions.set("lasso", { capture: () => this.capture(), apply: (s) => this.restore(s) });
    v.on("gpu-restored", () => this.applyBits());
    this.bindPointer();
  }

  commands() {
    const list = [{ title: "Lasso select objects", icon: "wand", shortcut: "X", keywords: "select region polygon", run: () => this.arm(true) }];
    if (this.selection) {
      list.push(
        { title: "Export selection as CSV", icon: "download", run: () => this.exportCsv() },
        { title: "Clear lasso selection", icon: "close", run: () => this.clear() },
      );
    }
    return list;
  }

  onKey(e) {
    if ((e.key === "x" || e.key === "X") && !e.metaKey && !e.ctrlKey) { this.arm(!this.armed); return true; }
    if (e.key === "Escape" && this.armed) { this.arm(false); return true; }
    return false;
  }

  arm(on) {
    this.armed = on;
    this.btn.setAttribute("aria-pressed", String(on));
    this.ui.root.dataset.lasso = String(on);
    this.viewer.controls.enabled = !on;
    if (on) this.ui.toast("Draw around the objects to select", { icon: icon("wand"), ms: 1600 });
  }

  bindPointer() {
    const canvas = this.viewer.canvas;
    let pts = null;
    const local = (e) => {
      const r = canvas.getBoundingClientRect();
      return [e.clientX - r.left, e.clientY - r.top];
    };
    canvas.addEventListener("pointerdown", (e) => {
      if (!this.armed || e.button !== 0) return;
      e.stopPropagation();
      canvas.setPointerCapture(e.pointerId);
      pts = [local(e)];
    }, true);
    canvas.addEventListener("pointermove", (e) => {
      if (!pts) return;
      const p = local(e);
      const last = pts[pts.length - 1];
      if (Math.hypot(p[0] - last[0], p[1] - last[1]) > 3) pts.push(p);
      this.path.setAttribute("d", `M${pts.map((q) => q.join(",")).join("L")}Z`);
    });
    const end = () => {
      if (!pts) return;
      const poly = pts;
      pts = null;
      this.path.setAttribute("d", "");
      this.arm(false);
      if (poly.length > 3) this.select(poly);
    };
    canvas.addEventListener("pointerup", end);
    canvas.addEventListener("pointercancel", end);
  }

  select(poly) {
    const v = this.viewer;
    const cam = v.renderer.camera;
    const sel = new Map();
    let total = 0;
    const scr = [0, 0, 0];
    for (const trace of v.traces) {
      if (!trace.points || !v.state.traces[trace.key]?.visible) continue;
      const batch = v.points.byKey.get(trace.key);
      if (!batch) continue;
      const hits = [];
      for (let i = 0; i < trace.points.count; i++) {
        if (batch.state[i] & 2) continue; // hidden by the filter
        const p = v.objectPosition(trace.key, i);
        const s = p && cam.project(p, scr);
        if (s && inside(s[0], s[1], poly)) hits.push(i);
      }
      if (hits.length) { sel.set(trace.key, Uint32Array.from(hits)); total += hits.length; }
    }
    if (!total) {
      this.ui.toast("Nothing inside the lasso");
      return;
    }
    this.selection = sel;
    this.isolate = false;
    this.applyBits();
    this.renderBar();
  }

  clear() {
    this.selection = null;
    this.applyBits();
    this.renderBar();
  }

  applyBits() {
    const v = this.viewer;
    for (const batch of v.points.batches) {
      if (!this.selection) { batch.setStateBits(LASSO_BIT | 16, null); continue; }
      const bits = new Uint8Array(batch.count).fill(this.isolate ? 16 : LASSO_BIT);
      const idx = this.selection.get(batch.key);
      if (idx) for (const i of idx) bits[i] = 0;
      batch.setStateBits(LASSO_BIT | 16, bits);
    }
    v.invalidatePick();
    v.renderer.invalidate();
  }

  count() {
    let n = 0;
    for (const idx of this.selection?.values() || []) n += idx.length;
    return n;
  }

  renderBar() {
    const bar = this.bar;
    bar.replaceChildren();
    bar.hidden = !this.selection;
    if (!this.selection) return;
    bar.append(
      h("span", { class: "ov-selbar-count" }, icon("wand"), `${this.count().toLocaleString()} selected`),
      h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.frame() }, icon("target"), "Frame"),
      h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-pressed": String(this.isolate), onclick: () => { this.isolate = !this.isolate; this.applyBits(); this.renderBar(); } }, "Isolate"),
      h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.exportCsv() }, icon("download"), "CSV"),
      iconButton("close", "Clear selection", () => this.clear()),
    );
  }

  frame() {
    const v = this.viewer;
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const [key, idx] of this.selection) {
      for (const i of idx) {
        const p = v.objectPosition(key, i);
        if (!p) continue;
        for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); }
      }
    }
    if (!Number.isFinite(lo[0])) return;
    const pose = clonePose(v.pose);
    if (v.state.view.mode === "sky") return;
    pose.target = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
    const r = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2;
    pose.distance = Math.max(40, r / Math.tan((pose.fov * Math.PI) / 360) * 1.3);
    v.animateTo(pose, { duration: 1000 });
  }

  async exportCsv() {
    const v = this.viewer;
    const t = v.timeline.time;
    const rows = [["name", "layer", "age_now_myr", "age_at_t_myr", "n_stars", "l_deg", "b_deg", "dist_pc", "ra_deg", "dec_deg", `x_pc_at_t${t}`, `y_pc_at_t${t}`, `z_pc_at_t${t}`]];
    for (const [key, idx] of this.selection) {
      await v.objectMeta(key);
      for (const i of idx) {
        const d = v.describeObject(key, i);
        const p = d.position || [NaN, NaN, NaN];
        rows.push([d.name, d.traceName, d.ageNow, d.ageAt, d.nStars, d.l, d.b, d.dist, d.ra, d.dec, p[0], p[1], p[2]]);
      }
    }
    const csv = rows.map((r) => r.map(csvCell).join(",")).join("\n") + "\n";
    downloadBlob(new Blob([csv], { type: "text/csv" }), `oviz-selection-${rows.length - 1}.csv`);
    this.ui.toast(`Exported ${rows.length - 1} objects`, { icon: icon("download") });
  }

  capture() {
    if (!this.selection) return null;
    const out = {};
    for (const [k, idx] of this.selection) out[k] = Array.from(idx);
    return { selection: out, isolate: this.isolate };
  }

  restore(s) {
    if (!s || !s.selection) { if (this.selection) this.clear(); return; }
    this.selection = new Map(Object.entries(s.selection).map(([k, v]) => [k, Uint32Array.from(v)]));
    this.isolate = !!s.isolate;
    this.applyBits();
    this.renderBar();
  }
}

function inside(x, y, poly) {
  let c = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const [xi, yi] = poly[i], [xj, yj] = poly[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) c = !c;
  }
  return c;
}

function csvCell(v) {
  if (v == null || (typeof v === "number" && !Number.isFinite(v))) return "";
  if (typeof v === "number") return String(+v.toPrecision(10));
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
