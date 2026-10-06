// Lasso selection: draw around objects to select them (Shift+drag at any
// time, or L to keep lassoing until L or Esc, as in classic Oviz). The rest
// hide, like the classic filter, or dim. Volumes are selected with the
// objects: the dust outside the lasso hides too (Display settings can turn
// that off), and a lasso around dust alone selects just that dust. The
// selection bar offers framing, dimming instead of hiding and a CSV export.

import { h, icon, iconButton, downloadBlob } from "./dom.js";
import { toggle } from "./controls.js";
import { clonePose } from "../engine/camera.js";
import { birthFadeAt } from "../app/timeline.js";
import { selectionKey } from "../app/state.js";

const LASSO_BIT = 8;
// A lasso fade cut short (its State change interrupted) finishes over this.
const SETTLE_MS = 400;

export class LassoPlugin {
  constructor() {
    this.name = "lasso";
    this.armed = false;
    this.selection = null; // Map traceKey → Uint32Array of indices
    this.isolate = true; // hide the rest (classic); false dims it
    this.mask = null; // the lasso in NDC + the camera it was drawn with (clips volumes)
    this.restoreToken = 0; // a newer State wins over a classic selection still resolving
    this.blend = null; // while the drawing fades between selections (a newer fade replaces it)
    this.pendingSettle = null;
    this.stopSettle = null;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    this.svg.setAttribute("class", "ov-lasso");
    this.path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    this.svg.append(this.path);
    v.container.append(this.svg);
    this.btn = iconButton("wand", "Lasso select", () => this.arm(!this.armed), { cls: "ov-glass", shortcut: "L", pressed: false });
    ui.skyExtras.parentElement.querySelector(".ov-corner-btns:last-child")?.prepend(this.btn);
    this.bar = h("div", { class: "ov-selbar ov-glass", hidden: true, role: "region", "aria-label": "Lasso selection" });
    ui.ui.append(this.bar);
    v.stateExtensions.set("lasso", {
      capture: () => this.capture(),
      apply: (s) => this.restore(s),
      transition: (from, to) => this.transition(from, to),
    });
    v.on("gpu-restored", () => this.applyBits());
    // The bar names the dust in the selection, which follows the volumes shown.
    v.on("volume", () => { if (this.selection) this.renderBar(); });
    this.bindPointer();
  }

  commands() {
    const list = [{ title: "Lasso select objects", icon: "wand", shortcut: "L", keywords: "select region polygon", run: () => this.arm(true) }];
    if (this.selection) {
      list.push(
        { title: "Export selection as CSV", icon: "download", run: () => this.exportCsv() },
        { title: "Clear lasso selection", icon: "close", run: () => this.clear() },
      );
    }
    return list;
  }

  settingsRows() {
    const v = this.viewer;
    return [toggle({
      label: "Lasso selects volumes too", hint: "Dust outside a lasso hides with the other objects",
      checked: v.state.global.lassoVolumes !== false,
      onChange: (x) => { v.setGlobal({ lassoVolumes: x }); this.applyMask(); this.renderBar(); },
    })];
  }

  onKey(e) {
    if (e.metaKey || e.ctrlKey) return false;
    const k = e.key.toLowerCase();
    // L arms the lasso (classic key; X is an alias), C toggles its filter.
    if ((k === "l" && !e.shiftKey) || k === "x") { this.arm(!this.armed); return true; }
    if (k === "c" && !e.shiftKey && this.selection) { this.toggleFilter(); return true; }
    if (e.key === "Escape" && this.armed) { this.arm(false); return true; }
    return false;
  }

  /** Show or hide the dimming while keeping the selection (classic "C"). */
  toggleFilter() {
    this.ui.pushSelectionUndo?.();
    this.filterOff = !this.filterOff;
    this.applyBits();
    this.renderBar();
    this.ui.toast(`Selection filter: ${this.filterOff ? "off" : "on"}`, { ms: 1100 });
  }

  arm(on) {
    this.armed = on;
    this.btn.setAttribute("aria-pressed", String(on));
    this.ui.root.dataset.lasso = String(on);
    this.viewer.controls.enabled = !on;
    // Touch screens have no L or Esc: there a lasso is one stroke.
    const touch = typeof matchMedia === "function" && matchMedia("(pointer: coarse)").matches;
    if (on) this.ui.toast(touch ? "Draw around objects with a finger to select them" : "Draw around objects to select · L or Esc to stop lassoing", { icon: icon("wand"), ms: 2200 });
  }

  bindPointer() {
    const canvas = this.viewer.canvas;
    let pts = null;
    const local = (e) => {
      const r = canvas.getBoundingClientRect();
      return [e.clientX - r.left, e.clientY - r.top];
    };
    let shiftStroke = false;
    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      // Shift+drag lassos without arming (classic). The press still reaches
      // the figure, so a Shift+click without dragging measures as before;
      // only the orbit is held off for the stroke.
      const shift = !this.armed && e.shiftKey && e.pointerType !== "touch" && this.viewer.state.view.mode !== "sky";
      if (!this.armed && !shift) return;
      shiftStroke = shift;
      if (shift) this.viewer.controls.enabled = false;
      else e.stopPropagation();
      try { canvas.setPointerCapture(e.pointerId); } catch (_) { /* synthetic or finished pointer */ }
      pts = [local(e)];
    }, true);
    canvas.addEventListener("pointermove", (e) => {
      if (!pts) return;
      const p = local(e);
      const last = pts[pts.length - 1];
      if (Math.hypot(p[0] - last[0], p[1] - last[1]) > 3) pts.push(p);
      this.path.setAttribute("d", `M${pts.map((q) => q.join(",")).join("L")}Z`);
    });
    const end = (e) => {
      if (!pts) return;
      const poly = pts;
      pts = null;
      this.path.setAttribute("d", "");
      // An armed lasso stays armed for the next stroke (L or Esc ends it);
      // a finger's stroke ends it, since touch screens have neither key.
      if (shiftStroke) this.viewer.controls.enabled = !this.armed;
      shiftStroke = false;
      if (this.armed && e?.pointerType === "touch") this.arm(false);
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
    const t = v.timeline.time;
    const g = v.state.global;
    for (const trace of v.traces) {
      // Samples of a surface or model (not pickable) are not objects to select.
      if (!trace.points || trace.pickable === false || !v.state.traces[trace.key]?.visible) continue;
      const batch = v.points.byKey.get(trace.key);
      if (!batch) continue;
      const ages = v.data.get(trace.key)?.ageNow;
      const hits = [];
      for (let i = 0; i < trace.points.count; i++) {
        if (batch.state[i] & 18) continue; // hidden by the filter or an isolated lasso
        // Not yet born (or long faded) at this time: invisible, so not selectable.
        if (ages && birthFadeAt(t, ages[i], g.fadeTime, g.fadeInOut) <= 0) continue;
        const p = v.objectPosition(trace.key, i);
        const s = p && cam.project(p, scr);
        if (s && inside(s[0], s[1], poly)) hits.push(i);
      }
      if (hits.length) { sel.set(trace.key, Uint32Array.from(hits)); total += hits.length; }
    }
    // A lasso around dust alone still selects that dust.
    if (!total && !this.lassoedVolumes().length) {
      this.ui.toast("Nothing inside the lasso");
      return;
    }
    this.ui.pushSelectionUndo?.();
    this.restoreToken++;
    this.selection = sel;
    this.filterOff = false;
    // The volumes are clipped to the lasso as it was drawn: its outline in
    // normalised device coordinates and the camera of that moment.
    const W = cam.width || 1, H = cam.height || 1;
    this.mask = { polygon: poly.map(([x, y]) => [(x / W) * 2 - 1, 1 - (y / H) * 2]), vp: Array.from(cam.viewProj) };
    this.applyBits();
    this.renderBar();
  }

  clear() {
    this.restoreToken++;
    this.selection = null;
    this.mask = null;
    this.applyBits();
    this.renderBar();
  }

  applyBits() {
    const v = this.viewer;
    // Bits for the CPU's readers (Sky members, AR, the next lasso); the
    // shaders draw from the levels, which a State change can blend.
    for (const batch of v.points.batches) {
      if (!this.selection || this.filterOff) batch.setStateBits(LASSO_BIT | 16, null);
      else {
        const bits = new Uint8Array(batch.count).fill(this.isolate ? 16 : LASSO_BIT);
        const idx = this.selection.get(batch.key);
        if (idx) for (const i of idx) bits[i] = 0;
        batch.setStateBits(LASSO_BIT | 16, bits);
      }
    }
    // A fade cut short (its State change interrupted) settles onto the
    // selection from where it was; otherwise the drawing lands at once.
    if (this.fading()) this.settleSoon();
    else this.land();
    v.invalidatePick();
    v.renderer.invalidate();
  }

  /** Draw the selection as it is, ending any fade. */
  land() {
    const v = this.viewer;
    const levels = this.levelsFor(this.current());
    this.stopSettle?.();
    this.blend = null;
    this.pendingSettle = null;
    for (const batch of v.points.batches) batch.setLasso(levels?.get(batch.key) || null);
    v.points.setLassoMix(1, 1);
    this.applyMask();
  }

  /** Whether the drawing is part-way through a fade between selections. */
  fading() {
    const mix = this.viewer.points.lassoMix;
    return !!this.blend && (mix[0] < 1 || mix[1] < 1);
  }

  /**
   * Start fading the drawing from what it shows now (a selection, or a
   * fade part-way) to `levels` and `clip` (see `levelsFor`, `clipFor`).
   * Returns the fade's stepper, `(out, inn)` progress, which goes inert
   * once a newer fade takes over.
   */
  fadeTo(levels, clip) {
    const v = this.viewer;
    const styles = v.points.params?.styles;
    const mix = v.points.lassoMix;
    for (const batch of v.points.batches) {
      const dim = styles?.get(batch.key)?.dimOpacity ?? 0.16;
      batch.setLasso(levels?.get(batch.key) || null, batch.lassoDrawn(mix, dim));
    }
    v.points.setLassoMix(0, 0);
    v.volumes.blendMask(clip);
    this.stopSettle?.();
    const blend = (this.blend = {});
    return (out, inn) => {
      if (this.blend !== blend) return;
      v.points.setLassoMix(out, inn);
      v.volumes.setMaskMix(out, inn);
    };
  }

  /**
   * Settle a fade cut short on the next frame, unless a State change starts
   * its own fade from here first (next/previous pressed mid-flight).
   */
  settleSoon() {
    const v = this.viewer;
    if (window.matchMedia?.("(prefers-reduced-motion: reduce)").matches) { this.land(); return; }
    const pending = (this.pendingSettle = {});
    const off = v.addAnimator(() => {
      off();
      if (this.pendingSettle !== pending) return;
      this.pendingSettle = null;
      this.settle();
    });
    v.renderer.invalidate();
  }

  /** Finish a fade cut short over a moment, from where it was to the selection now. */
  settle() {
    const v = this.viewer;
    const cur = this.current();
    const step = this.fadeTo(this.levelsFor(cur), this.clipFor(cur));
    const blend = this.blend;
    const start = performance.now();
    let off = null;
    const stop = () => {
      off?.();
      off = null;
      if (this.stopSettle === stop) this.stopSettle = null;
      v.renderer.release("lasso-settle");
    };
    off = v.addAnimator((now) => {
      if (this.blend !== blend) { stop(); return; }
      const k = Math.min(Math.max((now - start) / SETTLE_MS, 0), 1);
      const e = k * k * (3 - 2 * k);
      step(e, e);
      if (k >= 1) this.land();
    });
    this.stopSettle = stop;
    v.renderer.hold("lasso-settle");
  }

  /** The live selection in the shape States keep (null when there is none). */
  current() {
    return this.selection ? { selection: this.selection, isolate: this.isolate, filterOff: !!this.filterOff, mask: this.mask } : null;
  }

  /**
   * Per trace, each object's lasso level for a selection `s` (a State's
   * lasso or `current()`): 255 shown, 128 dimmed, 0 hidden. Null when it
   * selects nothing (everything shown).
   */
  levelsFor(s) {
    if (!s?.selection || s.filterOff) return null;
    const v = this.viewer;
    const sel = s.selection instanceof Map ? s.selection : new Map(Object.entries(s.selection));
    const out = new Map();
    for (const batch of v.points.batches) {
      const l = new Uint8Array(batch.count).fill(s.isolate === false ? 128 : 0);
      for (const i of sel.get(batch.key) || []) if (i < l.length) l[i] = 255;
      out.set(batch.key, l);
    }
    return out;
  }

  /** The clip a selection puts on the volumes (null when none). */
  clipFor(s) {
    const v = this.viewer;
    if (!s?.selection || !s.mask || s.filterOff || v.state.global.lassoVolumes === false) return null;
    return { ...s.mask, outside: s.isolate === false ? 0.25 : 0 };
  }

  /**
   * A State change from one lasso selection to another plays as a fade
   * instead of a jump at arrival: what leaves the selection (objects, and
   * dust outside the new outline) fades out over the first three quarters
   * of the flight, what joins fades in over the last three quarters (an
   * overlapping crossfade, so the view never blinks), and the camera lands
   * as the new selection completes. The State's own apply then lands
   * exactly. The fade starts from what is drawn, so a change that
   * interrupts another carries on from where that one was.
   */
  transition(from, to) {
    // A classic selection (by name or outline) resolves at arrival; nothing to blend.
    if (to && !to.selection && (to.names?.length || to.mask)) return null;
    const levels = this.levelsFor(to), clip = this.clipFor(to);
    const same = JSON.stringify(from ?? null) === JSON.stringify(to ?? null);
    const none = !levels && !clip && !this.levelsFor(from) && !this.clipFor(from);
    if (!this.fading() && (same || none)) return null;
    // This change's fade carries on from a fade it cut short.
    this.pendingSettle = null;
    let step = null;
    return {
      // In real time (raw), not the camera's eased progress, which would
      // bunch both fades into the middle of the flight. The fade starts on
      // the first frame, so an instant change never starts one.
      step: (t, raw = t) => {
        step ||= this.fadeTo(levels, clip);
        const [out, inn] = lassoFade(raw);
        step(out, inn);
      },
    };
  }

  /** Clip (or dim) the volumes outside the lasso while its filter is on. */
  applyMask() {
    const v = this.viewer;
    v.volumes.setMask?.(this.clipFor(this.current()));
    v.renderer.invalidate();
  }

  count() {
    let n = 0;
    for (const idx of this.selection?.values() || []) n += idx.length;
    return n;
  }

  /** The visible volumes a lasso selects with the objects (none when lassos leave volumes alone). */
  lassoedVolumes() {
    const v = this.viewer;
    if (v.state.global.lassoVolumes === false) return [];
    const seen = new Set();
    return (v.manifest.volumes || []).filter((s) => v.state.volumes[s.stateKey]?.visible && !seen.has(s.stateKey) && seen.add(s.stateKey));
  }

  renderBar() {
    const bar = this.bar;
    bar.replaceChildren();
    bar.hidden = !this.selection;
    // Phones lift the key above the bar while it shows.
    this.ui.root.dataset.lassoSelection = String(!!this.selection);
    if (!this.selection) return;
    const n = this.count();
    const dust = this.mask ? this.lassoedVolumes() : [];
    const dustName = dust.length === 1 ? dust[0].name : `${dust.length} volumes`;
    const what = h("span", { class: "ov-selbar-what" },
      n ? h("span", null, `${n.toLocaleString()} selected`) : null,
      dust.length ? h("span", { class: "ov-selbar-dust" }, n ? `+ ${dustName}` : `${dustName} selected`) : null,
      n || dust.length ? null : h("span", null, "Lasso area selected"));
    // Frame and CSV need objects; a lasso around dust alone has none.
    bar.append(...[
      h("span", { class: "ov-selbar-count", title: dust.length ? "The dust outside the lasso is hidden with the other objects" : "" }, icon("wand"), what),
      n ? h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-label": "Frame the selection", onclick: () => this.frame() }, icon("target"), h("span", { class: "ov-btn-text" }, "Frame")) : null,
      h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-pressed": String(!this.isolate), title: this.isolate ? "Show the rest dimmed" : "Hide the rest", onclick: () => { this.isolate = !this.isolate; this.applyBits(); this.renderBar(); } }, h("span", { class: "ov-btn-text" }, "Dim others"), h("span", { class: "ov-btn-short" }, "Dim")),
      n ? h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-label": "Export the selection as CSV", onclick: () => this.exportCsv() }, icon("download"), h("span", { class: "ov-btn-text" }, "CSV")) : null,
      iconButton("close", "Clear selection", () => this.clear()),
    ].filter(Boolean));
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
    return { selection: out, isolate: this.isolate, filterOff: !!this.filterOff, mask: this.mask ? { polygon: this.mask.polygon, vp: this.mask.vp } : null };
  }

  restore(s) {
    const token = ++this.restoreToken;
    if (!s || (!s.selection && !s.names?.length && !s.mask)) {
      if (this.selection) this.clear();
      else if (this.blend) this.applyBits(); // a fade toward no selection still lands
      return;
    }
    if (!s.selection) {
      // A classic State: its clusters by name (or only the lasso outline).
      this.resolveLegacy(s).then((selection) => {
        if (token !== this.restoreToken) return;
        if (selection) this.restoreSelection({ ...s, selection });
        else if (this.selection) this.clear();
      });
      return;
    }
    this.restoreSelection(s);
  }

  restoreSelection(s) {
    this.selection = new Map(Object.entries(s.selection).map(([k, v]) => [k, Uint32Array.from(v)]));
    this.isolate = !!s.isolate;
    this.filterOff = !!s.filterOff;
    this.mask = s.mask && Array.isArray(s.mask.polygon) && Array.isArray(s.mask.vp) ? { polygon: s.mask.polygon, vp: s.mask.vp } : null;
    this.applyBits();
    this.renderBar();
  }

  /** What a classic lasso selected: its clusters by name, else the objects inside its outline. */
  async resolveLegacy(s) {
    const v = this.viewer;
    const out = {};
    let n = 0;
    const traces = v.traces.filter((t) => t.points && t.pickable !== false);
    if (s.names?.length) {
      const want = new Set(s.names.map(selectionKey));
      for (const t of traces) {
        const names = (await v.objectMeta(t.key))?.name || [];
        const hits = [];
        names.forEach((name, i) => { if (want.has(selectionKey(name))) hits.push(i); });
        if (hits.length) { out[t.key] = hits; n += hits.length; }
      }
    } else if (s.mask) {
      const { polygon, vp } = s.mask;
      for (const t of traces) {
        if (!v.state.traces[t.key]?.visible) continue;
        const hits = [];
        for (let i = 0; i < t.points.count; i++) {
          const q = projectNdc(v.objectPosition(t.key, i), vp);
          if (q && inside(q[0], q[1], polygon)) hits.push(i);
        }
        if (hits.length) { out[t.key] = hits; n += hits.length; }
      }
    }
    return n || s.mask ? out : null;
  }
}

/** A point through a column-major view-projection matrix, in normalised device coordinates. */
function projectNdc(p, m) {
  if (!p) return null;
  const x = m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12];
  const y = m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13];
  const w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
  return w > 1e-9 ? [x / w, y / w] : null;
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

/**
 * A State change's lasso fades at linear progress `raw` (0–1): what leaves
 * fades out over the first three quarters, what joins fades in over the
 * last three, each eased (smoothstep). Returns [fade-out, fade-in].
 */
export function lassoFade(raw) {
  const ease = (x) => { const t = Math.min(Math.max(x, 0), 1); return t * t * (3 - 2 * t); };
  return [ease(raw / 0.75), ease((raw - 0.25) / 0.75)];
}
