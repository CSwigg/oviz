// Annotate: draw labels, curves and arrows, bubbles and shells into the
// figure, in 3D. A quiet tool bar picks what to draw; a click places it (on
// a cluster it rides along with the cluster's orbit, on the dust it lands
// in the dust, elsewhere on the plane through the orbit centre). Selected
// annotations show handles to move, bend and resize them; Shift-drag moves
// along the line of sight. What you draw joins a group in the key on its
// own: a label put on a shell names it, a run of red arrows is one entry.
// Annotations are part of the figure (saved with it and in every view, and
// carried by view links), and view transitions morph them.

import { h, icon, iconButton, clear, localStorageGet, localStorageSet } from "./dom.js";
import { slider, toggle, miniSeg, numberInput, colorPicker, select } from "./controls.js";
import { figureId } from "./story.js";
import { AnnotationsLayer, ANNOT_SEG_FLOATS, ANNOT_ARROW_FLOATS, ANNOT_SPHERE_FLOATS } from "../layers/annotations.js";
import {
  ANNOT_DEFAULTS, emptyAnnotations, normalizeAnnotations, notesToAnnotations, resolveAnnotPoint,
  annotationsTimeDependent, curveSamples, sphereRotation, sphereAxes, rotationQuat, polylineDistance,
  annotKind, annotKindName, annotGroupName, annotGroupLook, assignGroup, lerpAnnotations, annotId,
} from "./annotate-model.js";
import { parseColor } from "../core/color.js";
import { presentDayFade } from "../app/timeline.js";
import { cloneJson } from "../app/state.js";

const TOOLS = [
  { id: "select", icon: "pointer", label: "Select and edit", hint: "Click an annotation to edit it · drag empty space to orbit" },
  { id: "text", icon: "text", label: "Label", hint: "Click to place a label · on a cluster it follows its orbit" },
  { id: "curve", icon: "curve", label: "Curve", hint: "Click to add points · double-click or Enter to finish" },
  { id: "arrow", icon: "arrow", label: "Arrow", hint: "Drag from where the arrow starts to where it points" },
  { id: "bubble", icon: "bubble", label: "Bubble", hint: "Drag out from the centre to size a bubble" },
  { id: "shell", icon: "shell", label: "Shell", hint: "Drag out from the centre to size a shell" },
];
const ACCENT = [1, 0.82, 0.48];
const HISTORY = 120;
const AXIS_COLORS = ["#ff8a80", "#8fe39b", "#82b9ff"];

export class AnnotatePlugin {
  constructor() {
    this.name = "annotate";
    this.doc = emptyAnnotations();
    this.view = null; // the document part-way through a State transition
    this.active = false;
    this.tool = "select";
    this.sel = null; // { id, point }
    this.hoverId = null;
    this.history = [];
    this.future = [];
    this.version = 0;
    this.sceneVersion = 0;
    this.geom = new Map();
    this.drawing = null; // a curve being clicked out
    this.gesture = null; // a press–drag that creates (arrow, bubble, shell)
    this.drag = null; // a handle or body being dragged
    this.editing = null; // the inline label editor
    this.styles = {}; // the last look used per tool
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    ui.annotate = this;
    this.layer = v.renderer.add(new AnnotationsLayer(v.gl));
    v.on("gpu-restored", () => {
      this.layer = v.renderer.add(new AnnotationsLayer(v.gl));
      this._sceneKey = "";
      v.renderer.invalidate();
    });
    v.renderer.beforeRender.push(() => this.frame());
    v.renderer.overlays.push({ update: (frame) => this.updateOverlay(frame) });
    v.stateExtensions.set("annotations", {
      capture: () => cloneJson(this.doc),
      // A view saved before annotations existed leaves them as they are.
      apply: (d) => {
        this.view = null;
        if (d && typeof d === "object") this.setDoc(normalizeAnnotations(d), { history: false, keepSelection: true });
        else this.refresh();
      },
      transition: (a, b) => this.transition(a, b),
    });
    // Notes from earlier figures and views land as labels.
    v.stateExtensions.set("notes", { capture: () => undefined, apply: (n) => { if (Array.isArray(n) && n.length) this.mergeDoc(notesToAnnotations(n)); } });
    let doc = normalizeAnnotations(v.manifest.annotations);
    const legacy = v.manifest.initialState?.manual_labels;
    if (Array.isArray(legacy) && legacy.length) doc = mergeInto(doc, notesToAnnotations(legacy.map(legacyLabel).filter(Boolean)));
    this.doc = doc;
    this.base = JSON.stringify(doc);
    this.buildBar();
    this.buildPanel();
    this.buildOverlay();
    this.button = iconButton("pen", "Annotate", () => this.setActive(!this.active), { cls: "ov-annot-btn", shortcut: "K" });
    this.button.setAttribute("aria-pressed", "false");
    ui.addToolbarButton?.(this.button);
    this.bindPointer();
    v.on("time", () => { if (this.timeDependent) v.renderer.invalidate(); });
    v.on("viewmode", () => { this._sceneKey = ""; });
    v.on("state-applied", () => this.renderPanel());
    // Drafts wait for the Story plugin, which says whether this figure keeps them.
    setTimeout(() => this.restoreDraft(), 0);
    this.refresh();
  }

  // ------------------------------------------------------------ document

  get drawn() {
    return this.view || this.doc;
  }

  item(id) {
    return this.doc.items.find((i) => i.id === id) || null;
  }

  group(id) {
    return this.doc.groups.find((g) => g.id === id) || null;
  }

  /** Replace the document (from a State, a link, a draft or undo). */
  setDoc(doc, { history = true, keepSelection = false } = {}) {
    if (history) this.pushHistory();
    this.doc = doc;
    if (!keepSelection || (this.sel && !this.item(this.sel.id))) this.sel = null;
    this.changed({ history: false });
  }

  /** Add another document's groups and items (skipping ids already here). */
  mergeDoc(extra) {
    const merged = mergeInto(cloneJson(this.doc), extra);
    if (merged.items.length === this.doc.items.length) return;
    this.doc = merged;
    this.changed({ history: false });
  }

  pushHistory() {
    this.history.push(JSON.stringify(this.doc));
    if (this.history.length > HISTORY) this.history.shift();
    this.future = [];
  }

  /** Call after every edit: redraws, refreshes the key, panel and draft. */
  changed({ history = true } = {}) {
    if (history) this.pushHistory();
    this.version++;
    this.timeDependent = annotationsTimeDependent(this.doc);
    this._sceneKey = "";
    this.refresh();
    this.saveDraft();
  }

  refresh() {
    this.viewer.renderer.invalidate();
    this.ui.layers?.render();
    this.ui.legend?.render();
    this.renderPanel();
    this.syncBar();
  }

  undo() {
    if (!this.history.length) { this.ui.toast("Nothing to undo"); return; }
    this.future.push(JSON.stringify(this.doc));
    this.doc = JSON.parse(this.history.pop());
    if (this.sel && !this.item(this.sel.id)) this.sel = null;
    this.changed({ history: false });
  }

  redo() {
    if (!this.future.length) return;
    this.history.push(JSON.stringify(this.doc));
    this.doc = JSON.parse(this.future.pop());
    if (this.sel && !this.item(this.sel.id)) this.sel = null;
    this.changed({ history: false });
  }

  /** Add a new item, grouped by the smart rules, and select it. */
  addItem(item, { select = true } = {}) {
    this.pushHistory();
    item.id = annotId("i");
    item.group = assignGroup(this.doc, item, {
      resolve: (p) => this.resolve(p),
      screenDistance: (curve, pos) => this.curveScreenDistance(curve, pos),
    });
    this.doc.items.push(item);
    if (select) this.sel = { id: item.id, point: null };
    this.changed({ history: false });
    return item;
  }

  removeItem(id) {
    const it = this.item(id);
    if (!it) return;
    this.pushHistory();
    this.doc.items = this.doc.items.filter((i) => i.id !== id);
    const used = new Set(this.doc.items.map((i) => i.group));
    this.doc.groups = this.doc.groups.filter((g) => used.has(g.id));
    for (const g of this.doc.groups) if (g.nameFrom === id) delete g.nameFrom;
    if (this.sel?.id === id) this.sel = null;
    this.changed({ history: false });
  }

  duplicate(id) {
    const it = this.item(id);
    if (!it) return;
    const copy = cloneJson(it);
    // Offset the copy by a few pixels' worth, so it shows as a new one.
    const cam = this.viewer.renderer.camera;
    const ref = this.anchorOf(it);
    const s = ref ? cam.pixelScale(Math.max(1e-6, dist3(ref, cam.eye))) * 24 : 0;
    const d = [cam.right[0] * s, cam.right[1] * s, cam.right[2] * s];
    const shift = (p) => (p.pos ? { pos: add3(p.pos, d) } : { ...p, offset: add3(p.offset || [0, 0, 0], d) });
    if (copy.kind === "text") copy.at = shift(copy.at);
    else if (copy.kind === "sphere") copy.center = shift(copy.center);
    else copy.points = copy.points.map(shift);
    this.pushHistory();
    copy.id = annotId("i");
    this.doc.items.push(copy);
    this.sel = { id: copy.id, point: null };
    this.changed({ history: false });
  }

  /** Change an item; `commit` records an undo step (the end of a slider drag). */
  patch(id, fields, { commit = true } = {}) {
    const it = this.item(id);
    if (!it) return;
    if (commit) this.pushHistory();
    Object.assign(it, fields);
    // The next one drawn looks like this one.
    const tool = it.kind === "text" ? "text" : it.kind === "sphere" ? it.style === "shell" ? "shell" : "bubble" : it.arrow !== "none" ? "arrow" : "curve";
    const keep = ["color", "opacity", "width", "dash", "size", "weight", "bg", "style"];
    this.styles[tool] = { ...(this.styles[tool] || {}), ...Object.fromEntries(Object.entries(fields).filter(([k]) => keep.includes(k))) };
    this.version++;
    this._sceneKey = "";
    this.viewer.renderer.invalidate();
    if (commit) {
      this.ui.layers?.render();
      this.ui.legend?.render();
      this.saveDraft();
    }
  }

  patchGroup(id, fields) {
    const g = this.group(id);
    if (!g) return;
    this.pushHistory();
    Object.assign(g, fields);
    this.changed({ history: false });
  }

  // ------------------------------------------------------------ geometry

  resolve(p) {
    const v = this.viewer;
    return resolveAnnotPoint(p, (t, i) => v.objectPosition(t, i));
  }

  /** A world point that stands for an item (its anchor, centre or middle). */
  anchorOf(it) {
    if (it.kind === "text") return this.resolve(it.at);
    if (it.kind === "sphere") return this.resolve(it.center);
    const pts = it.points.map((p) => this.resolve(p)).filter(Boolean);
    if (!pts.length) return null;
    const c = [0, 0, 0];
    for (const p of pts) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
    return [c[0] / pts.length, c[1] / pts.length, c[2] / pts.length];
  }

  project(p) {
    const s = p && this.viewer.renderer.camera.project(p, [0, 0, 0]);
    return s ? [s[0], s[1]] : null;
  }

  /** Pixel distance on screen from world point `pos` to a curve. */
  curveScreenDistance(curve, pos) {
    const at = this.project(pos);
    if (!at) return Infinity;
    const pts = curve.points.map((p) => this.resolve(p));
    if (pts.some((p) => !p)) return Infinity;
    const s = curveSamples(pts, { smooth: curve.smooth, closed: curve.closed, perSegment: 8 }).points.map((p) => this.project(p));
    return polylineDistance(at[0], at[1], s);
  }

  // ------------------------------------------------------------ drawing

  frame() {
    const v = this.viewer;
    const doc = this.drawn;
    const key = `${this.version}|${this.view ? this.sceneVersion : 0}|${this.timeDependent || this.view ? v.timeline.frame : 0}|${v.state.view.mode}|${this.sel?.id || ""}|${this.hoverId || ""}|${this.editing?.id || ""}|${this.previewKey || ""}`;
    if (key === this._sceneKey && this.layer.scene) return;
    this._sceneKey = key;
    this.layer.setScene(this.buildScene(doc));
  }

  buildScene(doc) {
    const v = this.viewer;
    const time = v.timeline.time;
    const groupOn = new Map(doc.groups.map((g) => [g.id, g.visible !== false]));
    const segs = [], arrows = [], spheres = [], texts = [];
    this.geom = new Map();
    const pushCurve = (id, pts, it, alpha, selected) => {
      const { points: sp, arc } = curveSamples(pts, { smooth: it.smooth, closed: it.closed });
      const total = arc[arc.length - 1] || 1;
      const dash = it.dash === "dash" ? [total / 30, total / 46] : it.dash === "dot" ? [total / 220, total / 70] : [0, 0];
      const c = parseColor(it.color);
      const width = it.width;
      if (selected || id === this.hoverId) {
        // A soft accent underlay marks what is selected (or under the pointer).
        const a = selected ? 0.38 : 0.2;
        for (let i = 1; i < sp.length; i++) segs.push(...sp[i - 1], ...sp[i], arc[i - 1], arc[i], ...ACCENT, a, width + 7, 0, 0, 0);
      }
      for (let i = 1; i < sp.length; i++) segs.push(...sp[i - 1], ...sp[i], arc[i - 1], arc[i], c[0], c[1], c[2], alpha, width, dash[0], dash[1], 0);
      const head = Math.max(14, width * 5);
      const back = (end) => {
        // The arrowhead points along the last stretch at least a head long on screen.
        const tipS = this.project(sp[end]);
        const step = end === 0 ? 1 : -1;
        let k = end + step;
        while (k > 0 && k < sp.length - 1 && tipS) {
          const s = this.project(sp[k]);
          if (s && Math.hypot(s[0] - tipS[0], s[1] - tipS[1]) >= head * 0.8) break;
          k += step;
        }
        return sp[Math.max(0, Math.min(sp.length - 1, k))];
      };
      if (!it.closed && (it.arrow === "end" || it.arrow === "both")) arrows.push(...sp[sp.length - 1], ...back(sp.length - 1), c[0], c[1], c[2], alpha, head, width);
      if (!it.closed && (it.arrow === "start" || it.arrow === "both")) arrows.push(...sp[0], ...back(0), c[0], c[1], c[2], alpha, head, width);
      this.geom.set(id, { kind: "curve", samples: sp, points: pts });
    };
    for (const it of doc.items) {
      if (it.visible === false || !groupOn.get(it.group)) continue;
      let alpha = (it.opacity ?? 1) * (it.fade ?? 1);
      if (it.present) alpha *= presentDayFade(time);
      if (alpha <= 0.002) continue;
      const selected = this.sel?.id === it.id && !this.view;
      if (it.kind === "curve") {
        const pts = it.points.map((p) => this.resolve(p));
        if (pts.some((p) => !p)) continue;
        pushCurve(it.id, pts, it, alpha, selected);
      } else if (it.kind === "sphere") {
        const c3 = this.resolve(it.center);
        if (!c3) continue;
        const q = rotationQuat(sphereRotation(it.rot));
        const c = parseColor(it.color);
        const style = it.style === "shell" ? 1 : it.style === "wire" ? 2 : 0;
        spheres.push(...c3, ...it.radii, ...q, c[0], c[1], c[2], alpha, style, selected || it.id === this.hoverId ? 1 : 0);
        this.geom.set(it.id, { kind: "sphere", center: c3, axes: sphereAxes(it) });
      } else {
        const at = this.resolve(it.at);
        if (!at) continue;
        this.geom.set(it.id, { kind: "text", at });
        if (this.editing?.id === it.id) continue;
        texts.push({ id: it.id, pos: at, text: it.text, size: it.size, weight: it.weight, color: it.color, bg: it.bg, dx: it.dx, dy: it.dy, opacity: alpha });
      }
    }
    // What is being drawn right now.
    const pv = this.preview();
    if (pv?.kind === "curve") pushCurve("__preview", pv.points, pv, 0.85, false);
    if (pv?.kind === "sphere") {
      const c = parseColor(pv.color);
      spheres.push(...pv.center, ...pv.radii, 0, 0, 0, 1, c[0], c[1], c[2], pv.opacity, pv.style === "shell" ? 1 : 0, 1);
    }
    return {
      version: ++this.sceneVersion,
      segments: new Float32Array(segs), segmentCount: segs.length / ANNOT_SEG_FLOATS,
      arrows: new Float32Array(arrows), arrowCount: arrows.length / ANNOT_ARROW_FLOATS,
      spheres: new Float32Array(spheres), sphereCount: spheres.length / ANNOT_SPHERE_FLOATS,
      texts,
    };
  }

  /** The shape a drawing gesture would make now (resolved positions), or null. */
  preview() {
    if (this.drawing && this.drawing.points.length) {
      const pts = this.drawing.points.map((p) => this.resolve(p));
      if (this.drawing.hover) pts.push(this.drawing.hover);
      if (pts.length >= 2 && pts.every(Boolean)) return { kind: "curve", points: pts, ...this.look("curve") };
    }
    const g = this.gesture;
    if (g?.end && g.tool === "arrow") return { kind: "curve", points: [g.startPos, g.end], ...this.look("arrow") };
    if (g?.radius > 0) return { kind: "sphere", center: g.startPos, radii: [g.radius, g.radius, g.radius], ...this.look(g.tool) };
    return null;
  }

  /** The look the next item of `tool` gets: its defaults and the last one edited. */
  look(tool) {
    return { ...ANNOT_DEFAULTS[tool], ...(this.styles[tool] || {}) };
  }

  // ------------------------------------------------------------ picking

  /** The annotation under (x, y) CSS px: labels first, then curves, then shells. */
  hitTest(x, y) {
    const doc = this.doc;
    const boxes = this.layer.textBoxes();
    const groupOn = new Map(doc.groups.map((g) => [g.id, g.visible !== false]));
    const shown = (it) => it.visible !== false && groupOn.get(it.group) && this.geom.has(it.id);
    for (let i = doc.items.length - 1; i >= 0; i--) {
      const it = doc.items[i];
      if (it.kind !== "text" || !shown(it)) continue;
      const b = boxes.get(it.id);
      if (b && x >= b.x - 3 && x <= b.x + b.w + 3 && y >= b.y - 3 && y <= b.y + b.h + 3) return { id: it.id };
    }
    let best = null;
    for (const it of doc.items) {
      if (it.kind !== "curve" || !shown(it)) continue;
      const g = this.geom.get(it.id);
      const d = polylineDistance(x, y, g.samples.map((p) => this.project(p)));
      if (d <= Math.max(7, it.width / 2 + 5) && (!best || d < best.d)) best = { id: it.id, d };
    }
    if (best) return best;
    let inside = null;
    for (const it of doc.items) {
      if (it.kind !== "sphere" || !shown(it)) continue;
      const s = this.sphereScreen(it);
      if (!s) continue;
      const d = Math.hypot(x - s.c[0], y - s.c[1]);
      if (Math.abs(d - s.r) <= 9) return { id: it.id };
      if (d < s.r && (!inside || s.r < inside.r)) inside = { id: it.id, r: s.r };
    }
    return inside;
  }

  /** Centre and approximate silhouette radius of a shell on screen. */
  sphereScreen(it) {
    const g = this.geom.get(it.id);
    if (!g) return null;
    const c = this.project(g.center);
    if (!c) return null;
    const cam = this.viewer.renderer.camera;
    const r = Math.max(...it.radii);
    const e = this.project(add3(g.center, [cam.right[0] * r, cam.right[1] * r, cam.right[2] * r]));
    return e ? { c, r: Math.hypot(e[0] - c[0], e[1] - c[1]) } : null;
  }

  /** Where a click lands: on a cluster (riding with it), the dust, or the plane through the orbit centre. */
  placeAt(x, y, e, { plane = null } = {}) {
    const v = this.viewer;
    if (!e?.altKey) {
      const hit = v.pick(x, y);
      if (hit && hit.kind !== "member" && hit.trace != null) {
        const pos = v.objectPosition(hit.trace, hit.index);
        if (pos) return { point: { trace: hit.trace, index: hit.index }, pos, object: hit };
      }
      if (!plane) {
        const surf = v.pickSurface(x, y);
        if (surf) return { point: { pos: surf }, pos: surf };
      }
    }
    const pos = plane ? this.onPlane(x, y, plane) : this.focusPoint(x, y);
    return pos ? { point: { pos }, pos } : null;
  }

  /** The point under (x, y) on the plane through `through` facing the camera. */
  onPlane(x, y, through) {
    const cam = this.viewer.renderer.camera;
    const { origin: o, dir: d } = cam.ray(x, y);
    const f = cam.forward;
    const den = d[0] * f[0] + d[1] * f[1] + d[2] * f[2];
    if (Math.abs(den) < 1e-9) return null;
    const t = ((through[0] - o[0]) * f[0] + (through[1] - o[1]) * f[1] + (through[2] - o[2]) * f[2]) / den;
    return t > 0 ? [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t] : null;
  }

  focusPoint(x, y) {
    const v = this.viewer;
    const cam = v.renderer.camera;
    if (v.pose.distance > 1e-3) return this.onPlane(x, y, v.pose.target);
    // Sky view looks out from the Sun: place at a typical cluster distance.
    const { origin: o, dir: d } = cam.ray(x, y);
    return [o[0] + d[0] * 300, o[1] + d[1] * 300, o[2] + d[2] * 300];
  }

  // ------------------------------------------------------------ pointer

  bindPointer() {
    const c = this.viewer.canvas;
    const xy = (e) => { const r = c.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    const swallow = (e) => { e.stopImmediatePropagation(); e.preventDefault(); };
    c.addEventListener("pointerdown", (e) => {
      if (!this.active || e.button !== 0) return;
      const [x, y] = xy(e);
      this.commitEditor();
      if (this.tool === "select") {
        const hit = this.hitTest(x, y);
        if (hit) {
          this.select(hit.id);
          this.startDrag(e, { type: "move", id: hit.id }, x, y, c);
          swallow(e);
        } else if (this.sel) {
          this.select(null);
        }
        return;
      }
      if (this.tool === "arrow" || this.tool === "bubble" || this.tool === "shell") {
        const place = this.placeAt(x, y, e);
        if (!place) return;
        this.gesture = { tool: this.tool, start: place.point, startPos: place.pos, x0: x, y0: y, radius: 0, end: null };
        try { c.setPointerCapture(e.pointerId); } catch (_) { /* synthetic */ }
        swallow(e);
        return;
      }
      // Labels and curves are placed by clicks; a drag still orbits.
      this.pending = { x, y };
    }, true);
    c.addEventListener("pointermove", (e) => {
      if (!this.active) return;
      const [x, y] = xy(e);
      if (this.drag) { this.updateDrag(x, y, e); swallow(e); return; }
      if (this.gesture) { this.updateGesture(x, y, e); swallow(e); return; }
      if (this.pending && Math.hypot(x - this.pending.x, y - this.pending.y) > 5) this.pending = null;
      this.hover(x, y, e);
    }, true);
    c.addEventListener("pointerup", (e) => {
      if (!this.active) return;
      const [x, y] = xy(e);
      if (this.drag) { this.endDrag(); swallow(e); return; }
      if (this.gesture) { this.finishGesture(x, y, e); swallow(e); return; }
      if (this.pending && e.button === 0) {
        this.pending = null;
        if (this.tool === "text") { this.placeText(x, y, e); swallow(e); }
        else if (this.tool === "curve") { this.addCurvePoint(x, y, e); swallow(e); }
      }
    }, true);
    c.addEventListener("dblclick", (e) => {
      if (!this.active) return;
      swallow(e);
      if (this.tool === "curve") { this.finishCurve(); return; }
      const [x, y] = xy(e);
      const hit = this.hitTest(x, y);
      if (hit && this.item(hit.id)?.kind === "text") this.openEditor(hit.id);
    }, true);
    c.addEventListener("pointerleave", () => { if (this.drawing) { this.drawing.hover = null; this.bumpPreview(); } });
  }

  hover(x, y, e) {
    const c = this.viewer.canvas;
    if (this.tool === "select") {
      const hit = this.hitTest(x, y);
      const id = hit?.id || null;
      c.style.cursor = id ? "move" : "";
      if (id !== this.hoverId) { this.hoverId = id; this.viewer.renderer.invalidate(); }
      return;
    }
    c.style.cursor = "crosshair";
    if (this.drawing) {
      const last = this.resolve(this.drawing.points[this.drawing.points.length - 1]);
      const place = this.placeAt(x, y, e, { plane: last });
      this.drawing.hover = place?.pos || null;
      this.bumpPreview();
    }
    // A ring marks the cluster a click would attach to.
    const hit = !e.altKey ? this.viewer.pick(x, y) : null;
    this.snapTarget = hit && hit.kind !== "member" && hit.trace != null ? hit : null;
    this.overlayKey = "";
    this.viewer.renderer.invalidate();
  }

  bumpPreview() {
    this.previewKey = (this.previewKey || 0) + 1;
    this.viewer.renderer.invalidate();
  }

  async placeText(x, y, e) {
    const place = this.placeAt(x, y, e);
    if (!place) return;
    const look = this.look("text");
    const onObject = place.point.trace != null;
    // A label on a cluster starts as the cluster's name (names load on demand).
    if (onObject) await this.viewer.objectMeta(place.point.trace).catch(() => null);
    const name = onObject ? this.viewer.describeObject(place.point.trace, place.point.index)?.name : "";
    const it = this.addItem({ kind: "text", at: place.point, text: name || "Label", ...look, dy: onObject ? -16 : 0 });
    this.openEditor(it.id, { fresh: true });
  }

  addCurvePoint(x, y, e) {
    const prev = this.drawing?.points.length ? this.resolve(this.drawing.points[this.drawing.points.length - 1]) : null;
    const place = this.placeAt(x, y, e, { plane: prev });
    if (!place) return;
    if (!this.drawing) this.drawing = { points: [], hover: null };
    const pts = this.drawing.points;
    // The second click of a double-click (which finishes) adds no point.
    const last = pts.length ? this.project(this.resolve(pts[pts.length - 1])) : null;
    if (last && Math.hypot(last[0] - x, last[1] - y) < 5) return;
    // A click on the first point closes the loop.
    if (pts.length >= 3) {
      const first = this.project(this.resolve(pts[0]));
      if (first && Math.hypot(first[0] - x, first[1] - y) < 10) { this.finishCurve({ closed: true }); return; }
    }
    pts.push(place.point);
    this.bumpPreview();
    if (pts.length === 1) this.ui.toast("Keep clicking to add points · double-click or Enter to finish", { ms: 2200 });
  }

  finishCurve({ closed = false } = {}) {
    const d = this.drawing;
    this.drawing = null;
    this.bumpPreview();
    if (!d || d.points.length < 2) return;
    this.addItem({ kind: "curve", points: d.points, ...this.look("curve"), closed: closed && d.points.length > 2 });
  }

  updateGesture(x, y, e) {
    const g = this.gesture;
    if (g.tool === "arrow") {
      const place = this.placeAt(x, y, e, { plane: g.startPos });
      g.end = place?.pos || null;
      g.endPoint = place?.point || null;
    } else {
      const p = this.onPlane(x, y, g.startPos);
      g.radius = p ? dist3(p, g.startPos) : 0;
    }
    this.bumpPreview();
  }

  finishGesture(x, y, e) {
    const g = this.gesture;
    this.gesture = null;
    this.bumpPreview();
    const moved = Math.hypot(x - g.x0, y - g.y0);
    if (g.tool === "arrow") {
      if (moved < 6 || !g.endPoint) { this.ui.toast("Drag from where the arrow starts to where it points", { ms: 2000 }); return; }
      this.addItem({ kind: "curve", points: [g.start, g.endPoint], ...this.look("arrow") });
      return;
    }
    let r = g.radius;
    if (moved < 4 || !(r > 0)) {
      // A click: a shell about 80 px across.
      const cam = this.viewer.renderer.camera;
      r = cam.pixelScale(Math.max(1e-6, dist3(g.startPos, cam.eye))) * 40;
    }
    this.addItem({ kind: "sphere", center: g.start, radii: [r, r, r], ...this.look(g.tool) });
  }

  // ------------------------------------------------------------ dragging

  startDrag(e, d, x, y, target) {
    const it = this.item(d.id);
    if (!it) return;
    const ref = d.ref || this.anchorOf(it);
    if (!ref) return;
    const plane = ref.slice();
    this.drag = { ...d, x0: x, y0: y, ref: plane, start: this.onPlane(x, y, plane) || plane, before: JSON.stringify(this.doc), origin: cloneJson(it), moved: false };
    try { target.setPointerCapture(e.pointerId); } catch (_) { /* synthetic */ }
  }

  /** World displacement for a drag to (x, y): across the screen, or along the line of sight with Shift. */
  dragDelta(x, y, e) {
    const d = this.drag;
    const cam = this.viewer.renderer.camera;
    if (e.shiftKey) {
      const s = cam.pixelScale(Math.max(1e-6, dist3(d.ref, cam.eye)));
      const k = -(y - d.y0) * s * 2;
      return [cam.forward[0] * k, cam.forward[1] * k, cam.forward[2] * k];
    }
    const p = this.onPlane(x, y, d.ref);
    return p ? sub3(p, d.start) : [0, 0, 0];
  }

  updateDrag(x, y, e) {
    const d = this.drag;
    const it = this.item(d.id);
    if (!it) return;
    if (!d.moved && Math.hypot(x - d.x0, y - d.y0) < 2) return;
    d.moved = true;
    const o = d.origin;
    const move = (p, delta) => (p.pos ? { pos: add3(p.pos, delta) } : { ...p, offset: add3(p.offset || [0, 0, 0], delta) });
    if (d.type === "move") {
      const delta = this.dragDelta(x, y, e);
      if (it.kind === "text") it.at = move(o.at, delta);
      else if (it.kind === "sphere") it.center = move(o.center, delta);
      else it.points = o.points.map((p) => move(p, delta));
    } else if (d.type === "point") {
      // ⌘ / Ctrl snaps the point onto a cluster or the dust under the pointer.
      if (e.metaKey || e.ctrlKey) {
        const place = this.placeAt(x, y, e);
        if (place) it.points = o.points.map((p, k) => (k === d.index ? place.point : p));
      } else {
        const delta = this.dragDelta(x, y, e);
        it.points = o.points.map((p, k) => (k === d.index ? move(p, delta) : p));
      }
    } else if (d.type === "radius") {
      const c = this.resolve(it.center);
      const p = c && this.onPlane(x, y, d.ref);
      if (c && p) {
        const k = Math.max(1e-3, dist3(p, c)) / Math.max(...o.radii);
        it.radii = o.radii.map((r) => Math.max(1e-3, r * k));
      }
    } else if (d.type === "axis") {
      const c = this.resolve(it.center);
      const p = c && this.onPlane(x, y, d.ref);
      if (c && p) {
        const ax = sphereAxes(o)[d.axis];
        const len = Math.hypot(...ax) || 1;
        const u = [ax[0] / len, ax[1] / len, ax[2] / len];
        const r = Math.abs(dot3(sub3(p, c), u));
        it.radii = o.radii.map((x0, k) => (k === d.axis ? Math.max(1e-3, r) : x0));
      }
    } else if (d.type === "size") {
      const anchor = this.project(this.resolve(it.at));
      if (anchor) {
        const r0 = Math.max(4, Math.hypot(d.x0 - anchor[0], d.y0 - anchor[1]));
        const r1 = Math.hypot(x - anchor[0], y - anchor[1]);
        it.size = Math.min(96, Math.max(6, Math.round(o.size * (r1 / r0))));
      }
    }
    this.version++;
    this._sceneKey = "";
    this.viewer.renderer.invalidate();
  }

  endDrag() {
    const d = this.drag;
    this.drag = null;
    if (!d) return;
    if (d.moved) {
      this.history.push(d.before);
      if (this.history.length > HISTORY) this.history.shift();
      this.future = [];
      this.changed({ history: false });
    }
  }

  // ------------------------------------------------------------ handles

  buildOverlay() {
    const ns = "http://www.w3.org/2000/svg";
    this.svg = document.createElementNS(ns, "svg");
    this.svg.setAttribute("class", "ov-annot-handles");
    this.svg.setAttribute("aria-hidden", "true");
    this.viewer.container.appendChild(this.svg);
    this.svg.addEventListener("pointerdown", (e) => {
      const el = e.target.closest?.("[data-handle]");
      if (!el || !this.sel) return;
      e.stopPropagation();
      e.preventDefault();
      const r = this.viewer.canvas.getBoundingClientRect();
      const x = e.clientX - r.left, y = e.clientY - r.top;
      const it = this.item(this.sel.id);
      if (!it) return;
      const kind = el.dataset.handle;
      if (kind === "insert") {
        // Bend a curve: a new point where the "+" sits, dragged at once.
        const k = Number(el.dataset.index);
        const pos = this.geom.get(it.id)?.mids?.[k];
        if (!pos) return;
        this.pushHistory();
        it.points.splice(k + 1, 0, { pos: pos.slice() });
        this.sel.point = k + 1;
        this.version++;
        this.startDrag(e, { type: "point", id: it.id, index: k + 1, ref: pos.slice() }, x, y, el);
        this.drag.before = this.history.pop(); // one undo step for insert + drag
        this.drag.moved = true;
        this.renderPanel();
        return;
      }
      const index = Number(el.dataset.index);
      if (kind === "point") {
        this.sel.point = index;
        this.renderPanel();
        this.startDrag(e, { type: "point", id: it.id, index, ref: this.resolve(it.points[index]) }, x, y, el);
      } else if (kind === "axis") {
        const c = this.resolve(it.center);
        this.startDrag(e, { type: "axis", id: it.id, axis: index, ref: c && add3(c, sphereAxes(it)[index]) }, x, y, el);
      } else {
        this.startDrag(e, { type: kind, id: it.id }, x, y, el);
      }
    });
    this.svg.addEventListener("pointermove", (e) => {
      if (!this.drag) return;
      const r = this.viewer.canvas.getBoundingClientRect();
      this.updateDrag(e.clientX - r.left, e.clientY - r.top, e);
    });
    const end = () => { if (this.drag) this.endDrag(); };
    this.svg.addEventListener("pointerup", end);
    this.svg.addEventListener("pointercancel", end);
    this.svg.addEventListener("dblclick", (e) => {
      const el = e.target.closest?.("[data-handle='point']");
      if (!el || !this.sel) return;
      e.stopPropagation();
      const it = this.item(this.sel.id);
      if (it?.kind === "curve" && it.points.length > 2) {
        this.pushHistory();
        it.points.splice(Number(el.dataset.index), 1);
        if (it.points.length < 3) it.closed = false;
        this.sel.point = null;
        this.changed({ history: false });
      }
    });
  }

  updateOverlay(frame) {
    const cam = frame.camera;
    const it = this.active && this.sel && !this.view ? this.item(this.sel.id) : null;
    const key = `${cam.version}|${this.sceneVersion}|${this.active}|${this.sel?.id}|${this.sel?.point}|${this.snapTarget?.trace}:${this.snapTarget?.index}|${this.tool}`;
    if (key === this.overlayKey) return;
    this.overlayKey = key;
    const ns = "http://www.w3.org/2000/svg";
    const svg = this.svg;
    while (svg.firstChild) svg.firstChild.remove();
    const add = (tag, attrs) => {
      const el = document.createElementNS(ns, tag);
      for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
      svg.appendChild(el);
      return el;
    };
    if (this.active && this.tool !== "select" && this.snapTarget) {
      const p = this.project(this.viewer.objectPosition(this.snapTarget.trace, this.snapTarget.index));
      if (p) add("circle", { cx: p[0], cy: p[1], r: 11, class: "ov-annot-snap" });
    }
    if (!it) return;
    const g = this.geom.get(it.id);
    if (!g) return;
    const handle = (p, kind, index, extra = {}) => {
      const s = this.project(p);
      if (!s) return null;
      return add("circle", { cx: s[0], cy: s[1], r: extra.r || 6, class: `ov-annot-handle ${extra.cls || ""}`, "data-handle": kind, "data-index": index ?? "", ...(extra.fill ? { style: `fill:${extra.fill}` } : {}) });
    };
    if (it.kind === "curve") {
      // "+" between points bends the curve; the points themselves move.
      const mids = [];
      const n = it.points.length;
      const per = it.smooth > 0 ? 24 : 1;
      for (let k = 0; k < (it.closed ? n : n - 1); k++) {
        const s = g.samples;
        const mid = per > 1 ? s[Math.min(s.length - 1, k * per + per / 2)] : g.points[k].map((x, j) => (x + g.points[(k + 1) % n][j]) / 2);
        mids.push(mid);
        const el = handle(mid, "insert", k, { r: 5, cls: "ov-annot-handle--insert" });
        if (el) el.appendChild(document.createElementNS(ns, "title")).textContent = "Drag to bend";
      }
      g.mids = mids;
      g.points.forEach((p, k) => handle(p, "point", k, { cls: this.sel.point === k ? "is-active" : "" }));
    } else if (it.kind === "sphere") {
      handle(g.center, "move", null, { r: 6 });
      const s = this.sphereScreen(it);
      if (s) {
        const cam2 = this.viewer.renderer.camera;
        const r = Math.max(...it.radii);
        handle(add3(g.center, [cam2.right[0] * r, cam2.right[1] * r, cam2.right[2] * r]), "radius", null, { r: 6, cls: "ov-annot-handle--radius" });
      }
      g.axes.forEach((ax, k) => handle(add3(g.center, ax), "axis", k, { r: 4.5, cls: "ov-annot-handle--axis", fill: AXIS_COLORS[k] }));
    } else {
      const b = this.layer.textBoxes().get(it.id);
      if (b) {
        add("rect", { x: b.x - 4, y: b.y - 3, width: b.w + 8, height: b.h + 6, rx: 6, class: "ov-annot-box" });
        add("circle", { cx: b.x + b.w + 4, cy: b.y + b.h + 3, r: 5, class: "ov-annot-handle ov-annot-handle--size", "data-handle": "size", "data-index": "" });
      }
      const a = this.project(g.at);
      if (a && (it.dx || it.dy)) add("circle", { cx: a[0], cy: a[1], r: 3, class: "ov-annot-anchor" });
    }
  }

  // ------------------------------------------------------------ the label editor

  openEditor(id, { fresh = false } = {}) {
    const it = this.item(id);
    if (!it || it.kind !== "text") return;
    this.commitEditor();
    const at = this.project(this.resolve(it.at));
    if (!at) return;
    const ta = h("textarea", { class: "ov-annot-editor", rows: String(Math.max(1, it.text.split("\n").length)), spellcheck: "false", "aria-label": "Label text" });
    ta.value = it.text;
    Object.assign(ta.style, {
      left: `${at[0] + (it.dx || 0)}px`, top: `${at[1] + (it.dy || 0)}px`,
      fontSize: `${it.size}px`, fontWeight: String(it.weight), color: it.color,
    });
    const before = this.history.length;
    this.editing = { id, el: ta, fresh, before, text: it.text };
    this.viewer.container.appendChild(ta);
    const fit = () => { ta.style.width = "0"; ta.style.width = `${Math.max(60, ta.scrollWidth + 4)}px`; ta.rows = Math.max(1, ta.value.split("\n").length); };
    fit();
    ta.addEventListener("input", () => { it.text = ta.value; fit(); this.version++; });
    ta.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); this.commitEditor(); }
      if (e.key === "Escape") { e.preventDefault(); this.commitEditor(); }
    });
    ta.addEventListener("blur", () => this.commitEditor());
    requestAnimationFrame(() => { ta.focus(); ta.select(); });
    this._sceneKey = "";
    this.viewer.renderer.invalidate();
  }

  commitEditor() {
    const ed = this.editing;
    if (!ed) return;
    this.editing = null;
    ed.el.remove();
    const it = this.item(ed.id);
    if (!it) return;
    const text = it.text.replace(/\s+$/, "");
    if (!text.trim()) {
      // An emptied label is removed (one undo brings it back).
      it.text = ed.text;
      this.removeItem(it.id);
      return;
    }
    it.text = text;
    if (text !== ed.text && !ed.fresh) {
      const cur = it.text;
      it.text = ed.text;
      this.pushHistory();
      it.text = cur;
    }
    this.changed({ history: false });
  }

  // ------------------------------------------------------------ the tool bar

  buildBar() {
    this.toolButtons = new Map();
    const tools = h("div", { class: "ov-annot-tools", role: "radiogroup", "aria-label": "Draw" });
    for (const t of TOOLS) {
      const b = h("button", { class: "ov-annot-tool", type: "button", role: "radio", "aria-label": t.label, "data-tip": `${t.label}\n${t.hint}`, onclick: () => this.setTool(t.id) }, icon(t.icon));
      this.toolButtons.set(t.id, b);
      tools.append(b);
    }
    this.undoBtn = iconButton("undo", "Undo", () => this.undo(), { shortcut: "⌘Z" });
    this.redoBtn = iconButton("redo", "Redo", () => this.redo(), { shortcut: "⇧⌘Z" });
    this.hintEl = h("span", { class: "ov-annot-hint", "aria-live": "polite" });
    this.bar = h("div", { class: "ov-annot-bar ov-glass", role: "toolbar", "aria-label": "Annotate", hidden: true },
      tools, h("span", { class: "ov-annot-sep" }), this.undoBtn, this.redoBtn, h("span", { class: "ov-annot-sep" }),
      h("button", { class: "ov-annot-done", type: "button", onclick: () => this.setActive(false) }, "Done"));
    this.ui.ui.append(this.bar, h("div", { class: "ov-annot-hint-row" }, this.hintEl));
    this.syncBar();
  }

  syncBar() {
    if (!this.bar) return;
    this.bar.hidden = !this.active;
    for (const [id, b] of this.toolButtons) b.setAttribute("aria-checked", String(id === this.tool));
    this.undoBtn.disabled = !this.history.length;
    this.redoBtn.disabled = !this.future.length;
    const t = TOOLS.find((x) => x.id === this.tool);
    this.hintEl.textContent = this.active ? t?.hint || "" : "";
    this.hintEl.parentElement.hidden = !this.active;
    this.button?.setAttribute("aria-pressed", String(this.active));
  }

  setActive(on) {
    on = !!on;
    if (on === this.active) return;
    this.active = on;
    this.ui.root.dataset.annotating = String(on);
    this.ui.setOverlay?.("annotate", on);
    if (!on) {
      this.commitEditor();
      this.drawing = null;
      this.gesture = null;
      this.sel = null;
      this.hoverId = null;
      this.snapTarget = null;
      this.viewer.canvas.style.cursor = "";
      this.ui.idle?.wake();
    }
    this.overlayKey = "";
    this._sceneKey = "";
    this.syncBar();
    this.renderPanel();
    this.viewer.renderer.invalidate();
  }

  setTool(id) {
    if (!TOOLS.some((t) => t.id === id)) return;
    if (!this.active) this.setActive(true);
    if (this.drawing) this.finishCurve();
    this.tool = id;
    if (id !== "select") this.select(null);
    this.viewer.canvas.style.cursor = id === "select" ? "" : "crosshair";
    this.syncBar();
  }

  select(id) {
    const same = (this.sel?.id || null) === (id || null);
    this.sel = id ? { id, point: same ? this.sel.point : null } : null;
    if (id && !this.active) this.setActive(true);
    if (id && this.tool !== "select") { this.tool = "select"; this.syncBar(); }
    this.overlayKey = "";
    this._sceneKey = "";
    this.renderPanel();
    this.viewer.renderer.invalidate();
  }

  // ------------------------------------------------------------ the style panel

  buildPanel() {
    this.panelBody = h("div", { class: "ov-panel-body" });
    this.panelTitle = h("span", { class: "ov-panel-title" });
    this.panel = h("section", { class: "ov-panel ov-annot-panel ov-glass ov-chrome", "aria-label": "Annotation", hidden: true },
      h("div", { class: "ov-panel-head" }, this.panelTitle, iconButton("close", "Done editing", () => this.select(null))),
      this.panelBody);
    this.ui.ui.append(this.panel);
  }

  renderPanel() {
    if (!this.panel) return;
    const it = this.active && this.sel ? this.item(this.sel.id) : null;
    this.panel.hidden = !it;
    if (!it) return;
    const body = this.panelBody;
    clear(body);
    const g = this.group(it.group);
    this.panelTitle.textContent = annotKindName(it);
    const commit = (fields) => this.patch(it.id, fields);
    const live = (fields) => this.patch(it.id, fields, { commit: false });
    const done = () => this.changed({ history: false });
    const pct = (x) => `${Math.round(x * 100)}%`;
    const fields = [];
    if (it.kind === "text") {
      const input = h("textarea", { class: "ov-input ov-annot-text", rows: "2", spellcheck: "false", "aria-label": "Text" });
      input.value = it.text;
      input.addEventListener("keydown", (e) => e.stopPropagation());
      input.addEventListener("input", () => live({ text: input.value }));
      input.addEventListener("change", () => { if (!input.value.trim()) return; this.patch(it.id, { text: input.value }); });
      fields.push(h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, "Text")), input));
      fields.push(slider({ label: "Size", min: 8, max: 64, step: 1, value: it.size, format: (x) => `${Math.round(x)} px`, onInput: (x) => live({ size: Math.round(x) }), onChange: done }));
      fields.push(miniSeg({ label: "Weight", options: [{ value: 400, label: "Regular" }, { value: 600, label: "Bold" }], value: it.weight, onChange: (x) => commit({ weight: x }) }));
      fields.push(toggle({ label: "Tag", hint: "A dark backing behind the text", checked: it.bg, onChange: (x) => commit({ bg: x }) }));
    } else if (it.kind === "curve") {
      fields.push(slider({ label: "Width", min: 0.5, max: 12, step: 0.5, value: it.width, format: (x) => `${x.toFixed(1)} px`, onInput: (x) => live({ width: x }), onChange: done }));
      fields.push(slider({ label: "Smoothness", min: 0, max: 1, step: 0.05, value: it.smooth, format: (x) => (x <= 0 ? "Straight" : pct(x)), onInput: (x) => live({ smooth: x }), onChange: done }));
      fields.push(miniSeg({ label: "Line", options: [{ value: "solid", label: "Solid" }, { value: "dash", label: "Dashed" }, { value: "dot", label: "Dotted" }], value: it.dash, onChange: (x) => commit({ dash: x }) }));
      fields.push(miniSeg({ label: "Arrowheads", options: [{ value: "none", label: "None" }, { value: "end", label: "End" }, { value: "start", label: "Start" }, { value: "both", label: "Both" }], value: it.arrow, onChange: (x) => commit({ arrow: x }) }));
      if (it.points.length > 2) fields.push(toggle({ label: "Closed loop", checked: it.closed, onChange: (x) => commit({ closed: x }) }));
    } else {
      fields.push(miniSeg({ label: "Look", options: [{ value: "bubble", label: "Bubble" }, { value: "shell", label: "Shell" }, { value: "wire", label: "Wire" }], value: it.style, onChange: (x) => commit({ style: x }) }));
      const r = Math.max(...it.radii);
      fields.push(slider({ label: "Radius", min: Math.max(0.1, r / 20), max: r * 4, value: r, scale: "log", format: (x) => fmtPc(x), onInput: (x) => live({ radii: it.radii.map((k) => (k / r) * x) }), onChange: done }));
      const axes = h("div", { class: "ov-annot-xyz" }, ...["x", "y", "z"].map((a, k) => numberInput({ label: `R${a} (pc)`, value: round(it.radii[k]), onChange: (x) => { if (x > 0) commit({ radii: it.radii.map((y0, j) => (j === k ? x : y0)) }); } })));
      fields.push(axes);
      const rot = h("div", { class: "ov-annot-xyz" }, ...["x", "y", "z"].map((a, k) => numberInput({ label: `Tilt ${a} (°)`, value: round(it.rot[k]), onChange: (x) => commit({ rot: it.rot.map((y0, j) => (j === k ? x : y0)) }) })));
      fields.push(rot);
    }
    fields.push(h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, "Colour")),
      colorPicker({ value: it.color, onChange: (c) => commit({ color: c }) })));
    fields.push(slider({ label: "Opacity", min: 0.05, max: 1, step: 0.01, value: it.opacity, format: pct, onInput: (x) => live({ opacity: x }), onChange: done }));
    // Where it is: the point being edited, the centre, or the label's anchor.
    const pt = it.kind === "text" ? it.at : it.kind === "sphere" ? it.center : it.points[this.sel.point ?? -1] || null;
    if (pt) {
      const setPt = (p) => {
        if (it.kind === "text") commit({ at: p });
        else if (it.kind === "sphere") commit({ center: p });
        else commit({ points: it.points.map((q, k) => (k === this.sel.point ? p : q)) });
        this.renderPanel();
      };
      if (pt.trace != null) {
        const name = this.viewer.describeObject(pt.trace, pt.index)?.name || "a cluster";
        fields.push(h("div", { class: "ov-annot-attached" }, icon("follow"), h("span", { class: "ov-grow" }, `Follows ${name}`),
          h("button", { class: "ov-btn ov-btn--ghost", type: "button", onclick: () => { const p = this.resolve(pt); if (p) setPt({ pos: p.map(round) }); } }, "Fix in place")));
      } else {
        const label = it.kind === "curve" ? `Point ${this.sel.point + 1}` : "Position";
        fields.push(h("div", { class: "ov-field-head ov-annot-sub" }, h("span", { class: "ov-field-label" }, `${label} · pc`)));
        fields.push(h("div", { class: "ov-annot-xyz" }, ...["x", "y", "z"].map((a, k) => numberInput({ label: a, value: round(pt.pos[k]), onChange: (x) => setPt({ pos: pt.pos.map((y0, j) => (j === k ? x : y0)) }) }))));
      }
    } else if (it.kind === "curve") {
      fields.push(h("p", { class: "ov-field-hint" }, "Drag a point to move it (Shift: along the line of sight, ⌘: onto a cluster or the dust). Drag a “+” to bend; double-click a point to remove it."));
    }
    fields.push(toggle({ label: "Present day only", hint: "Fades out away from t = 0, like the dust maps", checked: !!it.present, onChange: (x) => commit({ present: x || undefined }) }));
    // Group: the key entry it belongs to.
    const groups = this.doc.groups.map((x) => ({ value: x.id, label: annotGroupName(this.doc, x) }));
    fields.push(select({ label: "Key entry", options: [...groups, { value: "__new", label: "New entry…" }], value: it.group, onChange: (val) => this.moveToGroup(it.id, val) }));
    if (g) {
      const name = h("input", { class: "ov-input", type: "text", value: annotGroupName(this.doc, g), "aria-label": "Key entry name" });
      name.addEventListener("keydown", (e) => { e.stopPropagation(); if (e.key === "Enter") name.blur(); });
      name.addEventListener("change", () => this.renameGroup(g.id, name.value));
      fields.push(h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, "Entry name")), name));
    }
    fields.push(h("div", { class: "ov-annot-actions" },
      h("button", { class: "ov-btn ov-btn--ghost", type: "button", onclick: () => this.duplicate(it.id) }, icon("copy"), "Duplicate"),
      h("button", { class: "ov-btn ov-btn--ghost ov-annot-delete", type: "button", onclick: () => this.removeItem(it.id) }, icon("trash"), "Delete")));
    body.append(...fields);
  }

  moveToGroup(id, groupId) {
    const it = this.item(id);
    if (!it) return;
    this.pushHistory();
    if (groupId === "__new") {
      const g = { id: annotId("g"), name: "", auto: true, visible: true };
      this.doc.groups.push(g);
      groupId = g.id;
    }
    it.group = groupId;
    const used = new Set(this.doc.items.map((i) => i.group));
    this.doc.groups = this.doc.groups.filter((g) => used.has(g.id));
    this.changed({ history: false });
  }

  renameGroup(id, name) {
    const g = this.group(id);
    if (!g) return;
    name = String(name || "").trim();
    this.pushHistory();
    if (name) { g.name = name; g.auto = false; delete g.nameFrom; } else { g.name = ""; g.auto = true; }
    this.changed({ history: false });
  }

  // ------------------------------------------------------------ the key and the layers panel

  /** Annotation groups as key entries (after the data layers). */
  toggleItems() {
    return this.doc.groups.map((g) => ({ id: `annotation:${g.id}`, kind: "annotation", key: g.id, name: annotGroupName(this.doc, g) }));
  }

  groupVisible(id) {
    return this.group(id)?.visible !== false;
  }

  setGroupVisible(id, on) {
    const g = this.group(id);
    if (!g || (g.visible !== false) === !!on) return;
    g.visible = !!on;
    this.version++;
    this._sceneKey = "";
    this.viewer.renderer.invalidate();
    this.ui.legend?.sync();
    this.ui.layers?.sync?.();
    this.saveDraft();
  }

  /** The key's dot for a group: its first shape's mark in its colour. */
  legendDot(dot, groupId) {
    const g = this.group(groupId);
    if (!g) return;
    const look = annotGroupLook(this.doc, g);
    dot.dataset.kind = "annotation";
    dot.dataset.mark = look.kind;
    dot.style.background = "";
    dot.style.color = look.color;
    if (dot.dataset.markDrawn !== look.kind) {
      dot.innerHTML = MARKS[look.kind] || MARKS.text;
      dot.dataset.markDrawn = look.kind;
    }
  }

  /** The Annotations section of the layers panel. */
  layersSection() {
    const doc = this.doc;
    const sec = h("div", { class: "ov-section ov-annot-section" });
    const head = h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Annotations"),
      iconButton("pen", this.active ? "Stop annotating" : "Annotate", () => this.setActive(!this.active), { shortcut: "K", cls: "ov-section-eye" }));
    sec.append(head);
    if (!doc.groups.length) {
      sec.append(h("p", { class: "ov-annot-empty" }, "Draw labels, curves, arrows, bubbles and shells in 3D. Press ", h("kbd", null, "K"), " to start."));
      return sec;
    }
    for (const g of doc.groups) {
      const items = doc.items.filter((i) => i.group === g.id);
      const dot = h("span", { class: "ov-legend-dot" });
      this.legendDot(dot, g.id);
      const eye = h("button", { class: "ov-swatch-btn", type: "button", "aria-label": `Show or hide ${annotGroupName(doc, g)}` }, dot);
      eye.addEventListener("click", () => { this.setGroupVisible(g.id, !this.groupVisible(g.id)); this.ui.layers?.render(); });
      const name = h("button", { class: "ov-row-name", type: "button", title: "Click to show or hide · double-click to rename" }, annotGroupName(doc, g));
      name.addEventListener("click", () => { this.setGroupVisible(g.id, !this.groupVisible(g.id)); this.ui.layers?.render(); });
      name.addEventListener("dblclick", () => {
        const input = h("input", { class: "ov-input", type: "text", value: annotGroupName(doc, g) });
        input.addEventListener("keydown", (e) => { e.stopPropagation(); if (e.key === "Enter") input.blur(); if (e.key === "Escape") { input.value = annotGroupName(doc, g); input.blur(); } });
        input.addEventListener("blur", () => this.renameGroup(g.id, input.value));
        name.replaceWith(input);
        input.focus();
        input.select();
      });
      const row = h("div", { class: "ov-row ov-annot-row", "data-visible": String(g.visible !== false) }, eye, name, h("span", { class: "ov-annot-count" }, String(items.length)));
      sec.append(row);
      for (const it of items) {
        const b = h("button", { class: "ov-annot-item", type: "button", "aria-pressed": String(this.sel?.id === it.id), onclick: () => this.select(it.id) },
          h("span", { class: "ov-annot-mark", style: { color: it.color }, html: MARKS[annotKind(it)] || MARKS.text }),
          h("span", { class: "ov-grow" }, it.kind === "text" ? it.text.split("\n")[0] : annotKindName(it)));
        sec.append(b);
      }
    }
    return sec;
  }

  // ------------------------------------------------------------ keys and commands

  /** Keys while annotating (called before the viewer's own shortcuts). */
  handleKey(e) {
    if (!this.active) return false;
    const mod = e.metaKey || e.ctrlKey;
    const k = e.key;
    if (mod && k.toLowerCase() === "z") { if (e.shiftKey) this.redo(); else this.undo(); return true; }
    if (mod && k.toLowerCase() === "y") { this.redo(); return true; }
    if (mod && k.toLowerCase() === "d" && this.sel) { this.duplicate(this.sel.id); return true; }
    if (mod) return false;
    if (k === "Escape") {
      if (this.drawing) this.finishCurve();
      else if (this.sel) this.select(null);
      else if (this.tool !== "select") this.setTool("select");
      else this.setActive(false);
      return true;
    }
    if (k === "Enter" && this.drawing) { this.finishCurve(); return true; }
    if (k === "Enter" && this.sel && this.item(this.sel.id)?.kind === "text") { this.openEditor(this.sel.id); return true; }
    if ((k === "Delete" || k === "Backspace") && this.sel) {
      const it = this.item(this.sel.id);
      if (it?.kind === "curve" && this.sel.point != null && it.points.length > 2) {
        this.pushHistory();
        it.points.splice(this.sel.point, 1);
        if (it.points.length < 3) it.closed = false;
        this.sel.point = null;
        this.changed({ history: false });
      } else this.removeItem(this.sel.id);
      return true;
    }
    return false;
  }

  /** Global keys: K starts and stops annotating. */
  onKey(e) {
    if (e.metaKey || e.ctrlKey || e.altKey) return false;
    if (e.key === "k" || e.key === "K") { this.setActive(!this.active); return true; }
    return false;
  }

  helpGroups() {
    const mod = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.platform) ? "⌘" : "Ctrl";
    return [["Annotate", [
      ["Start / stop annotating", "K"], ["Finish a curve", "Enter"], ["Move along the line of sight", "⇧ Drag"],
      ["Snap a point onto a cluster or the dust", `${mod} Drag`], ["Undo / redo", `${mod} Z / ⇧${mod} Z`],
      ["Duplicate", `${mod} D`], ["Delete", "⌫"], ["Deselect · done", "Esc"],
    ]]];
  }

  commands() {
    const n = this.doc.items.length;
    const hidden = this.doc.groups.filter((g) => g.visible === false).length;
    return [
      { title: this.active ? "Stop annotating" : "Annotate: draw in 3D", icon: "pen", shortcut: "K", keywords: "draw sketch annotate label note text arrow curve bubble shell sphere", run: () => this.setActive(!this.active) },
      ...TOOLS.filter((t) => t.id !== "select").map((t) => ({ title: `Draw: ${t.label}`, icon: t.icon, keywords: `annotate ${t.label.toLowerCase()} add`, run: () => this.setTool(t.id) })),
      ...(hidden ? [{ title: "Show every annotation", icon: "eye", run: () => { for (const g of this.doc.groups) g.visible = true; this.changed(); } }] : []),
      ...(n ? [{ title: "Remove all annotations", sub: `${n} item${n === 1 ? "" : "s"} · ⌘Z brings them back`, icon: "trash", run: () => this.setDoc(emptyAnnotations()) }] : []),
    ];
  }

  // ------------------------------------------------------------ State transitions

  transition(a, b) {
    const A = a && typeof a === "object" ? normalizeAnnotations(a) : this.doc;
    const B = b && typeof b === "object" ? normalizeAnnotations(b) : A;
    if (JSON.stringify(A) === JSON.stringify(B)) return null;
    const resolve = (p) => this.resolve(p);
    return {
      step: (t) => {
        this.view = lerpAnnotations(A, B, t, resolve);
        this.version++;
        this._sceneKey = "";
      },
    };
  }

  // ------------------------------------------------------------ drafts and export

  draftKey() {
    return `oviz.annotations.${figureId(this.viewer.manifest)}`;
  }

  keepsDrafts() {
    const story = this.ui.plugins.find((p) => p.name === "states");
    return !(story?.readOnly || story?.project?.autosave === false);
  }

  restoreDraft() {
    if (!this.keepsDrafts()) return;
    const raw = localStorageGet(this.draftKey());
    if (!raw || raw === this.base || raw === JSON.stringify(this.doc)) return;
    let doc;
    try { doc = normalizeAnnotations(JSON.parse(raw)); } catch (_) { return; }
    const backup = this.doc;
    this.setDoc(doc, { history: false });
    this.ui.toast(`Your annotations are back (${doc.items.length})`, {
      icon: icon("pen"), ms: 6000,
      action: { label: "Discard", run: () => { this.setDoc(backup, { history: false }); localStorageSet(this.draftKey(), ""); } },
    });
  }

  saveDraft() {
    if (!this.keepsDrafts()) return;
    clearTimeout(this._draftTimer);
    this._draftTimer = setTimeout(() => {
      const s = JSON.stringify(this.doc);
      localStorageSet(this.draftKey(), s === this.base ? "" : s);
    }, 400);
  }

  /** Saved figures keep what was drawn (it is part of the figure). */
  exportManifest(m) {
    if (this.doc.items.length) m.annotations = cloneJson(this.doc);
    else delete m.annotations;
  }
}

// ------------------------------------------------------------------ helpers

/** Glyphs for the key: a label, a curve, an arrow, a filled bubble, a ring. */
const MARKS = {
  text: '<svg viewBox="0 0 12 12" width="12" height="12"><path d="M2.5 2.5h7M6 2.5v7.5" stroke="currentColor" stroke-width="1.8" fill="none" stroke-linecap="round"/></svg>',
  curve: '<svg viewBox="0 0 12 12" width="12" height="12"><path d="M1.5 9.5C3 3 8 10 10.5 2.5" stroke="currentColor" stroke-width="1.6" fill="none" stroke-linecap="round"/></svg>',
  arrow: '<svg viewBox="0 0 12 12" width="12" height="12"><path d="M2 10 9.5 2.5M5 2.5h4.5V7" stroke="currentColor" stroke-width="1.6" fill="none" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  bubble: '<svg viewBox="0 0 12 12" width="12" height="12"><circle cx="6" cy="6" r="4.6" fill="currentColor" fill-opacity="0.45" stroke="currentColor" stroke-width="1.2"/></svg>',
  shell: '<svg viewBox="0 0 12 12" width="12" height="12"><circle cx="6" cy="6" r="4.4" fill="none" stroke="currentColor" stroke-width="1.8"/></svg>',
  wire: '<svg viewBox="0 0 12 12" width="12" height="12"><circle cx="6" cy="6" r="4.5" fill="none" stroke="currentColor" stroke-width="1.1"/><ellipse cx="6" cy="6" rx="4.5" ry="1.7" fill="none" stroke="currentColor" stroke-width="0.9"/><ellipse cx="6" cy="6" rx="1.7" ry="4.5" fill="none" stroke="currentColor" stroke-width="0.9"/></svg>',
};

function mergeInto(doc, extra) {
  const ids = new Set(doc.items.map((i) => i.id));
  const items = extra.items.filter((i) => !ids.has(i.id));
  if (!items.length) return doc;
  const groups = new Set(items.map((i) => i.group));
  for (const g of extra.groups) if (groups.has(g.id) && !doc.groups.some((x) => x.id === g.id)) doc.groups.push(g);
  doc.items.push(...items);
  return doc;
}

/** A classic manual label as a note. */
function legacyLabel(l, i) {
  if (!l || typeof l !== "object") return null;
  const x = Number(l.x ?? l.position?.x), y = Number(l.y ?? l.position?.y), z = Number(l.z ?? l.position?.z);
  if (![x, y, z].every(Number.isFinite)) return null;
  const size = Number(l.size ?? l.font_size);
  return { id: String(l.id || `legacy-${i}`), text: String(l.text || "Note"), anchor: { pos: [x, y, z] }, color: String(l.color || "#ffffff"), ...(size > 0 ? { size } : {}), ...(l.visible === false ? { visible: false } : {}) };
}

function add3(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function sub3(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function dot3(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function dist3(a, b) { return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }
function round(x) { return Math.round(x * 100) / 100; }
function fmtPc(x) { return x >= 1000 ? `${(x / 1000).toFixed(2)} kpc` : `${x >= 100 ? Math.round(x) : x.toFixed(1)} pc`; }
