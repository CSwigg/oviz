// Birth tree (classic Oviz "Dendrogram"): each cluster is linked to its
// nearest older cluster at its birth, either to where the older cluster was
// at that moment ("birth to older track") or to where it was born ("birth
// to birth"), within a distance threshold or a birth-age threshold. Hovering
// a branch highlights its descendants in 3D (the rest dims); clicking pins
// it. Built from the trace data already in the figure.

import { h, iconButton, clear } from "./dom.js";
import { select, numberInput } from "./controls.js";
import { clamp } from "../core/math.js";
import { rafThrottle } from "../core/emitter.js";

export const TREE_BIT = 32; // per-object state bit: dimmed by the birth tree

/**
 * Link clusters into a birth tree (the classic `buildDendrogramModel`).
 * `entries`: `{ key, name, age, birthTime, birth: [x, y, z], at(t) }` where
 * `at(t)` is the cluster's position at time t (Myr) or null. Returns
 * `{ nodes, roots, branches, leafCount, maxAge, maxAxis, mode, links }`;
 * nodes carry `parent`, `children`, `axis` (birth age, or cumulative birth
 * distance in "birth_age_myr" mode) and a leaf `order` for plotting.
 */
export function buildBirthTree(entries, { threshold = 100, mode = "distance_pc", links = "birth_to_older_track" } = {}) {
  const byAge = mode === "birth_age_myr";
  const toBirth = links === "birth_to_birth";
  const limit = Math.max(Number(threshold) || 0, 0);
  const nodes = (entries || [])
    .filter((e) => e && e.key && Number.isFinite(e.age) && Number.isFinite(e.birthTime) && e.birth && e.birth.every(Number.isFinite))
    .map((e) => ({ ...e, parent: null, dist: NaN, gap: NaN, axis: 0, children: [], order: 0 }))
    .sort((a, b) => b.age - a.age || String(a.name).localeCompare(String(b.name)));
  const byKey = new Map(nodes.map((n) => [n.key, n]));
  for (let i = 0; i < nodes.length; i++) {
    const node = nodes[i];
    let best = null, bestD = Infinity, bestGap = Infinity;
    for (let j = 0; j < i; j++) {
      const cand = nodes[j];
      const gap = Math.max(0, cand.age - node.age);
      if (byAge && !(gap <= limit)) continue;
      const p = toBirth ? cand.birth : cand.at?.(node.birthTime);
      if (!p) continue;
      const d = Math.hypot(node.birth[0] - p[0], node.birth[1] - p[1], node.birth[2] - p[2]);
      if (!byAge && !(d <= limit)) continue;
      if (d >= bestD) continue;
      best = cand; bestD = d; bestGap = gap;
    }
    if (best) {
      node.parent = best.key;
      node.dist = bestD;
      node.gap = bestGap;
      best.children.push(node.key);
    }
    node.axis = byAge ? 0 : Math.max(node.age, 0);
  }
  if (byAge) {
    // The height is the distance travelled from the root, birth to birth.
    const memo = new Map();
    const cumulative = (n) => {
      if (!n) return 0;
      if (memo.has(n.key)) return memo.get(n.key);
      memo.set(n.key, 0); // guards against a cycle
      const v = (n.parent ? cumulative(byKey.get(n.parent)) : 0) + (Number.isFinite(n.dist) ? Math.max(n.dist, 0) : 0);
      memo.set(n.key, v);
      return v;
    };
    for (const n of nodes) n.axis = cumulative(n);
  }
  const order = (a, b) => b.age - a.age || String(a.name).localeCompare(String(b.name));
  const roots = nodes.filter((n) => !n.parent).sort(order);
  let leaf = 0;
  const place = (n) => {
    const kids = n.children.map((k) => byKey.get(k)).filter(Boolean).sort(order);
    if (!kids.length) { n.order = leaf++; return n.order; }
    const o = kids.map(place);
    n.order = o.reduce((s, x) => s + x, 0) / o.length;
    return n.order;
  };
  roots.forEach(place);
  const desc = new Map();
  const descendants = (n) => {
    if (desc.has(n.key)) return desc.get(n.key);
    const keys = [n.key];
    desc.set(n.key, keys);
    for (const k of n.children) { const c = byKey.get(k); if (c) keys.push(...descendants(c)); }
    return keys;
  };
  const branches = [];
  for (const n of nodes) {
    const parent = n.parent ? byKey.get(n.parent) : null;
    if (!parent) continue;
    const keys = descendants(n);
    branches.push({ key: `${n.key}->${parent.key}`, child: n, parent, label: `${n.name} from ${parent.name}`, keys, count: keys.length, distance: n.dist, gap: n.gap });
  }
  const axes = nodes.map((n) => n.axis).filter((x) => Number.isFinite(x) && x >= 0);
  return {
    nodes, roots, branches, byKey,
    leafCount: Math.max(leaf, 1),
    maxAge: nodes.reduce((m, n) => Math.max(m, n.age), 0),
    maxAxis: axes.length ? Math.max(...axes) : 0,
    mode: byAge ? "birth_age_myr" : "distance_pc",
    links: toBirth ? "birth_to_birth" : "birth_to_older_track",
  };
}

export class BirthTreePlugin {
  constructor() {
    this.name = "birth-tree";
    this.opts = null;
    this.pinned = null; // hit region key
    this.hovered = null;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    const cfg = v.manifest.widgets?.birthTree || {};
    this.title = cfg.title || "Birth tree";
    this.opts = {
      trace: null,
      mode: cfg.mode === "birth_age_myr" ? "birth_age_myr" : "distance_pc",
      links: cfg.links === "birth_to_birth" ? "birth_to_birth" : "birth_to_older_track",
      thresholdPc: Number.isFinite(cfg.thresholdPc) ? cfg.thresholdPc : 100,
      thresholdAgeMyr: Number.isFinite(cfg.thresholdAgeMyr) ? cfg.thresholdAgeMyr : 5,
    };
    this.defaultTrace = cfg.defaultTrace;
    ui.widgets.register({
      id: "birth-tree", title: this.title, icon: "trail",
      description: "Which clusters were born near older ones: a dendrogram of birth sites",
      size: { w: 560, h: 380 }, minSize: { w: 380, h: 260 },
      available: () => this.traceOptions().length > 0,
      mount: (panel) => this.mount(panel),
    });
    v.on("gpu-restored", () => this.applyHighlight());
  }

  /** Point traces in the current group whose objects have ages (and more than one object). */
  traceOptions() {
    const v = this.viewer;
    return v.traces.filter((t) => t.points && t.showInLegend && t.points.count > 1 && t.points.ageNow?.blob && v.state.traces[t.key]?.inGroup !== false);
  }

  currentTrace() {
    const opts = this.traceOptions();
    if (!opts.length) return null;
    const want = this.opts.trace || this.defaultTrace;
    const hit = opts.find((t) => t.key === want);
    if (hit) return hit;
    // Otherwise the largest visible trace.
    const vis = opts.filter((t) => this.viewer.state.traces[t.key]?.visible);
    return (vis.length ? vis : opts).slice().sort((a, b) => b.points.count - a.points.count)[0];
  }

  // ---------------------------------------------------------------- data

  entriesFor(trace) {
    const v = this.viewer;
    const tl = v.timeline;
    const d = v.data.get(trace.key);
    const meta = v.objectMetaSync?.(trace.key) || {};
    if (!d?.ageNow) return [];
    const lo = Math.min(tl.times[0], tl.times[tl.count - 1]), hi = Math.max(tl.times[0], tl.times[tl.count - 1]);
    const out = [];
    for (let i = 0; i < trace.points.count; i++) {
      const age = d.ageNow[i];
      if (!(age === age) || age < 0) continue;
      const birthTime = -age;
      // Born before the timeline starts: use its earliest position (classic).
      const bt = clamp(birthTime, lo, hi);
      const birth = v.objectPosition(trace.key, i, tl.timeToFrame(bt));
      if (!birth) continue;
      out.push({
        key: `${trace.key}:${i}`, index: i,
        name: String(meta.name?.[i] || `#${i + 1}`).replace(/_/g, " "),
        age, birthTime, birth: birth.slice(),
        at: (t) => v.objectPosition(trace.key, i, tl.timeToFrame(clamp(t, lo, hi))),
      });
    }
    return out;
  }

  model() {
    const trace = this.currentTrace();
    if (!trace) return null;
    const o = this.opts;
    const threshold = o.mode === "birth_age_myr" ? o.thresholdAgeMyr : o.thresholdPc;
    const sig = `${trace.key}|${o.mode}|${o.links}|${threshold}`;
    if (this._model?.sig === sig) return this._model;
    const t0 = performance.now();
    const m = buildBirthTree(this.entriesFor(trace), { threshold, mode: o.mode, links: o.links });
    m.sig = sig;
    m.trace = trace;
    m.ms = performance.now() - t0;
    this._model = m;
    return m;
  }

  // ---------------------------------------------------------------- widget

  mount(panel) {
    const v = this.viewer;
    this.panel = panel;
    const body = panel.body;
    body.classList.add("ov-tree");
    this.controls = h("div", { class: "ov-tree-controls" });
    this.canvas = h("canvas", { class: "ov-tree-canvas", "aria-label": "Birth tree" });
    this.tip = h("div", { class: "ov-tree-tip", hidden: true });
    this.foot = h("div", { class: "ov-tree-foot" });
    body.append(this.controls, h("div", { class: "ov-tree-plot" }, this.canvas, this.tip), this.foot);
    this.renderControls();
    const redraw = rafThrottle(() => this.draw());
    this.redraw = redraw;
    const offs = [
      v.on("time", redraw),
      v.on("style", () => { this.renderControls(); redraw(); }),
      v.on("theme", redraw),
    ];
    this.canvas.addEventListener("pointermove", (e) => this.hover(e));
    this.canvas.addEventListener("pointerleave", () => { this.hovered = null; this.tip.hidden = true; this.applyHighlight(); redraw(); });
    this.canvas.addEventListener("click", (e) => this.click(e));
    // Names come with the object metadata; draw again once it is in.
    const withMeta = () => {
      const t = this.currentTrace();
      if (t && !v.objectMetaSync?.(t.key)?.name) v.objectMeta(t.key).then(() => { this._model = null; redraw(); });
    };
    return {
      onOpen: () => { withMeta(); this.renderControls(); redraw(); this.applyHighlight(); },
      onClose: () => { this.hovered = null; this.clearHighlight(); },
      onResize: () => redraw(),
      capture: () => ({ ...this.opts, trace: this.currentTrace()?.key || null, pinned: this.pinned }),
      apply: (s) => {
        if (!s || typeof s !== "object") return;
        for (const k of ["trace", "mode", "links", "thresholdPc", "thresholdAgeMyr"]) if (s[k] != null) this.opts[k] = s[k];
        this.pinned = s.pinned || null;
        this.renderControls();
        redraw();
        this.applyHighlight();
      },
      destroy: () => offs.forEach((off) => off?.()),
    };
  }

  renderControls() {
    const c = this.controls;
    if (!c) return;
    clear(c);
    const traces = this.traceOptions();
    const cur = this.currentTrace();
    const o = this.opts;
    const age = o.mode === "birth_age_myr";
    const set = (patch) => {
      Object.assign(this.opts, patch);
      this.pinned = null;
      this.renderControls();
      const t = this.currentTrace();
      if (t && !this.viewer.objectMetaSync?.(t.key)?.name) this.viewer.objectMeta(t.key).then(() => { this._model = null; this.redraw?.(); });
      this.redraw?.();
      this.applyHighlight();
    };
    c.append(
      select({ label: "Trace", value: cur?.key, options: traces.map((t) => ({ value: t.key, label: `${t.name} (${t.points.count})` })), onChange: (k) => set({ trace: k }) }),
      select({ label: "Links", value: o.links, options: [{ value: "birth_to_older_track", label: "Birth to older track" }, { value: "birth_to_birth", label: "Birth to birth" }], onChange: (x) => set({ links: x }) }),
      select({ label: "Mode", value: o.mode, options: [{ value: "distance_pc", label: "Distance threshold" }, { value: "birth_age_myr", label: "Birth-age threshold" }], onChange: (x) => set({ mode: x }) }),
      numberInput({ label: age ? "Threshold (Myr)" : "Threshold (pc)", value: age ? o.thresholdAgeMyr : o.thresholdPc, step: age ? "0.5" : "5", onChange: (x) => set(age ? { thresholdAgeMyr: Math.max(0, x) } : { thresholdPc: Math.max(0, x) }) }),
    );
  }

  /** Canvas layout of the current model (used for drawing and hit tests). */
  layout() {
    const m = this.model();
    const cv = this.canvas;
    const W = cv.clientWidth || 1, H = cv.clientHeight || 1;
    const margin = { left: 56, right: 14, top: 14, bottom: 22 };
    const pw = Math.max(40, W - margin.left - margin.right), ph = Math.max(40, H - margin.top - margin.bottom);
    if (!m) return { m, W, H, margin, pw, ph };
    const age = m.mode === "birth_age_myr";
    const axisMax = age ? Math.max(m.maxAxis * 1.08, 25) : Math.max(m.maxAge + 5, 1);
    const x = (order) => (m.leafCount <= 1 ? margin.left + pw / 2 : margin.left + (order / Math.max(m.leafCount - 1, 1)) * pw);
    const y = (val) => margin.top + (1 - clamp(val / axisMax, 0, 1)) * ph;
    return { m, W, H, margin, pw, ph, axisMax, x, y, age };
  }

  draw() {
    const cv = this.canvas;
    if (!cv || !cv.isConnected) return;
    const L = this.layout();
    const { m, W, H, margin, pw, ph } = L;
    const dpr = Math.max(window.devicePixelRatio || 1, 1);
    if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const css = getComputedStyle(this.ui.root);
    const axisColor = css.getPropertyValue("--ov-text-3").trim() || "rgba(236,240,246,0.5)";
    const textColor = css.getPropertyValue("--ov-text-2").trim() || "rgba(236,240,246,0.74)";
    if (!m || !m.nodes.length) {
      ctx.fillStyle = textColor;
      ctx.font = "12px system-ui, sans-serif";
      ctx.fillText("No clusters with ages in this trace.", margin.left, margin.top + 20);
      this.regions = [];
      this.foot.textContent = "";
      return;
    }
    const v = this.viewer;
    const ts = v.state.traces[m.trace.key] || {};
    const color = ts.color || m.trace.legendColor || m.trace.color || "#ffffff";
    // Axes and ticks.
    ctx.strokeStyle = axisColor;
    ctx.fillStyle = axisColor;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top); ctx.lineTo(margin.left, margin.top + ph); ctx.lineTo(margin.left + pw, margin.top + ph);
    ctx.stroke();
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (const f of [0, 0.33, 0.66, 1]) {
      const val = f * L.axisMax, py = L.y(val);
      ctx.globalAlpha = 0.14;
      ctx.beginPath(); ctx.moveTo(margin.left, py); ctx.lineTo(margin.left + pw, py); ctx.stroke();
      ctx.globalAlpha = 1;
      ctx.fillText(`${compact(val)} ${L.age ? "pc" : "Myr"}`, margin.left - 6, py);
    }
    ctx.save();
    ctx.translate(12, margin.top + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText(L.age ? "Birth distance" : "Birth age", 0, 0);
    ctx.restore();
    // Formation age now: a dashed line at −t (distance mode).
    const formation = Math.max(0, -v.timeline.time);
    if (!L.age) {
      const py = L.y(Math.min(formation, L.axisMax));
      ctx.save(); ctx.setLineDash([5, 5]); ctx.globalAlpha = 0.8;
      ctx.beginPath(); ctx.moveTo(margin.left, py); ctx.lineTo(margin.left + pw, py); ctx.stroke();
      ctx.restore();
    }
    const active = this.activeKeys();
    const regions = [];
    for (const b of m.branches) {
      const cx = L.x(b.child.order), cy = L.y(b.child.axis), px = L.x(b.parent.order), py = L.y(b.parent.axis);
      const pinned = this.pinned === b.key, hov = !pinned && this.hovered === b.key;
      const on = active && b.keys.some((k) => active.has(k));
      ctx.strokeStyle = color;
      ctx.globalAlpha = pinned ? 1 : hov ? 0.96 : on ? 0.86 : active ? 0.25 : 0.62;
      ctx.lineWidth = pinned ? 3.4 : hov ? 2.8 : on ? 2.2 : 1.5;
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx, py); ctx.lineTo(px, py); ctx.stroke();
      regions.push({ type: "branch", key: b.key, label: `${b.label} · ${b.count} cluster${b.count === 1 ? "" : "s"} · ${Number.isFinite(b.distance) ? `${b.distance.toFixed(0)} pc` : ""}${b.gap > 0 ? ` · ${b.gap.toFixed(1)} Myr older` : ""}`, keys: b.keys, segs: [[cx, cy, cx, py], [cx, py, px, py]] });
    }
    ctx.globalAlpha = 1;
    for (const n of m.nodes) {
      const nx = L.x(n.order), ny = L.y(n.axis);
      const rk = `node:${n.key}`;
      const pinned = this.pinned === rk, hov = !pinned && this.hovered === rk;
      const on = active && active.has(n.key);
      const r = pinned ? 5 : hov ? 4.6 : on ? 3.8 : 3.1;
      ctx.fillStyle = color;
      ctx.globalAlpha = pinned || hov ? 1 : on ? 0.92 : active ? 0.3 : 0.82;
      ctx.beginPath(); ctx.arc(nx, ny, r, 0, Math.PI * 2); ctx.fill();
      regions.push({ type: "node", key: rk, label: `${n.name} · ${n.age.toFixed(1)} Myr`, keys: [n.key], node: n, cx: nx, cy: ny, r: Math.max(7, r + 2.5) });
    }
    ctx.globalAlpha = 1;
    this.regions = regions;
    const linked = m.branches.length;
    this.foot.textContent = `${m.nodes.length.toLocaleString()} clusters · ${linked.toLocaleString()} links · ${m.roots.length.toLocaleString()} roots · hover a branch to highlight it in 3D, click to pin`;
  }

  hitAt(e) {
    const r = this.canvas.getBoundingClientRect();
    const x = e.clientX - r.left, y = e.clientY - r.top;
    let best = null, bestD = Infinity;
    for (const g of this.regions || []) {
      if (g.type === "node") {
        const d = Math.hypot(x - g.cx, y - g.cy);
        if (d <= g.r && d < bestD) { best = g; bestD = d - 1; }
      } else {
        for (const [x0, y0, x1, y1] of g.segs) {
          const dx = x1 - x0, dy = y1 - y0, l2 = dx * dx + dy * dy;
          const u = l2 > 0 ? clamp(((x - x0) * dx + (y - y0) * dy) / l2, 0, 1) : 0;
          const d = Math.hypot(x0 + dx * u - x, y0 + dy * u - y);
          if (d <= 7 && d < bestD) { best = g; bestD = d; }
        }
      }
    }
    return best ? { ...best, x, y } : null;
  }

  hover(e) {
    const hit = this.hitAt(e);
    const key = hit?.key || null;
    if (hit) {
      this.tip.textContent = hit.label;
      this.tip.hidden = false;
      this.tip.style.left = `${clamp(hit.x + 12, 4, (this.canvas.clientWidth || 200) - 200)}px`;
      this.tip.style.top = `${Math.max(4, hit.y - 28)}px`;
    } else this.tip.hidden = true;
    if (key !== this.hovered) {
      this.hovered = key;
      this.applyHighlight();
      this.redraw();
    }
  }

  click(e) {
    const hit = this.hitAt(e);
    this.pinned = hit && hit.key !== this.pinned ? hit.key : null;
    // A node click also selects that cluster in the figure.
    if (hit?.type === "node") this.ui.select({ trace: this.model().trace.key, index: hit.node.index });
    this.applyHighlight();
    this.redraw();
  }

  /** Keys of the highlighted clusters (the pinned or hovered branch/node), or null. */
  activeKeys() {
    const key = this.hovered || this.pinned;
    if (!key) return null;
    const g = (this.regions || []).find((r) => r.key === key);
    return g ? new Set(g.keys) : null;
  }

  applyHighlight() {
    const v = this.viewer;
    const open = this.ui.widgets.isOpen("birth-tree");
    const active = open ? this.activeKeys() : null;
    const m = active ? this.model() : null;
    for (const batch of v.points.batches) {
      if (!active || !m) { batch.setStateBits(TREE_BIT, null); continue; }
      // Everything but the highlighted clusters dims, other traces included.
      const bits = new Uint8Array(batch.count).fill(TREE_BIT);
      if (batch.key === m.trace.key) for (const k of active) bits[Number(k.slice(k.lastIndexOf(":") + 1))] = 0;
      batch.setStateBits(TREE_BIT, bits);
    }
    v.invalidatePick();
    v.renderer.invalidate();
  }

  clearHighlight() {
    for (const batch of this.viewer.points.batches) batch.setStateBits(TREE_BIT, null);
    this.viewer.renderer.invalidate();
  }
}

function compact(v) {
  const a = Math.abs(v);
  if (a >= 1000) return `${(v / 1000).toFixed(a >= 10000 ? 0 : 1)}k`;
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(0);
  return v.toFixed(1);
}
