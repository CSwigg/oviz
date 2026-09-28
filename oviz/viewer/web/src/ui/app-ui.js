// The application shell: wires panels, dock, overlays, shortcuts and modes.

import { h, icon, iconButton, kbd, isEditable, MOD, downloadBlob, copyText, clear, localStorageGet, localStorageSet } from "./dom.js";
import { slider, toggle, miniSeg, installTooltips, createToasts } from "./controls.js";
import { LayersPanel } from "./layers.js";
import { TimelineDock } from "./dock.js";
import { Inspector, HoverCard } from "./inspector.js";
import { Palette } from "./palette.js";
import { lutGradient } from "../core/color.js";
import { rafThrottle } from "../core/emitter.js";
import { formatAngle, formatDistance, clamp } from "../core/math.js";
import { clonePose } from "../engine/camera.js";
import { cpuFramePosition, frameOffset } from "../engine/frames.js";
import { SkyPlugin } from "../sky/sky.js";
import { encodeViewHash, decodeViewHash } from "../app/viewhash.js";
import { StoryPlugin } from "./story.js";
import { RecorderPlugin } from "./recorder.js";
import { FilterPlugin } from "./filter.js";

export function mountUI(root, viewer) {
  const ui = new AppUI(root, viewer);
  ui.use(new SkyPlugin());
  ui.use(new FilterPlugin());
  ui.use(new StoryPlugin());
  ui.use(new RecorderPlugin());
  ui.layers.render();
  return ui;
}

export class AppUI {
  constructor(root, viewer) {
    this.root = root;
    this.viewer = viewer;
    this.manifest = viewer.manifest;
    this.luts = new Map();
    this.trails = new Map();
    this.following = null;
    this.selection = null;
    this.plugins = [];
    for (const [name, entry] of Object.entries(this.manifest.colormaps || {})) {
      const lut = viewer.store.peek(entry.blob);
      if (lut) this.luts.set(name, lut);
    }
    root.dataset.view = viewer.state.view.mode;
    this.ui = h("div", { class: "ov-ui" });
    root.append(this.ui);
    this.toasts = createToasts(this.ui);
    this.tooltips = installTooltips(root);
    this.buildTop();
    this.layers = new LayersPanel(this);
    this.inspector = new Inspector(this);
    this.dock = new TimelineDock(this);
    this.hovercard = new HoverCard(this);
    this.palette = new Palette(this);
    this.buildBottom();
    this.ui.append(this.layers.el, this.inspector.el, this.hovercard.el);
    this.bindPointer();
    this.bindKeys();
    this.bindViewer();
    this.applyTheme(localStorageGet("oviz.theme") || document.documentElement.dataset.ovizTheme || "dark");
    this.setLayersOpen(window.innerWidth > 720);
  }

  // ------------------------------------------------------------- building

  buildTop() {
    const m = this.manifest;
    const title = m.title || "Oviz figure";
    this.brand = h("div", { class: "ov-brand ov-glass" },
      h("span", { class: "ov-brand-mark", "aria-hidden": "true" }),
      h("span", { class: "ov-brand-name" }, "oviz"),
      m.title ? h("span", { class: "ov-brand-sep" }) : null,
      m.title ? h("span", { class: "ov-brand-title", title }, title) : null,
      h("button", { class: "ov-brand-info", type: "button", "aria-label": "About this figure", "data-tip": "About this figure", onclick: () => this.showAbout() }, icon("info")),
    );
    this.search = h("button", { class: "ov-search-trigger ov-glass", type: "button", "aria-label": "Search", onclick: () => this.palette.show() },
      icon("search"), h("span", { class: "ov-grow" }, "Search clusters, layers, actions…"), kbd(`${MOD === "⌘" ? "⌘" : "Ctrl"}`), kbd("K"));
    this.layersBtn = iconButton("layers", "Layers", () => this.setLayersOpen(!this.layersOpen), { shortcut: "L" });
    this.statesBtn = iconButton("bookmark", "Views & story", () => this.plugins.find((p) => p.name === "states")?.toggle(), { shortcut: "S" });
    this.shotBtn = iconButton("camera", "Screenshot", (e) => this.captureMenu(e.currentTarget), { shortcut: "P" });
    this.shareBtn = iconButton("share", "Share & export", (e) => this.shareMenu(e.currentTarget));
    this.settingsBtn = iconButton("sliders", "Display settings", (e) => this.settingsMenu(e.currentTarget));
    this.helpBtn = iconButton("help", "Shortcuts & help", () => this.showHelp(), { shortcut: "?" });
    this.toolbar = h("div", { class: "ov-toolbar ov-glass", role: "toolbar", "aria-label": "Figure tools" },
      this.layersBtn, this.statesBtn, h("span", { class: "ov-toolbar-sep" }),
      this.shotBtn, this.shareBtn, this.settingsBtn, this.helpBtn);
    this.top = h("header", { class: "ov-top ov-chrome" }, this.brand, this.search, this.toolbar);
    this.ui.append(this.top);
  }

  buildBottom() {
    const v = this.viewer;
    this.scale = h("div", { class: "ov-scale ov-glass" },
      this.scaleLabel = h("div", { class: "ov-scale-label" }),
      this.scaleBar = h("div", { class: "ov-scale-bar" }),
      this.scaleSub = h("div", { class: "ov-scale-sub" }));
    this.colorbars = h("div", { class: "ov-colorbars" });
    this.stream = h("div", { class: "ov-stream ov-glass", "data-done": "true" }, h("span", { class: "ov-spinner" }), h("span", null, "Loading volumes…"));
    const left = h("div", { class: "ov-bottom-left" }, this.stream, this.colorbars, this.scale);
    const hasSky = !!this.manifest.sky?.enabled;
    this.viewSeg = h("div", { class: "ov-seg ov-seg--view ov-glass", role: "group", "aria-label": "View" });
    this.viewThumb = h("span", { class: "ov-seg-thumb" });
    this.view3d = h("button", { type: "button", "data-value": "3d", "data-tip": "Galactic 3D view  V" }, icon("cube"), "3D");
    this.viewSky = h("button", { type: "button", "data-value": "sky", "data-tip": "Sky view from the Sun  V" }, icon("globe"), "Sky");
    this.view3d.addEventListener("click", () => this.setViewMode("3d"));
    this.viewSky.addEventListener("click", () => this.setViewMode("sky"));
    this.viewSeg.append(this.viewThumb, this.view3d, this.viewSky);
    if (!hasSky) this.viewSeg.hidden = true;
    this.homeBtn = iconButton("home", "Reset view", () => v.resetView(), { shortcut: "R", cls: "ov-glass" });
    this.fsBtn = iconButton("expand", "Fullscreen", () => this.toggleFullscreen(), { shortcut: "F", cls: "ov-glass" });
    this.skyExtras = h("div", { class: "ov-corner-btns" });
    const right = h("div", { class: "ov-bottom-right" }, this.skyExtras, this.viewSeg, h("div", { class: "ov-corner-btns" }, this.homeBtn, this.fsBtn));
    this.bottom = h("footer", { class: "ov-bottom ov-chrome" }, left, this.dock.el, right);
    this.ui.append(this.bottom);
    this.syncViewSeg();
    this.renderColorbars();
  }

  // ------------------------------------------------------------- viewer wiring

  bindViewer() {
    const v = this.viewer;
    const updateScale = rafThrottle(() => this.updateScale());
    v.on("camera", updateScale);
    v.renderer.afterRender.push(() => updateScale());
    v.on("style", () => this.renderColorbars());
    v.on("viewmode", ({ mode }) => {
      this.root.dataset.view = mode;
      this.syncViewSeg();
      this.renderColorbars();
    });
    v.on("time", () => {
      if (this.following) this.updateFollow();
    });
    new ResizeObserver(() => this.syncViewSeg()).observe(this.viewSeg);
    v.on("gpu-restored", () => {
      // Trails and measurement lines were GPU objects; drop them cleanly.
      this.trails.clear();
      this.clearMeasure();
      this.inspector.refresh();
    });
  }

  onBackgroundLoad(event) {
    const pending = (this.manifest.volumes || []).length;
    if (!pending) return;
    if (event === "start") this.stream.dataset.done = "false";
    if (event === "volumes-done") this.stream.dataset.done = "true";
    if (event === "volume") this.layers.sync();
  }

  afterBoot() {
    for (const p of this.plugins) p.afterBoot?.();
    const hv = decodeViewHash(location.hash);
    if (hv) this.applyViewHash(hv);
    this.updateScale();
  }

  async applyViewHash(hv) {
    const v = this.viewer;
    if (hv.group && (this.manifest.groups?.order || []).includes(hv.group)) this.layers.setGroup(hv.group);
    if (hv.time != null) v.timeline.setTime(hv.time);
    if (hv.mode === "sky" && this.manifest.sky?.enabled) await v.setViewMode("sky", { animate: false });
    if (hv.pose) v.setPose({ ...v.pose, ...hv.pose });
  }

  use(plugin) {
    this.plugins.push(plugin);
    plugin.attach?.(this);
    return plugin;
  }

  // ------------------------------------------------------------- pointer

  bindPointer() {
    const v = this.viewer;
    const canvas = v.canvas;
    let down = null;
    const hoverAt = rafThrottle((x, y, cx, cy) => {
      if (v.controls.pointers.size) return;
      const hit = v.pick(x, y);
      v.setHover(hit);
      this.root.dataset.hovering = hit ? "true" : "false";
      this.hovercard.show(hit, cx, cy);
      this._lastHover = { hit, cx, cy };
    });
    canvas.addEventListener("pointermove", (e) => {
      if (e.pointerType === "touch") return;
      const r = canvas.getBoundingClientRect();
      const rr = this.root.getBoundingClientRect();
      hoverAt(e.clientX - r.left, e.clientY - r.top, e.clientX - rr.left, e.clientY - rr.top);
    });
    canvas.addEventListener("pointerleave", () => {
      v.setHover(null);
      this.hovercard.show(null);
      this.root.dataset.hovering = "false";
    });
    canvas.addEventListener("pointerdown", (e) => {
      down = { x: e.clientX, y: e.clientY, t: performance.now() };
      this.hovercard.show(null);
      canvas.focus({ preventScroll: true });
    });
    canvas.addEventListener("pointerup", (e) => {
      if (!down) return;
      const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
      const quick = performance.now() - down.t < 450;
      down = null;
      if (moved > 5 || !quick) return;
      const r = canvas.getBoundingClientRect();
      const hit = v.pick(e.clientX - r.left, e.clientY - r.top, { radius: e.pointerType === "touch" ? 18 : 6 });
      if (e.shiftKey && hit && this.selection) this.measure(this.selection, hit);
      else this.select(hit);
    });
    canvas.addEventListener("dblclick", (e) => {
      const r = canvas.getBoundingClientRect();
      const hit = v.pick(e.clientX - r.left, e.clientY - r.top);
      if (hit && hit.kind !== "member") {
        this.select(hit);
        v.flyToObject(hit.trace, hit.index);
      }
    });
  }

  refreshHover() {
    const l = this._lastHover;
    if (l) this.hovercard.show(l.hit, l.cx, l.cy);
  }

  get narrow() {
    return this.root.clientWidth <= 720;
  }

  select(hit) {
    const v = this.viewer;
    this.selection = hit;
    if (hit && this.narrow) this.setLayersOpen(false);
    this.inspector.show(hit);
    const labels = v.labels;
    labels.removeDynamic("selection");
    this.clearMeasure();
    if (hit && hit.kind !== "member") {
      v.objectMeta(hit.trace).then(() => {
        if (this.selection !== hit) return;
        const d = v.describeObject(hit.trace, hit.index);
        labels.setDynamic("selection", d?.name || "", () => v.objectPosition(hit.trace, hit.index), "ov-label--sel");
        v.renderer.invalidate();
      });
    } else if (hit && hit.kind === "member") {
      labels.setDynamic("selection", `${hit.member.cluster} member`, () => hit.member.position, "ov-label--sel");
    }
    if (!hit && this.following) this.toggleFollow(null);
    v.selected = hit;
    v.renderer.invalidate();
    this.viewer.emit("select", hit);
  }

  async selectAndFly(hit, fly = true) {
    await this.viewer.objectMeta(hit.trace);
    this.select(hit);
    if (fly) this.viewer.flyToObject(hit.trace, hit.index);
  }

  // ------------------------------------------------------------- measure

  measure(a, b) {
    const v = this.viewer;
    const pa = a.kind === "member" ? a.member.position : v.objectPosition(a.trace, a.index);
    const pb = b.kind === "member" ? b.member.position : v.objectPosition(b.trace, b.index);
    if (!pa || !pb) return;
    const d = Math.hypot(pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]);
    const mid = () => {
      const qa = a.kind === "member" ? a.member.position : v.objectPosition(a.trace, a.index);
      const qb = b.kind === "member" ? b.member.position : v.objectPosition(b.trace, b.index);
      return qa && qb ? [(qa[0] + qb[0]) / 2, (qa[1] + qb[1]) / 2, (qa[2] + qb[2]) / 2] : null;
    };
    let text = formatDistance(d);
    const eye = v.skyEye();
    const ua = norm([pa[0] - eye[0], pa[1] - eye[1], pa[2] - eye[2]]);
    const ub = norm([pb[0] - eye[0], pb[1] - eye[1], pb[2] - eye[2]]);
    const ang = Math.acos(clamp(ua[0] * ub[0] + ua[1] * ub[1] + ua[2] * ub[2], -1, 1)) * 180 / Math.PI;
    text += ` · ${formatAngle(ang)}`;
    v.labels.setDynamic("measure", text, mid, "ov-label--measure");
    this.setTrailPath("__measure", [pa, pb], { color: "#7cc4ff", width: 1.5, dash: "dash", live: [a, b] });
    this.toast(`Separation ${text}`, { icon: icon("ruler") });
  }

  clearMeasure() {
    this.viewer.labels.removeDynamic("measure");
    this.removeTrail("__measure");
  }

  // ------------------------------------------------------------- trails & follow

  hasTrail(hit) {
    return this.trails.has(`${hit.trace}:${hit.index}`);
  }

  toggleTrail(hit) {
    const key = `${hit.trace}:${hit.index}`;
    if (this.trails.has(key)) { this.removeTrail(key); return; }
    const v = this.viewer;
    const trace = v.traceByKey.get(hit.trace);
    const d = v.data.get(hit.trace);
    const F = trace.points.position.frames || 1;
    const pts = [];
    const tmp = [0, 0, 0];
    // Sample finely so the Catmull-Rom path matches the moving point.
    const steps = Math.min(4000, (F - 1) * 8);
    for (let s = 0; s <= steps; s++) {
      const f = (s / steps) * (F - 1);
      const p = cpuFramePosition(d.position, F, trace.points.count, hit.index, f, tmp);
      if (!p) continue;
      const off = frameOffset(d.offset, f);
      pts.push([p[0] + off[0], p[1] + off[1], p[2] + off[2]]);
    }
    const desc = v.describeObject(hit.trace, hit.index);
    this.setTrailPath(key, pts, { color: desc?.color || trace.color, width: 1.6 });
    this.toast("Orbit trail shown for the full time range", { icon: icon("trail"), ms: 1800 });
  }

  setTrailPath(key, pts, { color, width = 1.5, dash = "solid" } = {}) {
    const v = this.viewer;
    this.removeTrail(key);
    if (pts.length < 2) return;
    const S = pts.length - 1;
    const pos = new Float32Array(S * 6);
    const arc = new Float32Array(S * 2);
    let acc = 0;
    for (let i = 0; i < S; i++) {
      const a = pts[i], b = pts[i + 1];
      pos.set(a, i * 6);
      pos.set(b, i * 6 + 3);
      arc[i * 2] = acc;
      acc += Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
      arc[i * 2 + 1] = acc;
    }
    const trace = { key, name: key, color, showInLegend: false, lines: { count: S, position: { frames: 1 }, color, width, dash, opacity: 0.9 } };
    v.lines.addTrace(trace, { position: pos, arc });
    v.extraStyles.set(key, { visible: true, presence: 1, opacityScale: 1, sizeScale: 1 });
    this.trails.set(key, trace);
    v.renderer.invalidate();
  }

  removeTrail(key) {
    if (!this.trails.has(key)) return;
    const v = this.viewer;
    v.lines.removeTrace(key);
    v.extraStyles.delete(key);
    this.trails.delete(key);
    v.renderer.invalidate();
  }

  toggleFollow(hit) {
    const v = this.viewer;
    if (!hit || (this.following && this.following.trace === hit.trace && this.following.index === hit.index)) {
      this.following = null;
      v.followOffset = null;
      return;
    }
    this.following = hit;
    const pos = v.objectPosition(hit.trace, hit.index);
    if (pos) v.flyToObject(hit.trace, hit.index, { distance: Math.min(v.pose.distance, 600) });
    this.toast("Following — the camera tracks this object through time", { icon: icon("follow") });
  }

  updateFollow() {
    const v = this.viewer;
    const f = this.following;
    if (!f || v.tween || v.state.view.mode === "sky") return;
    const pos = v.objectPosition(f.trace, f.index);
    if (!pos) return;
    v.pose.target = pos;
    v.renderer.camera.pose.target = pos;
    v.renderer.invalidate();
  }

  // ------------------------------------------------------------- modes

  setViewMode(mode) {
    if (mode === "sky" && !this.manifest.sky?.enabled) return;
    this.viewer.setViewMode(mode);
  }

  syncViewSeg() {
    const mode = this.viewer.state.view.mode;
    this.view3d.setAttribute("aria-pressed", String(mode === "3d"));
    this.viewSky.setAttribute("aria-pressed", String(mode === "sky"));
    const active = mode === "sky" ? this.viewSky : this.view3d;
    this.viewThumb.style.width = `${active.offsetWidth}px`;
    this.viewThumb.style.transform = `translateX(${active.offsetLeft}px)`;
  }

  setLayersOpen(open) {
    if (open && this.narrow && this.selection) this.select(null);
    this.layersOpen = open;
    this.layers.el.dataset.open = String(open);
    this.layersBtn.setAttribute("aria-pressed", String(open));
  }

  setZen(on) {
    this.root.dataset.zen = String(on);
  }

  toggleFullscreen() {
    const el = this.root;
    if (document.fullscreenElement) document.exitFullscreen?.();
    else (el.requestFullscreen || el.webkitRequestFullscreen)?.call(el).catch?.(() => this.toast("Fullscreen is not available here"));
  }

  applyTheme(theme) {
    this.theme = theme === "light" ? "light" : "dark";
    document.documentElement.dataset.ovizTheme = this.theme;
    const v = this.viewer;
    v.renderer.clearColor = [0, 0, 0, 0];
    v.renderer.invalidate();
    localStorageSet("oviz.theme", this.theme);
  }

  focusCanvas() {
    this.viewer.canvas.focus({ preventScroll: true });
  }

  // ------------------------------------------------------------- HUD

  updateScale() {
    const sb = this.viewer.scaleBar(110);
    if (!sb) return;
    this.scaleBar.style.width = `${Math.max(24, sb.px).toFixed(0)}px`;
    if (sb.unit === "deg") {
      this.scaleLabel.textContent = formatAngle(sb.value);
      this.scaleSub.textContent = `Field ${formatAngle(this.viewer.renderer.camera.hfov)}`;
    } else {
      this.scaleLabel.textContent = formatDistance(sb.value);
      const p = this.viewer.pose;
      this.scaleSub.textContent = `Eye ${formatDistance(p.distance)} away`;
    }
  }

  renderColorbars() {
    const v = this.viewer;
    clear(this.colorbars);
    for (const trace of v.traces) {
      const st = v.state.traces[trace.key] || {};
      const cb = trace.colorBy;
      if (!cb || !st.visible) continue;
      if ((st.colorMode || cb.defaultMode) !== "by_value") continue;
      const lut = this.luts.get(st.colormap || cb.colormap);
      if (!lut) continue;
      this.colorbars.append(h("div", { class: "ov-colorbar ov-glass" },
        h("div", { class: "ov-colorbar-title" }, h("span", { class: "ov-colorbar-dot", style: { background: trace.legendColor } }), `${trace.name} · ${cb.label}`),
        h("span", { class: "ov-colorbar-min" }, fmtNum(st.cmin ?? cb.cmin)),
        h("div", { class: "ov-colorbar-ramp", style: { background: lutGradient(lut) } }),
        h("span", { class: "ov-colorbar-max" }, fmtNum(st.cmax ?? cb.cmax))));
      if (this.colorbars.children.length >= 3) break;
    }
  }

  // ------------------------------------------------------------- menus

  menu(anchor, items, { title, align = "right", above = false } = {}) {
    this.closeMenu();
    const pop = h("div", { class: "ov-pop ov-glass", role: "menu" });
    if (title) pop.append(h("div", { class: "ov-pop-title" }, title));
    for (const it of items) {
      if (it === "-") { pop.append(h("div", { class: "ov-pop-sep" })); continue; }
      if (it.body) { pop.append(it.body); continue; }
      const b = h("button", { class: "ov-pop-item", type: "button", role: "menuitem" },
        it.icon ? icon(it.icon) : it.checked !== undefined ? h("span", { class: "ov-icon" }, it.checked ? "✓" : "") : null,
        h("span", { class: "ov-grow" }, it.label),
        it.shortcut ? kbd(it.shortcut) : null);
      b.addEventListener("click", () => { if (!it.keepOpen) this.closeMenu(); it.run?.(); });
      pop.append(b);
    }
    this.ui.append(pop);
    const r = anchor.getBoundingClientRect();
    const rr = this.root.getBoundingClientRect();
    const w = pop.offsetWidth, ht = pop.offsetHeight;
    let x = align === "center" ? r.left + r.width / 2 - w / 2 : align === "left" ? r.left : r.right - w;
    x = clamp(x - rr.left, 8, rr.width - w - 8);
    let y = above ? r.top - rr.top - ht - 8 : r.bottom - rr.top + 8;
    if (!above && y + ht > rr.height - 8) y = r.top - rr.top - ht - 8;
    pop.style.left = `${x}px`;
    pop.style.top = `${Math.max(8, y)}px`;
    anchor.setAttribute("aria-expanded", "true");
    const onDown = (e) => {
      if (!pop.contains(e.target) && e.target !== anchor && !anchor.contains(e.target)) this.closeMenu();
    };
    setTimeout(() => this.root.addEventListener("pointerdown", onDown, true), 0);
    this._menu = { pop, anchor, onDown };
    pop.querySelector("button, input, select")?.focus({ preventScroll: true });
    return pop;
  }

  closeMenu() {
    const m = this._menu;
    if (!m) return false;
    m.pop.remove();
    m.anchor.setAttribute("aria-expanded", "false");
    this.root.removeEventListener("pointerdown", m.onDown, true);
    this._menu = null;
    return true;
  }

  settingsMenu(anchor) {
    const v = this.viewer;
    const g = v.state.global;
    const body = h("div", { class: "ov-pop-body" },
      miniSeg({ label: "Theme", value: this.theme, options: [{ value: "dark", label: "Dark" }, { value: "light", label: "Light" }], onChange: (t) => this.applyTheme(t) }),
      slider({ label: "Point size", min: 0.25, max: 4, value: g.pointSize, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setGlobal({ pointSize: x }) }),
      slider({ label: "Point brightness", min: 0, max: 2, value: g.pointOpacity, format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setGlobal({ pointOpacity: x }) }),
      slider({ label: "Star glow", min: 0, max: 4, value: g.glow, format: (x) => (x <= 0.02 ? "Markers" : x.toFixed(2)), onInput: (x) => v.setGlobal({ glow: x }) }),
      slider({ label: "Birth fade", min: 0, max: 30, step: 0.5, value: g.fadeTime, format: (x) => `${x.toFixed(1)} Myr`, onInput: (x) => v.setGlobal({ fadeTime: x }) }),
      toggle({ label: "Fade by opacity", hint: "Otherwise newborn clusters grow in size", checked: g.fadeByOpacity, onChange: (x) => v.setGlobal({ fadeByOpacity: x }) }),
      toggle({ label: "Fade out after birth", checked: g.fadeInOut, onChange: (x) => v.setGlobal({ fadeInOut: x }) }),
      toggle({ label: "Size by member count", checked: g.sizeByStars, onChange: (x) => v.setGlobal({ sizeByStars: x }) }),
      toggle({ label: "Auto-orbit camera", checked: !!v.controls.autoOrbit, onChange: (x) => this.setAutoOrbit(x) }),
      slider({ label: "Field of view", min: 10, max: 100, step: 1, value: v.pose.fov, format: (x) => `${Math.round(x)}°`, onInput: (x) => { v.pose.fov = x; v.renderer.invalidate(); } }),
      toggle({ label: "Performance overlay", checked: !!this.perfEl, onChange: (x) => this.setPerfOverlay(x) }),
    );
    this.menu(anchor, [{ body }], { title: "Display" });
  }

  setAutoOrbit(on) {
    this.viewer.controls.autoOrbit = on ? 0.06 : 0;
    this.viewer.state.global.autoOrbit = on;
    this.viewer.renderer.invalidate();
  }

  captureMenu(anchor) {
    const items = [
      { label: "Save PNG", icon: "camera", shortcut: "P", run: () => this.screenshot(1) },
      { label: "Save PNG at 2× resolution", icon: "camera", run: () => this.screenshot(2) },
      { label: "Save PNG at 4× resolution", icon: "camera", run: () => this.screenshot(4) },
      { label: "Copy image to clipboard", icon: "copy", run: () => this.screenshot(1, { clipboard: true }) },
      "-",
    ];
    const rec = this.plugins.find((p) => p.name === "recorder");
    if (rec) items.push(...rec.menuItems());
    this.menu(anchor, items, { title: "Capture" });
  }

  shareMenu(anchor) {
    const items = [
      { label: "Copy link to this view", icon: "share", run: () => this.copyViewLink() },
    ];
    const states = this.plugins.find((p) => p.name === "states");
    if (states) items.push("-", ...states.exportMenuItems());
    this.menu(anchor, items, { title: "Share" });
  }

  async screenshot(scale = 1, { clipboard = false } = {}) {
    const v = this.viewer;
    const sky = this.plugins.find((p) => p.name === "sky");
    try {
      const blob = await v.renderer.capture({
        scale,
        background: sky?.captureBackground ? (ctx, w, h) => sky.captureBackground(ctx, w, h) : "#04060a",
      });
      if (clipboard && navigator.clipboard?.write) {
        await navigator.clipboard.write([new ClipboardItem({ "image/png": blob })]);
        this.toast("Image copied", { icon: icon("copy") });
      } else {
        downloadBlob(blob, `${slug(this.manifest.title || "oviz-figure")}-${stamp()}.png`);
        this.toast(`Saved ${scale > 1 ? `${scale}× ` : ""}PNG`, { icon: icon("camera") });
      }
    } catch (err) {
      console.error(err);
      this.toast("Screenshot failed");
    }
  }

  copyViewLink() {
    const url = new URL(location.href);
    url.hash = encodeViewHash(this.viewer);
    copyText(url.toString()).then(() => this.toast("Link to this view copied", { icon: icon("share") }));
  }

  // ------------------------------------------------------------- palette sources

  commands() {
    const v = this.viewer;
    const tl = v.timeline;
    const list = [
      { title: tl.playing ? "Pause" : "Play time", icon: tl.playing ? "pause" : "play", shortcut: "Space", pinned: true, run: () => this.dock.togglePlay() },
      { title: "Reset view", icon: "home", shortcut: "R", pinned: true, run: () => v.resetView() },
      { title: "Go to present day (t = 0)", icon: "target", keywords: "now zero", run: () => tl.setTime(0) },
      { title: this.layersOpen ? "Hide layers panel" : "Show layers panel", icon: "layers", shortcut: "L", run: () => this.setLayersOpen(!this.layersOpen) },
      { title: v.state.global.grid ? "Hide Galactic grid" : "Show Galactic grid", icon: "grid", shortcut: "G", run: () => v.setGlobal({ grid: !v.state.global.grid }) },
      { title: "Save screenshot", icon: "camera", shortcut: "P", run: () => this.screenshot(1) },
      { title: "Copy link to this view", icon: "share", run: () => this.copyViewLink() },
      { title: this.theme === "dark" ? "Switch to light theme" : "Switch to dark theme", icon: this.theme === "dark" ? "sun" : "moon", run: () => this.applyTheme(this.theme === "dark" ? "light" : "dark") },
      { title: v.controls.autoOrbit ? "Stop auto-orbit" : "Auto-orbit camera", icon: "orbit", shortcut: "O", run: () => this.setAutoOrbit(!v.controls.autoOrbit) },
      { title: "Hide interface (zen)", icon: "expand", shortcut: "Z", run: () => this.setZen(this.root.dataset.zen !== "true") },
      { title: "Keyboard shortcuts", icon: "keyboard", shortcut: "?", run: () => this.showHelp() },
      { title: "Fullscreen", icon: "expand", shortcut: "F", run: () => this.toggleFullscreen() },
    ];
    if (this.manifest.sky?.enabled) {
      list.splice(2, 0, { title: v.state.view.mode === "sky" ? "Switch to Galactic 3D" : "Switch to Sky view", icon: v.state.view.mode === "sky" ? "cube" : "globe", shortcut: "V", pinned: true, run: () => this.setViewMode(v.state.view.mode === "sky" ? "3d" : "sky") });
    }
    for (const p of this.plugins) list.push(...(p.commands?.() || []));
    return list;
  }

  stateItems() {
    const p = this.plugins.find((x) => x.name === "states");
    return p ? p.paletteItems() : [];
  }

  layerItems() {
    const v = this.viewer;
    const items = [];
    for (const t of v.traces) {
      if (!t.showInLegend || v.state.traces[t.key]?.inGroup === false) continue;
      const on = !!v.state.traces[t.key]?.visible;
      items.push({ title: `${on ? "Hide" : "Show"} ${t.name}`, swatch: t.legendColor, run: () => v.setTraceStyle(t.key, { visible: !on }) });
    }
    const seen = new Set();
    for (const s of this.manifest.volumes || []) {
      if (seen.has(s.stateKey)) continue;
      seen.add(s.stateKey);
      const on = !!v.state.volumes[s.stateKey]?.visible;
      items.push({ title: `${on ? "Hide" : "Show"} ${s.name}`, icon: "layers", run: () => v.setVolumeState(s.stateKey, { visible: !on }) });
    }
    return items;
  }

  cmapLabel(name) {
    const e = this.manifest.colormaps?.[name];
    return (e?.label || name).replace(/~\d+$/, "");
  }

  formatTime(t) {
    const a = Math.abs(t);
    let s = t.toFixed(a >= 100 ? 0 : 1);
    if (/^-0(\.0)?$/.test(s)) s = s.slice(1);
    return `${s.replace("-", "−")} Myr`;
  }

  frameTrace(key) {
    const v = this.viewer;
    const trace = v.traceByKey.get(key);
    const d = v.data.get(key);
    if (!trace?.points || !d) return;
    const N = trace.points.count;
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    let n = 0;
    for (let i = 0; i < N; i++) {
      const p = v.objectPosition(key, i);
      if (!p) continue;
      n++;
      for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); }
    }
    if (!n) return;
    const c = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
    const r = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2;
    const p = clonePose(v.pose);
    p.target = c;
    p.distance = Math.max(60, r / Math.tan((p.fov * Math.PI) / 360) * 1.25);
    v.animateTo(p, { duration: 1100 });
  }

  toast(text, opts) {
    return this.toasts.show(text, opts);
  }

  // ------------------------------------------------------------- overlays

  sheet(title, content, { wide = false } = {}) {
    this.closeSheet();
    const scrim = h("div", { class: "ov-scrim", onclick: () => this.closeSheet() });
    const el = h("div", { class: "ov-sheet ov-glass", role: "dialog", "aria-label": title, style: wide ? { width: "min(900px, calc(100vw - 32px))" } : null },
      h("div", { class: "ov-sheet-head" }, h("div", { class: "ov-sheet-title" }, title), iconButton("close", "Close", () => this.closeSheet(), { shortcut: "Esc" })),
      h("div", { class: "ov-sheet-body" }, content));
    this.ui.append(scrim, el);
    this._sheet = { scrim, el };
    el.querySelector("button")?.focus({ preventScroll: true });
    return el;
  }

  closeSheet() {
    if (!this._sheet) return false;
    this._sheet.scrim.remove();
    this._sheet.el.remove();
    this._sheet = null;
    this.focusCanvas();
    return true;
  }

  showHelp() {
    const groups = [
      ["Navigate", [["Orbit", "Drag"], ["Pan", "⇧ Drag"], ["Zoom toward cursor", "Scroll"], ["Fly to object", "Double-click"], ["Reset view", "R"], ["Auto-orbit", "O"], ["3D ⇄ Sky", "V"]]],
      ["Time", [["Play / pause", "Space"], ["Step frame", "← →"], ["Step 5 frames", "⇧ ← →"], ["Slower / faster", "< >"], ["Present day", "0"]]],
      ["Inspect", [["Search anything", `${MOD} K`], ["Select object", "Click"], ["Measure separation", "⇧ Click"], ["Clear selection", "Esc"]]],
      ["Layers", [["Toggle layer 1–9", "1–9"], ["Solo layer", "⇧ 1–9"], ["Layers panel", "L"], ["Galactic grid", "G"]]],
      ["Views & story", [["Views panel", "S"], ["Save current view", "N"], ["Present", "⇧ P"], ["Next / previous view", "] ["]]],
      ["Capture", [["Screenshot", "P"], ["Hide interface", "Z"], ["Fullscreen", "F"], ["This help", "?"]]],
    ];
    const grid = h("div", { class: "ov-help-grid" }, groups.map(([title, rows]) => h("div", { class: "ov-help-group" },
      h("h4", null, title),
      rows.map(([label, keys]) => h("div", { class: "ov-help-item" }, h("span", null, label),
        h("span", { class: "ov-help-keys" }, keys.split(" ").map((k) => kbd(k))))))));
    this.sheet("Shortcuts", grid, { wide: true });
  }

  showAbout() {
    const m = this.manifest;
    const about = h("div", { class: "ov-about" });
    if (m.title) about.append(h("div", { style: { font: "650 16px/1.3 var(--ov-font-display)", color: "var(--ov-text)" } }, m.title));
    if (m.note) about.append(h("div", null, String(m.note)));
    const counts = [];
    const nPts = m.traces.filter((t) => t.points && t.showInLegend).reduce((s, t) => s + t.points.count, 0);
    if (nPts) counts.push(`${nPts.toLocaleString()} objects in ${m.traces.filter((t) => t.points && t.showInLegend).length} layers`);
    if (m.volumes?.length) counts.push(`${new Set(m.volumes.map((x) => x.stateKey)).size} volumes`);
    if (m.sky?.members) counts.push(`${m.sky.members.count.toLocaleString()} member stars`);
    const T = m.time?.values || [];
    if (T.length > 1) counts.push(`${T.length} time steps from ${this.formatTime(T[0])} to ${this.formatTime(T[T.length - 1])}`);
    about.append(h("h3", null, "Contents"), h("div", null, counts.join(" · ")));
    if (m.provenance) {
      about.append(h("h3", null, "Provenance"));
      for (const [k, val] of Object.entries(m.provenance)) about.append(h("div", null, h("b", null, `${k}: `), typeof val === "string" ? val : JSON.stringify(val)));
    }
    about.append(h("h3", null, "Viewer"), h("div", null, `Made with Oviz ${m.viewer?.version || ""} — WebGL2, ${this.viewer.renderer.caps.floatRT ? "float" : "8-bit"} volume targets. Sky imagery via Aladin Lite (CDS).`));
    this.sheet("About this figure", about);
  }

  setPerfOverlay(on) {
    if (!on) { this.perfEl?.remove(); this.perfEl = null; clearInterval(this._perfTimer); return; }
    this.perfEl = h("div", { class: "ov-stream ov-glass ov-mono", style: { position: "absolute", right: "16px", top: "70px", zIndex: 40 } });
    this.ui.append(this.perfEl);
    this._perfTimer = setInterval(() => {
      const r = this.viewer.renderer;
      this.perfEl.textContent = `frame ${r.fps().toFixed(2)} ms · volumes ${this.viewer.volumes.lastMarchMs.toFixed(1)} ms · ${r.frameCount} frames`;
    }, 500);
  }

  // ------------------------------------------------------------- keys

  bindKeys() {
    const v = this.viewer;
    const tl = v.timeline;
    window.addEventListener("keydown", (e) => {
      if (isEditable(e.target)) return;
      const mod = e.metaKey || e.ctrlKey;
      const k = e.key;
      if (mod && (k === "k" || k === "K")) { this.palette.toggle(); e.preventDefault(); return; }
      if (mod) {
        for (const p of this.plugins) if (p.onKey?.(e)) { e.preventDefault(); return; }
        return;
      }
      if (e.altKey) return;
      if (this.palette.open) return;
      for (const p of this.plugins) if (p.onKey?.(e)) { e.preventDefault(); return; }
      let handled = true;
      switch (k) {
        case " ": this.dock.togglePlay(); break;
        case "ArrowLeft": tl.pause(); tl.step(e.shiftKey ? -5 : -1); break;
        case "ArrowRight": tl.pause(); tl.step(e.shiftKey ? 5 : 1); break;
        case "0": tl.setTime(0); break;
        case "<": case ",": this.dock.cycleSpeed(-1); break;
        case ">": case ".": this.dock.cycleSpeed(1); break;
        case "/": this.palette.show(); break;
        case "?": this.showHelp(); break;
        case "r": case "R": v.resetView(); break;
        case "v": case "V": this.setViewMode(v.state.view.mode === "sky" ? "3d" : "sky"); break;
        case "l": case "L": this.setLayersOpen(!this.layersOpen); break;
        case "g": case "G": v.setGlobal({ grid: !v.state.global.grid }); break;
        case "o": case "O": this.setAutoOrbit(!v.controls.autoOrbit); break;
        case "z": case "Z": this.setZen(this.root.dataset.zen !== "true"); break;
        case "f": case "F": this.toggleFullscreen(); break;
        case "p": this.screenshot(1); break;
        case "t": case "T": this.layers.toggleAll(this.manifest.traces.filter((t) => t.showInLegend && v.state.traces[t.key]?.inGroup !== false)); break;
        case "Escape":
          if (this.closeMenu() || this.closeSheet()) break;
          if (this.root.dataset.zen === "true") { this.setZen(false); break; }
          if (this.selection) { this.select(null); break; }
          handled = false;
          break;
        default:
          if (/^[1-9]$/.test(k) || /^Digit[1-9]$/.test(e.code)) {
            const n = Number(e.code.startsWith("Digit") ? e.code.slice(5) : k) - 1;
            const traces = this.manifest.traces.filter((t) => t.showInLegend && v.state.traces[t.key]?.inGroup !== false);
            const t = traces[n];
            if (t) {
              if (e.shiftKey) this.layers.solo(t.key);
              else v.setTraceStyle(t.key, { visible: !v.state.traces[t.key]?.visible });
            }
          } else handled = false;
      }
      if (handled) e.preventDefault();
    });
  }
}

// ---------------------------------------------------------------- helpers

function norm(v) {
  const l = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / l, v[1] / l, v[2] / l];
}

function fmtNum(x) {
  if (!Number.isFinite(x)) return "—";
  const a = Math.abs(x);
  return (a >= 100 ? x.toFixed(0) : a >= 10 ? x.toFixed(1) : x.toFixed(2)).replace("-", "−");
}

function slug(s) {
  return String(s).trim().toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "oviz";
}

function stamp() {
  const d = new Date();
  const p = (n) => String(n).padStart(2, "0");
  return `${d.getFullYear()}${p(d.getMonth() + 1)}${p(d.getDate())}-${p(d.getHours())}${p(d.getMinutes())}${p(d.getSeconds())}`;
}
