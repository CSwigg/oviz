// The application shell: wires panels, dock, overlays, shortcuts and modes.

import { h, icon, iconButton, kbd, isEditable, controlConsumesKey, MOD, downloadBlob, copyTextLater, clear, localStorageGet, localStorageSet } from "./dom.js";
import { slider, toggle, miniSeg, installTooltips, createToasts } from "./controls.js";
import { LayersPanel } from "./layers.js";
import { TimelineDock } from "./dock.js";
import { Inspector } from "./inspector.js";
import { Palette } from "./palette.js";
import { lutGradient } from "../core/color.js";
import { rafThrottle } from "../core/emitter.js";
import { formatAngle, formatDistance, clamp } from "../core/math.js";
import { clonePose } from "../engine/camera.js";
import { anchorLabel } from "../engine/anchor.js";
import { autoOrbitRate, initialViewerState, skyStartView } from "../app/state.js";
import { keyMotion } from "../engine/controls.js";
import { cpuFramePosition, frameOffset } from "../engine/frames.js";
import { SkyPlugin } from "../sky/sky.js";
import { encodeViewHash, decodeViewHash, encodeStatePart, decodeStatePart, encodeViewsPart, decodeViewsPart } from "../app/viewhash.js";
import { captureState } from "../app/states.js";
import { StoryPlugin, figureId } from "./story.js";
import { RecorderPlugin } from "./recorder.js";
import { FilterPlugin } from "./filter.js";
import { NotesPlugin } from "./notes.js";
import { LassoPlugin } from "./lasso.js";
import { SelectionReticle } from "./reticle.js";
import { MODES, modeById, initialMode, applyMode } from "./modes.js";
import { reveal, conceal, SPRINGS, flickToDismiss } from "./motion.js";
import { LayoutManager } from "./layout.js";
import { CompactLegend } from "./legend.js";
import { IdleFade } from "./idle.js";
import { buildArModel } from "../ar/model.js";
import { WidgetHost } from "./widgets.js";
import { BirthTreePlugin } from "./birthtree.js";
import { SfhPlugin } from "./sfh.js";
import { ActionsPlugin } from "./actions.js";
import { AxesPlugin } from "./axes.js";

export function mountUI(root, viewer) {
  const ui = new AppUI(root, viewer);
  if (viewer.manifest.sky?.enabled) ui.use(new SkyPlugin());
  ui.use(new FilterPlugin());
  ui.use(new NotesPlugin());
  ui.use(new LassoPlugin());
  ui.use(new StoryPlugin());
  ui.use(new RecorderPlugin());
  ui.use(new BirthTreePlugin());
  ui.use(new SfhPlugin());
  ui.use(new ActionsPlugin());
  ui.use(new AxesPlugin());
  ui.layers.render();
  ui.applyLayout();
  return ui;
}

export class AppUI {
  constructor(root, viewer) {
    this.root = root;
    this.viewer = viewer;
    this.manifest = viewer.manifest;
    this.luts = new Map();
    this.trails = new Map();
    this.selection = null;
    this.plugins = [];
    for (const [name, entry] of Object.entries(this.manifest.colormaps || {})) {
      const lut = viewer.store.peek(entry.blob);
      if (lut) this.luts.set(name, lut);
    }
    root.dataset.view = viewer.state.view.mode;
    this.mode = applyMode(initialMode());
    this.ui = h("div", { class: "ov-ui" });
    root.append(this.ui);
    this.toasts = createToasts(this.ui);
    this.tooltips = installTooltips(root);
    this.buildTop();
    this.layers = new LayersPanel(this);
    this.inspector = new Inspector(this);
    this.dock = new TimelineDock(this);
    this.palette = new Palette(this);
    this.buildBottom();
    this.ui.append(this.layers.el, this.inspector.el);
    this.reticle = new SelectionReticle(this);
    this.buildExtras();
    // Floating tool panels that plugins register (Relative SFH, Birth tree…).
    this.widgets = new WidgetHost(this);
    this.widgets.onChange(() => { this.widgetsBtn.hidden = !this.widgets.available().length; });
    viewer.stateExtensions.set("widgets", { capture: () => this.widgets.capture(), apply: (s) => this.widgets.apply(s) });
    this.layout = new LayoutManager(this);
    this.overlays = new Set();
    this.idle = new IdleFade(this);
    // States restore the selection and the panels too.
    viewer.stateExtensions.set("ui", { capture: () => this.captureUiState(), apply: (s) => this.applyUiState(s) });
    // Phones: every bottom sheet swipes down to dismiss.
    const sheetHead = (el) => el.querySelector(".ov-panel-head");
    flickToDismiss(this.inspector.el, sheetHead(this.inspector.el), { axis: "y", sign: 1, enabled: () => this.narrow, onDismiss: () => this.select(null) });
    flickToDismiss(this.layers.el, sheetHead(this.layers.el), { axis: "y", sign: 1, enabled: () => this.narrow, onDismiss: () => this.setLayersOpen(false) });
    this.bindPointer();
    this.bindKeys();
    this.bindViewer();
    // The authored theme, unless someone switched theme on this figure.
    this.themeKey = `oviz.theme.${figureId(this.manifest)}`;
    document.documentElement.dataset.ovizFigureTheme ||= document.documentElement.dataset.ovizTheme || "dark";
    this.applyTheme(localStorageGet(this.themeKey) || document.documentElement.dataset.ovizTheme || "dark", { remember: false });
    this.setLayersOpen(window.innerWidth > 720 && this.modeConfig.layersOpen);
  }

  get modeConfig() {
    return modeById(this.mode);
  }

  /** Controls that only some layouts place (kept in a hidden holder). */
  buildExtras() {
    const v = this.viewer;
    this.viewToggle = iconButton("globe", "Sky view", () => this.setViewMode(v.state.view.mode === "sky" ? "3d" : "sky"), { cls: "ov-view-toggle", shortcut: "V" });
    this.moreBtn = iconButton("more", "More", (e) => this.moreMenu(e.currentTarget));
    this.modeBtn = iconButton("sidebar", "Detailed mode", () => this.toggleMode(), { cls: "ov-mode-btn", shortcut: "U" });
    this.legend = new CompactLegend(this);
    this.spare = h("div", { class: "ov-spare", hidden: true }, this.viewToggle, this.moreBtn, this.modeBtn, this.legend.el);
    this.ui.append(this.spare);
    this.syncViewToggle();
  }

  /** Arrange the components for the current mode and screen size. */
  applyLayout() {
    this.layout.apply(this.modeConfig, this.narrow);
    this.root.dataset.layout = this.layout.current.layout;
    this.modeBtn.setAttribute("aria-pressed", String(this.mode === "detailed"));
    this.syncViewToggle();
    this.legend.render();
    this.refit();
    requestAnimationFrame(() => {
      this.syncViewSeg();
      this.anchorLayers();
      this.reticle.update();
    });
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
    this.layersBtn = iconButton("layers", "Layers", () => this.setLayersOpen(!this.layersOpen), { shortcut: "⇧L" });
    this.statesBtn = iconButton("bookmark", "Views & story", () => this.plugins.find((p) => p.name === "states")?.toggle(), { shortcut: "Y" });
    this.widgetsBtn = iconButton("widgets", "Widgets", (e) => this.widgetsMenu(e.currentTarget));
    this.widgetsBtn.hidden = true;
    this.shotBtn = iconButton("camera", "Screenshot", (e) => this.captureMenu(e.currentTarget), { shortcut: "I" });
    this.shareBtn = iconButton("share", "Share & export", (e) => this.shareMenu(e.currentTarget));
    this.settingsBtn = iconButton("sliders", "Display settings", (e) => this.settingsMenu(e.currentTarget));
    this.helpBtn = iconButton("help", "Shortcuts & help", () => this.showHelp(), { shortcut: "?" });
    this.toolbar = h("div", { class: "ov-toolbar ov-glass", role: "toolbar", "aria-label": "Figure tools" },
      this.layersBtn, this.statesBtn, this.widgetsBtn, h("span", { class: "ov-toolbar-sep" }),
      this.shotBtn, this.shareBtn, this.settingsBtn, this.helpBtn);
    // Touch screens have no Z key: one tap hides the interface (the
    // "Show interface" pill brings it back).
    this.zenBtn = iconButton("eyeOff", "Hide interface", () => this.setZen(true), { cls: "ov-glass ov-zen-btn" });
    this.top = h("header", { class: "ov-top ov-chrome" }, this.brand, this.zenBtn, this.search, this.toolbar);
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
    // The scale's row also holds the Sky imagery credit.
    this.scaleRow = h("div", { class: "ov-scale-row" }, this.scale);
    const left = h("div", { class: "ov-bottom-left" }, this.stream, this.colorbars, this.scaleRow);
    this.bottomLeft = left;
    const hasSky = (this.hasSky = !!this.manifest.sky?.enabled);
    this.viewSeg = h("div", { class: "ov-seg ov-seg--view ov-glass", role: "group", "aria-label": "View" });
    this.viewThumb = h("span", { class: "ov-seg-thumb" });
    this.view3d = h("button", { type: "button", "data-value": "3d", "data-tip": "Galactic 3D view  V" }, icon("cube"), "3D");
    this.viewSky = h("button", { type: "button", "data-value": "sky", "data-tip": "Sky view from the Sun  V" }, icon("globe"), "Sky");
    this.view3d.addEventListener("click", () => this.setViewMode("3d"));
    this.viewSky.addEventListener("click", () => this.setViewMode("sky"));
    this.viewSeg.append(this.viewThumb, this.view3d, this.viewSky);
    if (!hasSky) this.viewSeg.hidden = true;
    this.homeBtn = iconButton("home", "Reset view", () => v.resetView(), { shortcut: "Home", cls: "ov-glass" });
    this.fsBtn = iconButton("expand", "Fullscreen", () => this.toggleFullscreen(), { shortcut: "M", cls: "ov-glass" });
    this.fsBtn.hidden = !this.canFullscreen;
    this.skyExtras = h("div", { class: "ov-corner-btns ov-sky-extras" });
    this.cornerBtns = h("div", { class: "ov-corner-btns" }, this.homeBtn, this.fsBtn);
    const right = h("div", { class: "ov-bottom-right" }, this.skyExtras, this.viewSeg, this.cornerBtns);
    this.bottomRight = right;
    this.bottom = h("footer", { class: "ov-bottom ov-chrome" }, left, this.dock.el, right);
    this.ui.append(this.bottom);
    this.syncViewSeg();
    this.renderColorbars();
  }

  // ------------------------------------------------------------- viewer wiring

  bindViewer() {
    const v = this.viewer;
    // Classic Oviz: any drag, scroll, key or reset takes the camera and
    // stops auto-orbit.
    v.on("user-camera", () => { if (v.state.global.autoOrbit) this.setAutoOrbit(false); });
    const updateScale = rafThrottle(() => this.updateScale());
    v.on("camera", updateScale);
    v.renderer.afterRender.push(() => updateScale());
    v.on("style", () => this.renderColorbars());
    const colorbars = rafThrottle(() => this.renderColorbars());
    v.on("volume", colorbars);
    v.on("time", colorbars);
    v.on("volume-ready", colorbars);
    v.on("viewmode", ({ mode }) => {
      // Member stars exist only in Sky view: a member selection ends with it.
      if (mode !== "sky" && this.selection?.kind === "member") this.select(null, { recordUndo: false });
      this.root.dataset.view = mode;
      this.syncViewToggle();
      this.syncViewSeg();
      this.renderColorbars();
      // The key drops the Sun in Sky view (the eye is at the Sun).
      this.legend.render();
    });
    // The measurement line and its label follow both objects through time.
    const measure = rafThrottle(() => this.updateMeasure());
    v.on("time", () => {
      if (this._measure) measure();
    });
    // The camera anchor shows in the scale readout and the inspector.
    v.on("anchor", ({ released } = {}) => {
      this.updateScale();
      if (this.selection) this.inspector.refresh();
      // The first time a pan or flight frees the camera, say how to get back.
      const home = v.defaultAnchor;
      if (released && !this._anchorHinted && home && home.kind !== "free") {
        this._anchorHinted = true;
        this.toast(`Camera free · Home re-anchors it to ${anchorLabel(home, this.anchorName(home))}`, { icon: icon("home"), ms: 3200 });
      }
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
    // First, so the opening State the plugins capture is the Sky view.
    this.applySkyStart();
    for (const p of this.plugins) p.afterBoot?.();
    // Figures that open orbiting (e.g. the website background profile).
    if (this.viewer.state.global.autoOrbit) this.setAutoOrbit(true);
    this.viewer.controls.speed = this.viewer.state.global.navSpeed ?? 1;
    const hv = decodeViewHash(location.hash);
    if (hv) this.applyViewHash(hv);
    else this.introFlight();
    this.updateScale();
    this.firstRunHint();
    let narrow = this.narrow;
    window.addEventListener("resize", rafThrottle(() => {
      if (this.narrow !== narrow) {
        narrow = this.narrow;
        this.applyLayout();
        this.setLayersOpen(!narrow && this.modeConfig.layersOpen);
      }
      this.anchorLayers();
    }));
  }

  /**
   * Figures saved looking out from the Sun (the classic "earth" view, e.g.
   * a sky-start export) open in Sky with that gaze and field of view. Home
   * and V lead to the 3D view the figure names for leaving Sky.
   */
  applySkyStart() {
    const v = this.viewer;
    const start = skyStartView(this.manifest, v.state);
    if (!start) return;
    const back = start.returnPose ? clonePose(start.returnPose) : null;
    if (back) v.homePose = clonePose(back);
    v.setViewMode("sky", { animate: false });
    v.setPose(start.pose);
    if (back) v.returnPose = back;
  }

  /** One gentle how-to line the first time someone opens a figure. */
  firstRunHint() {
    const story = this.plugins.find((p) => p.name === "states");
    if (story?.readOnly || localStorageGet("oviz.hint.v1")) return;
    localStorageSet("oviz.hint.v1", "1");
    const touch = window.matchMedia?.("(pointer: coarse)").matches;
    // Only figures with clusters (objects with ages) have details to click.
    const sun = this.manifest.world?.sunTrace;
    const clusters = this.manifest.traces.some((t) => t.showInLegend && t.key !== sun && t.points?.ageNow?.blob);
    const text = touch
      ? `Drag to orbit · pinch to zoom${clusters ? " · tap a cluster for details" : ""}`
      : `Drag to orbit · scroll to zoom${clusters ? " · click a cluster for details" : ""}`;
    setTimeout(() => this.toast(text, { icon: icon("info"), ms: 5200 }), 1600);
  }

  /**
   * Open on a short dolly into the home view, so the scene reads as 3D from
   * the first second. Skipped for present-only files, shared links, Sky
   * starts and reduced motion; any input takes the camera over at once.
   */
  introFlight() {
    const v = this.viewer;
    const story = this.plugins.find((p) => p.name === "states");
    if (story?.readOnly || this.manifest.states?.apply_first_on_load || v.state.view.mode !== "3d") return;
    if (window.matchMedia?.("(prefers-reduced-motion: reduce)").matches) return;
    const home = clonePose(v.pose);
    const start = clonePose(home);
    start.distance = home.distance * 1.55;
    start.yaw = home.yaw + 0.38;
    start.pitch = home.pitch * 0.8;
    v.setPose(start);
    v.animateTo(home, { duration: 2300, ease: (x) => 1 - Math.pow(1 - x, 3) });
  }

  async applyViewHash(hv) {
    const v = this.viewer;
    const story = this.plugins.find((p) => p.name === "states");
    // A link with saved views brings the sender's story; a presentation
    // link opens presenting it.
    if (hv.views && story?.homeState) {
      const shared = await decodeViewsPart(hv.views, story.homeState).catch(() => null);
      if (shared?.items?.length) {
        story.loadShared(shared);
        if (hv.present) {
          story.present(true, { startAt: Math.min(hv.present, shared.items.length) - 1, instant: true });
          return;
        }
      }
    }
    // A link with the whole State: rebuilt over the figure's opening State
    // and applied exactly, like a saved view.
    if (hv.state && story?.homeState) {
      const state = await decodeStatePart(hv.state, story.homeState).catch(() => null);
      if (state?.view) {
        await story.applyState(state, { instant: true });
        return;
      }
    }
    if (hv.group && (this.manifest.groups?.order || []).includes(hv.group)) this.layers.setGroup(hv.group);
    if (hv.time != null) v.timeline.setTime(hv.time);
    if (hv.mode === "sky" && this.manifest.sky?.enabled) await v.setViewMode("sky", { animate: false });
    // A 3D link into a figure that opens in Sky.
    else if (hv.mode === "3d" && v.state.view.mode === "sky") await v.setViewMode("3d", { animate: false });
    if (hv.pose) {
      // The link's anchor (or, for older links, the one its pose implies);
      // a centre that is not on it leaves the camera free.
      const anchor = hv.anchor && v.anchorAvailable(hv.anchor) ? hv.anchor : v.inferAnchor(hv.pose);
      v.setAnchor(anchor, { fly: false });
      v.setPose({ ...v.pose, ...hv.pose });
    }
    // The layers it showed (after the group preset), only if the list still matches.
    const items = this.layers.toggleItems("all");
    if (hv.layers?.length === items.length) items.forEach((it, i) => this.layers.setItemVisible(it, hv.layers[i]));
    // Its lasso selection: the objects by trace index, and the outline that clips the dust.
    const lasso = this.plugins.find((p) => p.name === "lasso");
    if (lasso && (hv.selection || hv.mask)) {
      const selection = {};
      for (const [ti, idx] of hv.selection || []) {
        const t = v.traces[ti];
        const ok = t?.points ? idx.filter((i) => i < t.points.count) : [];
        if (ok.length) selection[t.key] = ok;
      }
      lasso.restore({ selection, isolate: hv.isolate !== false, filterOff: !!hv.filterOff, mask: hv.mask || null });
    }
    // Its selected object.
    const t = hv.object ? v.traces[hv.object.trace] : null;
    if (t?.points && hv.object.index < t.points.count) {
      await v.objectMeta(t.key);
      this.select({ trace: t.key, index: hv.object.index }, { recordUndo: false });
    }
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
    // Hovering only highlights an object (and shows a pointer); its name
    // tag appears when it is clicked, so moving over the figure stays quiet.
    const hoverAt = rafThrottle((x, y) => {
      if (v.controls.pointers.size) return;
      const hit = v.pick(x, y);
      v.setHover(hit);
      this.root.dataset.hovering = hit ? "true" : "false";
    });
    canvas.addEventListener("pointermove", (e) => {
      // How far the press has wandered: a drag that comes back (a closed
      // lasso loop) is not a click.
      if (down) down.far = Math.max(down.far, Math.hypot(e.clientX - down.x, e.clientY - down.y));
      if (e.pointerType === "touch") return;
      const r = canvas.getBoundingClientRect();
      const rr = this.root.getBoundingClientRect();
      hoverAt(e.clientX - r.left, e.clientY - r.top, e.clientX - rr.left, e.clientY - rr.top);
    });
    // Keys work as soon as the pointer is over the figure (classic Oviz; it
    // matters in iframes and notebooks), unless a field is being typed in.
    canvas.addEventListener("pointerenter", (e) => {
      if (e.pointerType === "touch") return;
      const a = document.activeElement;
      if (a && a !== document.body && a !== canvas && (isEditable(a) || this.root.contains(a))) return;
      canvas.focus({ preventScroll: true });
    });
    canvas.addEventListener("pointerleave", () => {
      v.setHover(null);
      this.root.dataset.hovering = "false";
    });
    canvas.addEventListener("pointerdown", (e) => {
      down = { x: e.clientX, y: e.clientY, t: performance.now(), far: 0 };
      canvas.focus({ preventScroll: true });
    });
    canvas.addEventListener("pointerup", (e) => {
      if (!down) return;
      const moved = Math.max(down.far, Math.hypot(e.clientX - down.x, e.clientY - down.y));
      const quick = performance.now() - down.t < 450;
      down = null;
      if (moved > 5 || !quick) return;
      // A click on the figure dismisses a floating layers panel, like a menu.
      if (this.layersOpen && this.layersFloating) this.setLayersOpen(false);
      const r = canvas.getBoundingClientRect();
      const hit = v.pick(e.clientX - r.left, e.clientY - r.top, { radius: e.pointerType === "touch" ? 18 : 6 });
      if (e.shiftKey && hit && this.selection) this.measure(this.selection, hit);
      else this.select(hit);
    });
    // Double-click sets the camera's frame, keeping the zoom. On a cluster
    // the camera centres on it and moves with it through time; anywhere
    // else it goes back to the LSR, centred on the dust, line or Milky Way
    // plane under the pointer, or on the LSR itself over empty space. In
    // Sky view a double-click turns to an object and does nothing on empty sky.
    canvas.addEventListener("dblclick", (e) => {
      const r = canvas.getBoundingClientRect();
      const x = e.clientX - r.left, y = e.clientY - r.top;
      const hit = v.pick(x, y);
      const sky = v.state.view.mode === "sky";
      if (hit && hit.kind !== "member") {
        this.select(hit);
        if (sky) v.flyToObject(hit.trace, hit.index);
        else this.followObject(hit);
        return;
      }
      if (sky || hit) return;
      this.returnToLsr(v.pickSurface(x, y));
    });
  }

  get narrow() {
    return this.root.clientWidth <= 720;
  }

  select(hit, { recordUndo = true } = {}) {
    const v = this.viewer;
    if (recordUndo && hit !== this.selection && (hit || this.selection)) this.pushSelectionUndo();
    this.selection = hit;
    if (hit && this.narrow) this.setLayersOpen(false);
    // Focus mode opens the inspector only on request (and keeps it in step
    // with the selection once it is open).
    const detailsOpen = this.inspector.el.dataset.open === "true";
    this.inspector.show(this.modeConfig.selection === "panel" || (hit && detailsOpen) ? hit : null);
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
      labels.setDynamic("selection", `${hit.member.cluster} member`, () => this.selectionPosition(hit), "ov-label--sel");
    }
    // The camera's frame outlives the selection: a followed cluster stays
    // followed until a double-click elsewhere returns to the LSR.
    v.selected = hit;
    v.renderer.invalidate();
    this.viewer.emit("select", hit);
  }

  /**
   * The interface's part of a State: the selection, plus the panels where
   * they differ from the mode's own defaults (so a view saved in Detailed,
   * where the layers panel is open anyway, does not pop it open in Focus).
   */
  captureUiState() {
    const sel = this.selection;
    const detailsOpen = this.inspector.el.dataset.open === "true";
    const figureTheme = document.documentElement.dataset.ovizFigureTheme || "dark";
    return {
      selection: !sel ? null : sel.kind === "member" ? { kind: "member", index: sel.index } : { trace: sel.trace, index: sel.index },
      layers: this.narrow || !!this.layersOpen === !!this.modeConfig.layersOpen ? null : !!this.layersOpen,
      details: this.modeConfig.selection === "callout" && detailsOpen ? true : null,
      // Classic States restored the theme and zen mode too.
      theme: this.theme !== figureTheme ? this.theme : null,
      zen: this.root.dataset.zen === "true" && this.root.dataset.presenting !== "true" ? true : null,
    };
  }

  async applyUiState(s) {
    if (!s || typeof s !== "object") return; // older States leave the interface alone
    const v = this.viewer;
    const sel = s.selection;
    let hit = null;
    if (sel?.kind === "member") {
      if (v.state.view.mode === "sky" && this.sky) {
        await this.sky.ensureMembers?.();
        hit = this.sky.resolveMemberPick?.({ index: sel.index }) || null;
      }
    } else if (sel && v.traceByKey.has(sel.trace) && Number.isInteger(sel.index)) {
      hit = { trace: sel.trace, index: sel.index };
    } else if (sel?.name) {
      hit = await this.findObject(sel.traceName, sel.name);
    }
    const same = (a, b) => (!a && !b) || (!!a && !!b && (a.kind || "") === (b.kind || "") && a.trace === b.trace && a.index === b.index);
    if (!same(hit, this.selection)) this.select(hit, { recordUndo: false });
    if (!this.narrow) this.setLayersOpen(s.layers == null ? !!this.modeConfig.layersOpen : !!s.layers);
    if (this.modeConfig.selection === "callout") {
      if (hit && s.details) this.openDetails();
      else if (this.inspector.el.dataset.open === "true") this.inspector.show(null);
    }
    if (s.theme === "light" || s.theme === "dark") this.applyTheme(s.theme, { remember: false });
    if (this.root.dataset.presenting !== "true") this.setZen(!!s.zen);
  }

  /** Find an object by its (trace and) name, as classic States store it. */
  async findObject(traceName, name) {
    const v = this.viewer;
    const norm = (x) => String(x || "").replace(/_/g, " ").trim().toLowerCase();
    const want = norm(name);
    const pointTraces = v.traces.filter((t) => t.points);
    const named = pointTraces.filter((t) => traceName && (t.name === traceName || t.key === traceName));
    for (const t of named.length ? named : pointTraces) {
      const meta = await v.objectMeta(t.key);
      const i = (meta.name || []).findIndex((n) => norm(n) === want);
      if (i >= 0) return { trace: t.key, index: i };
    }
    return null;
  }

  /** Where a selection is now (member stars move through time too). */
  selectionPosition(hit) {
    if (!hit) return null;
    if (hit.kind === "member") return this.sky?.members?.starPosition(hit.index, this.viewer.timeline.time, true) || hit.member.position;
    return this.viewer.objectPosition(hit.trace, hit.index);
  }

  async selectAndFly(hit, fly = true) {
    await this.viewer.objectMeta(hit.trace);
    this.select(hit);
    if (fly) this.viewer.flyToObject(hit.trace, hit.index);
  }

  // ------------------------------------------------------------- measure

  measure(a, b) {
    this._measure = { a, b, key: "" };
    const text = this.updateMeasure();
    if (text) this.toast(`Separation ${text}`, { icon: icon("ruler") });
  }

  /** Redraw the measurement for the current time; returns its label. */
  updateMeasure() {
    const m = this._measure;
    if (!m) return "";
    const v = this.viewer;
    const pa = this.selectionPosition(m.a), pb = this.selectionPosition(m.b);
    if (!pa || !pb) return "";
    const key = [...pa, ...pb].map((x) => x.toFixed(3)).join(",");
    if (key === m.key) return m.text;
    m.key = key;
    const d = Math.hypot(pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]);
    // The angle as seen from the Sun (from the Sky camera in Sky view).
    const eye = v.state.view.mode === "sky" ? v.skyEye() : v.sunPosition();
    const ua = norm([pa[0] - eye[0], pa[1] - eye[1], pa[2] - eye[2]]);
    const ub = norm([pb[0] - eye[0], pb[1] - eye[1], pb[2] - eye[2]]);
    const ang = Math.acos(clamp(ua[0] * ub[0] + ua[1] * ub[1] + ua[2] * ub[2], -1, 1)) * 180 / Math.PI;
    m.text = `${formatDistance(d)} · ${formatAngle(ang)}`;
    const mid = () => {
      const qa = this.selectionPosition(m.a), qb = this.selectionPosition(m.b);
      return qa && qb ? [(qa[0] + qb[0]) / 2, (qa[1] + qb[1]) / 2, (qa[2] + qb[2]) / 2] : null;
    };
    v.labels.setDynamic("measure", m.text, mid, "ov-label--measure");
    this.setTrailPath("__measure", [pa, pb], { color: "#7cc4ff", width: 1.5, dash: "dash" });
    return m.text;
  }

  clearMeasure() {
    this._measure = null;
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

  /** The object the camera is anchored to (inspector "Follow"), if any. */
  get following() {
    const a = this.viewer.anchor;
    if (a.kind === "object") return { trace: a.trace, index: a.index };
    if (a.kind === "sun" && this.viewer.world.sunTrace) return { trace: this.viewer.world.sunTrace, index: 0 };
    return null;
  }

  /** Anchor the camera to an object (or stop): it moves with it through time. */
  toggleFollow(hit) {
    const v = this.viewer;
    const f = this.following;
    if (!hit || (f && f.trace === hit.trace && f.index === hit.index)) {
      if (f) this.stopFollowing();
      return;
    }
    const sun = hit.trace === v.world.sunTrace && hit.index === 0;
    v.setAnchor(sun ? { kind: "sun" } : { kind: "object", trace: hit.trace, index: hit.index }, { distance: Math.min(v.pose.distance, 600) });
    this.toast("Following · the camera moves with it through time", { icon: icon("follow") });
  }

  /**
   * Double-click on a cluster: its frame. The camera centres on it at the
   * same zoom and moves with it through time.
   */
  followObject(hit) {
    const v = this.viewer;
    const sun = hit.trace === v.world.sunTrace && hit.index === 0;
    const a = sun ? { kind: "sun" } : { kind: "object", trace: hit.trace, index: hit.index };
    if (!v.anchorAvailable(a)) { v.recenter(v.objectPosition(hit.trace, hit.index)); return; }
    const was = v.anchor;
    v.setAnchor(a);
    if (was.kind !== a.kind || was.trace !== a.trace || was.index !== a.index) {
      this.toast(`Following ${anchorLabel(a, this.anchorName(a))} · double-click empty space to return to the LSR`, { icon: icon("follow"), ms: 2600 });
    }
  }

  /**
   * Double-click off a cluster: back to the LSR frame, centred on `point`
   * (dust or a line under the pointer) or on the LSR itself. Figures
   * without an LSR re-centre on the point, or reset the view.
   */
  returnToLsr(point) {
    const v = this.viewer;
    const lsr = { kind: "lsr" };
    if (!v.anchorAvailable(lsr)) {
      if (point) v.recenter(point);
      else v.resetView();
      return;
    }
    const was = v.anchor.kind;
    v.setAnchor(lsr, { fly: !point });
    if (point) v.recenter(point);
    if (was !== "lsr") this.toast("Back to the LSR frame", { icon: icon("target"), ms: 1400 });
  }

  /** Stop moving with an object or layer: the LSR frame again, where the camera is. */
  stopFollowing() {
    const v = this.viewer;
    if (v.anchorAvailable({ kind: "lsr" })) v.setAnchor({ kind: "lsr" }, { fly: false });
    else v.releaseAnchor();
  }

  /** Anchor the camera to the LSR, the Sun or nothing ("free"). */
  setCameraAnchor(kind) {
    const v = this.viewer;
    if (kind === "free") v.releaseAnchor();
    else if (v.anchorAvailable({ kind })) v.setAnchor({ kind });
    else return;
    const name = anchorLabel(v.anchor, this.anchorName(v.anchor));
    this.toast(kind === "free" ? "Camera free · scroll zooms toward the pointer" : `Camera anchored to ${name}`, { icon: icon(kind === "free" ? "orbit" : "target"), ms: 1600 });
  }

  /** Display name of an object or layer anchor (empty for the others). */
  anchorName(a) {
    if (a?.kind === "layer") return this.viewer.traceByKey.get(a.trace)?.name || "";
    if (a?.kind !== "object") return "";
    return this.viewer.describeObject(a.trace, a.index)?.name || "";
  }

  /** Follow a whole layer: the camera moves with its median (classic focus group). */
  followLayer(traceKey) {
    const v = this.viewer;
    const a = { kind: "layer", trace: traceKey };
    if (v.anchor.kind === "layer" && v.anchor.trace === traceKey) {
      this.stopFollowing();
      this.toast("Stopped following the layer", { icon: icon("follow"), ms: 1400 });
      return;
    }
    if (!v.anchorAvailable(a)) return;
    v.setAnchor(a);
    this.toast(`Following ${anchorLabel(a, this.anchorName(a))} · the camera moves with the group through time`, { icon: icon("follow"), ms: 2400 });
  }

  anchorField() {
    const v = this.viewer;
    const a = v.anchor;
    const options = [];
    if (v.anchorAvailable({ kind: "lsr" })) options.push({ value: "lsr", label: "LSR" });
    if (v.anchorAvailable({ kind: "sun" })) options.push({ value: "sun", label: "Sun" });
    options.push({ value: "free", label: "Free" });
    if (a.kind === "object") options.push({ value: "object", label: "Object" });
    if (a.kind === "layer") options.push({ value: "layer", label: "Layer" });
    return miniSeg({
      label: "Camera anchor",
      hint: "What the camera orbits, zooms toward and moves with through time",
      value: a.kind,
      options,
      onChange: (k) => { if (k !== "object" && k !== "layer") this.setCameraAnchor(k); },
    });
  }

  // ------------------------------------------------------------- modes

  setViewMode(mode) {
    if (mode === "sky" && !this.manifest.sky?.enabled) return Promise.resolve();
    return this.viewer.setViewMode(mode);
  }

  syncViewSeg() {
    const mode = this.viewer.state.view.mode;
    this.view3d.setAttribute("aria-pressed", String(mode === "3d"));
    this.viewSky.setAttribute("aria-pressed", String(mode === "sky"));
    const active = mode === "sky" ? this.viewSky : this.view3d;
    this.viewThumb.style.width = `${active.offsetWidth}px`;
    this.viewThumb.style.transform = `translateX(${active.offsetLeft}px)`;
  }

  /** Layers float over the figure (Focus popover, phone sheet) instead of docking. */
  get layersFloating() {
    return this.narrow || this.layout?.current?.layout === "focus";
  }

  setLayersOpen(open) {
    if (open && this.narrow && this.selection) this.select(null);
    // A floating layers panel and the story drawer share the space above the bar.
    if (open && this.layersFloating && this.root.dataset.story === "true") {
      this.plugins.find((p) => p.name === "states")?.toggle(false);
    }
    this.layersOpen = open;
    this.layers.el.dataset.open = String(open);
    this.root.dataset.layers = String(open);
    this.idle?.wake();
    this.layersBtn.setAttribute("aria-pressed", String(open));
    if (open) this.anchorLayers();
  }

  /** Focus opens the layers as a popover above its button in the bar. */
  anchorLayers() {
    const el = this.layers.el;
    if (this.layout?.current?.layout !== "focus" || this.narrow) {
      el.style.removeProperty("--ov-anchor-x");
      el.style.removeProperty("--ov-anchor-bottom");
      return;
    }
    const b = this.layersBtn.getBoundingClientRect();
    const rr = this.root.getBoundingClientRect();
    el.style.setProperty("--ov-anchor-x", `${(b.left + b.width / 2 - rr.left).toFixed(0)}px`);
    el.style.setProperty("--ov-anchor-bottom", `${(rr.bottom - b.top + 12).toFixed(0)}px`);
  }

  /** Show the full details of the selection (Focus mode's callout). */
  openDetails() {
    if (this.selection) this.inspector.show(this.selection);
  }

  syncViewToggle() {
    if (!this.viewToggle) return;
    const sky = this.viewer.state.view.mode === "sky";
    const label = sky ? "Back to 3D view" : "Sky view from the Sun";
    this.viewToggle.replaceChildren(icon(sky ? "cube" : "globe"), h("span", { class: "ov-ibtn-label" }, sky ? "3D" : "Sky"));
    this.viewToggle.setAttribute("aria-label", label);
    this.viewToggle.dataset.tip = `${label}  V`;
    this.viewToggle.hidden = !this.hasSky;
  }

  /** Everything the Focus bar does not show directly. */
  moreMenu(anchor) {
    const v = this.viewer;
    const lasso = this.plugins.find((p) => p.name === "lasso");
    const hasTime = !this.dock.el.hidden;
    const keyboard = !window.matchMedia?.("(any-hover: none)").matches;
    this.menu(anchor, [
      { label: this.mode === "focus" ? "Detailed mode" : "Focus mode", icon: "sidebar", shortcut: "U", run: () => this.toggleMode() },
      { label: "Search", icon: "search", shortcut: "/", run: () => this.palette.show() },
      ...(!this.statesBtn.hidden ? [{ label: "Views & story", icon: "bookmark", shortcut: "Y", run: () => this.plugins.find((p) => p.name === "states")?.toggle() }] : []),
      { label: "Reset view", icon: "home", shortcut: "Home", run: () => v.resetView() },
      ...(this.canQuickLook ? this.arMenuItems() : []),
      ...(hasTime ? [{ label: `Playback speed · ${v.timeline.speed}×`, icon: "play", run: () => this.dock.speedMenu?.(anchor) }] : []),
      "-",
      { label: "Screenshot or video…", icon: "camera", shortcut: "I", run: () => this.captureMenu(anchor) },
      { label: "Share & export…", icon: "share", run: () => this.shareMenu(anchor) },
      { label: "Display settings…", icon: "sliders", run: () => this.settingsMenu(anchor) },
      ...(this.widgets.available().length ? [{ label: "Widgets…", icon: "widgets", run: () => this.widgetsMenu(anchor) }] : []),
      ...this.plugins.flatMap((p) => p.moreItems?.() || []),
      "-",
      ...(lasso ? [{ label: "Lasso select", icon: "wand", shortcut: "L", run: () => lasso.arm(true) }] : []),
      { label: "Hide interface", icon: "eyeOff", shortcut: "Z", run: () => this.setZen(true) },
      ...(this.canFullscreen ? [{ label: "Fullscreen", icon: "expand", shortcut: "M", run: () => this.toggleFullscreen() }] : []),
      ...(keyboard ? [{ label: "Keyboard shortcuts", icon: "keyboard", shortcut: "?", run: () => this.showHelp() }] : []),
      { label: "About this figure", icon: "info", run: () => this.showAbout() },
    ], { align: "right", above: true });
  }

  setZen(on) {
    const was = this.root.dataset.zen === "true";
    this.root.dataset.zen = String(on);
    if (on) this.tooltips?.hide();
    if (on && !was && !this._zenHinted && this.root.dataset.presenting !== "true" && !window.matchMedia?.("(hover: none)").matches) {
      this._zenHinted = true;
      this.toast("Interface hidden · press Z to bring it back", { icon: icon("eyeOff"), ms: 2400 });
    }
    // Touch screens have no Z or Esc: a small pill brings the interface back.
    if (on && !this.zenExit) {
      this.zenExit = h("button", { class: "ov-zen-exit ov-glass", type: "button", onclick: () => this.setZen(false) }, icon("collapse"), "Show interface");
      this.ui.append(this.zenExit);
    }
    this.refit();
  }

  /** Re-measure the stage now (zen mode changes its frame). */
  refit() {
    this.viewer.renderer.resize();
    try { this.plugins.find((p) => p.name === "sky")?.sky?.aladin?.view?.fixLayoutDimensions?.(); } catch (_) { /* Aladin not ready */ }
  }

  // ------------------------------------------------------------- AR

  /** Apple's AR Quick Look (iPhone and iPad Safari) can open USDZ links. */
  get canQuickLook() {
    try { return document.createElement("a").relList.supports("ar"); } catch (_) { return false; }
  }

  /** Whether the figure has a timeline (so AR can be a time-lapse or a still). */
  get arHasTime() {
    return this.viewer.timeline.count > 1 && this.manifest.time?.enabled !== false;
  }

  /**
   * Build the AR model of what is on screen: a time-lapse through the whole
   * timeline, or (`moment`) a still of the current time.
   */
  async arModel({ moment = false } = {}) {
    const key = moment ? `moment:${this.viewer.timeline.frame}` : "timelapse";
    if (this._arBuilding?.key === key) return this._arBuilding.promise;
    const promise = buildArModel(this.viewer, this, { moment }).finally(() => { if (this._arBuilding?.promise === promise) this._arBuilding = null; });
    this._arBuilding = { key, promise };
    return promise;
  }

  /**
   * Open the figure in AR Quick Look. The time-lapse opens a model shipped
   * next to the page when there is one (figure option `ar_model`); anything
   * else is built into a USDZ on the spot.
   */
  async viewInAr({ moment = false } = {}) {
    const hosted = moment ? null : this.manifest.viewer?.arModel;
    let href = null, objectUrl = null;
    if (hosted) {
      try { href = new URL(hosted, location.href).href; } catch (_) { href = null; }
    }
    if (!href) {
      this.toast(moment ? `Preparing AR at ${this.formatTime(this.viewer.timeline.time)}…` : "Preparing the AR model…", { icon: icon("ar"), ms: 2200 });
      try {
        const { blob } = await this.arModel({ moment });
        objectUrl = href = URL.createObjectURL(blob);
      } catch (err) {
        console.error(err);
        this.toast(err.message || "The AR model could not be built");
        return;
      }
    }
    // Quick Look opens links with rel="ar" that wrap an image.
    const a = h("a", { rel: "ar", href: `${href}#allowsContentScaling=1`, style: { display: "none" } }, h("img", { alt: "" }));
    if (objectUrl) a.download = `${slug(this.manifest.title || "oviz-figure")}.usdz`;
    this.ui.append(a);
    a.click();
    a.remove();
    if (objectUrl) setTimeout(() => URL.revokeObjectURL(objectUrl), 120000);
  }

  /** AR entries: the time-lapse, and (with a timeline) a still of this moment. */
  arMenuItems() {
    return [
      { label: this.arHasTime ? "View in AR (time-lapse)" : "View in AR", icon: "ar", run: () => this.viewInAr() },
      ...(this.arHasTime ? [{ label: `View this moment in AR · ${Math.abs(this.viewer.timeline.time) < 1e-6 ? "today" : this.formatTime(this.viewer.timeline.time)}`, icon: "ar", run: () => this.viewInAr({ moment: true }) }] : []),
    ];
  }

  /** Download the AR model as a USDZ (open it on an iPhone or iPad for AR). */
  async saveArModel({ moment = false } = {}) {
    this.toast("Building the 3D model…", { icon: icon("ar"), ms: 1400 });
    try {
      const { blob, summary } = await this.arModel({ moment });
      const when = moment ? `-${this.formatTime(this.viewer.timeline.time).replace(/[^0-9a-z.−-]+/gi, "").replace("−", "minus")}` : "";
      downloadBlob(blob, `${slug(this.manifest.title || "oviz-figure")}${when}.usdz`);
      this.toast(`Saved USDZ · ${summary.objects.toLocaleString()} objects${summary.dustSlices ? `, ${summary.dustSlices} dust slices` : ""}${summary.seconds ? `, ${Math.round(summary.seconds)} s time-lapse` : ""}`, { icon: icon("ar") });
    } catch (err) {
      console.error(err);
      this.toast(err.message || "The 3D model could not be built");
    }
  }

  /** Whether this page may go fullscreen (not iPhone Safari, not a locked iframe). */
  get canFullscreen() {
    const d = document;
    return !!(d.fullscreenEnabled || d.webkitFullscreenEnabled) && !!(this.root.requestFullscreen || this.root.webkitRequestFullscreen);
  }

  toggleFullscreen() {
    const d = document;
    if (d.fullscreenElement || d.webkitFullscreenElement) {
      (d.exitFullscreen || d.webkitExitFullscreen)?.call(d);
      return;
    }
    const el = this.root;
    const request = el.requestFullscreen || el.webkitRequestFullscreen;
    const unavailable = () => this.toast("Fullscreen is not available here");
    if (!this.canFullscreen || !request) { unavailable(); return; }
    try { Promise.resolve(request.call(el)).catch(unavailable); } catch (_) { unavailable(); }
  }

  applyTheme(theme, { remember = true } = {}) {
    this.theme = theme === "light" ? "light" : "dark";
    document.documentElement.dataset.ovizTheme = this.theme;
    this.viewer.emit("theme", { theme: this.theme });
    const v = this.viewer;
    v.renderer.clearColor = [0, 0, 0, 0];
    v.renderer.invalidate();
    // Remember an explicit choice for this figure only.
    if (remember) localStorageSet(this.themeKey, this.theme);
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
      // Sky view: where the view is centred, in Galactic coordinates.
      const cam = this.viewer.renderer.camera;
      const f = cam.forward;
      const l = ((Math.atan2(f[1], f[0]) * 180) / Math.PI + 360) % 360;
      const b = (Math.asin(clamp(f[2], -1, 1)) * 180) / Math.PI;
      this.scaleLabel.textContent = formatAngle(sb.value);
      this.scaleSub.textContent = `l ${l.toFixed(1)}° b ${b >= 0 ? "+" : "−"}${Math.abs(b).toFixed(1)}° · ${formatAngle(cam.hfov)} field`;
    } else {
      this.scaleLabel.textContent = formatDistance(sb.value);
      const v = this.viewer;
      const eye = formatDistance(v.pose.distance);
      // Anchored, the eye distance is measured from the anchor.
      const name = v.anchored ? anchorLabel(v.anchor, this.anchorName(v.anchor)) : "";
      this.scaleSub.textContent = name ? `Eye ${eye} from ${name}` : `Eye ${eye} away`;
      this.scaleSub.title = name ? `Camera anchored to ${name}: it orbits, zooms toward and moves with it through time` : "";
    }
  }

  renderColorbars() {
    const v = this.viewer;
    const bars = [];
    for (const trace of v.traces) {
      const st = v.state.traces[trace.key] || {};
      const cb = trace.colorBy;
      if (!cb || !st.visible) continue;
      if ((st.colorMode || cb.defaultMode) !== "by_value") continue;
      const name = st.colormap || cb.colormap;
      if (!this.luts.get(name)) continue;
      bars.push({ dot: trace.legendColor, title: `${trace.name} · ${cb.label}`, name, lo: st.cmin ?? cb.cmin, hi: st.cmax ?? cb.cmax });
    }
    // Volumes that carry a colour scale, while they are drawn.
    const drawn = new Set((v._volumeParams(v.timeline.frame, v.timeline.time).volumeDraws || []).filter((d) => d.fade > 0.5).map((d) => d.key));
    const seen = new Set();
    for (const spec of v.volumeSpecs || []) {
      const cb = spec.colorbar;
      if (!cb || !drawn.has(spec.key) || seen.has(spec.stateKey)) continue;
      seen.add(spec.stateKey);
      const vs = v.state.volumes[spec.stateKey] || {};
      const name = vs.colormap || spec.colormaps[0];
      if (!this.luts.get(name)) continue;
      let lo, hi;
      if (Array.isArray(cb.range)) [lo, hi] = cb.range;
      else {
        // A rate scale: per unit time, so a longer window means a lower rate.
        const k = (Number(cb.scale) || 1) * (cb.perWindow && spec.procedural ? spec.procedural.windowMyr / (Number(vs.kdeTimeWindowMyr) || spec.procedural.windowMyr) : 1);
        lo = Math.max(0, Number(vs.vmin) || 0) * k;
        hi = Math.max(lo, Number(vs.vmax) || 0) * k;
      }
      bars.push({ dot: spec.legendColor, title: `${cb.label}${cb.unit ? ` (${cb.unit})` : ""}`, name, lo, hi });
    }
    const sig = JSON.stringify(bars.slice(0, 3));
    if (sig === this._colorbarSig) return;
    this._colorbarSig = sig;
    clear(this.colorbars);
    for (const b of bars.slice(0, 3)) {
      this.colorbars.append(h("div", { class: "ov-colorbar ov-glass" },
        h("div", { class: "ov-colorbar-title" }, h("span", { class: "ov-colorbar-dot", style: { background: b.dot } }), b.title),
        h("span", { class: "ov-colorbar-min" }, fmtNum(b.lo)),
        h("div", { class: "ov-colorbar-ramp", style: { background: lutGradient(this.luts.get(b.name)) } }),
        h("span", { class: "ov-colorbar-max" }, fmtNum(b.hi))));
    }
  }

  // ------------------------------------------------------------- menus

  menu(anchor, items, { title, align = "right", above = false, cls = "" } = {}) {
    // Pressing the button of the open menu closes it.
    if (this._menu?.anchor === anchor) { this.closeMenu(); return null; }
    this.closeMenu();
    this.tooltips?.hide();
    // Menus from the bar would open on top of a floating layers panel.
    if (this.layersOpen && this.layersFloating && !this.layers.el.contains(anchor)) this.setLayersOpen(false);
    const pop = h("div", { class: `ov-pop ov-glass ${cls}`.trim(), role: "menu" });
    if (title) pop.append(h("div", { class: "ov-pop-title" }, title));
    for (const it of items) {
      if (it === "-") { pop.append(h("div", { class: "ov-pop-sep" })); continue; }
      if (it.body) { pop.append(it.body); continue; }
      const b = h("button", { class: "ov-pop-item", type: "button", role: "menuitem" },
        it.icon ? icon(it.icon) : it.checked !== undefined ? h("span", { class: "ov-icon" }, it.checked ? "✓" : "") : null,
        h("span", { class: "ov-grow" }, it.label),
        it.shortcut ? kbd(it.shortcut) : null);
      // Keyboard activation (detail 0) returns focus to the menu's button.
      b.addEventListener("click", (e) => { if (!it.keepOpen) this.closeMenu({ refocus: e.detail === 0 }); it.run?.(); });
      pop.append(b);
    }
    this.ui.append(pop);
    const r = anchor.getBoundingClientRect();
    const rr = this.root.getBoundingClientRect();
    const w = pop.offsetWidth;
    let ht = pop.offsetHeight;
    let x = align === "center" ? r.left + r.width / 2 - w / 2 : align === "left" ? r.left : r.right - w;
    x = clamp(x - rr.left, 8, rr.width - w - 8);
    // Open on the preferred side if it fits, else on the roomier side, and
    // scroll (rather than run off the screen) when neither side fits.
    const roomAbove = r.top - rr.top - 16, roomBelow = rr.bottom - r.bottom - 16;
    const fitsPreferred = ht <= (above ? roomAbove : roomBelow);
    const placeAbove = fitsPreferred ? above : roomAbove > roomBelow;
    const room = Math.max(120, placeAbove ? roomAbove : roomBelow);
    if (ht > room) { pop.style.maxHeight = `${room}px`; ht = pop.offsetHeight; }
    let y = placeAbove ? r.top - rr.top - ht - 8 : r.bottom - rr.top + 8;
    y = clamp(y, 8, Math.max(8, rr.height - ht - 8));
    pop.style.left = `${x}px`;
    pop.style.top = `${y}px`;
    anchor.setAttribute("aria-expanded", "true");
    const onDown = (e) => {
      if (!pop.contains(e.target) && e.target !== anchor && !anchor.contains(e.target)) this.closeMenu();
    };
    setTimeout(() => this.root.addEventListener("pointerdown", onDown, true), 0);
    this._menu = { pop, anchor, onDown };
    pop.querySelector("button, input, select")?.focus({ preventScroll: true });
    // The menu grows out of its button while the scene pulls back.
    reveal(pop, r, { items: ".ov-pop-title, .ov-pop-item, .ov-pop-sep, .ov-pop-body > *" });
    this.setOverlay("menu", true);
    return pop;
  }

  closeMenu({ refocus = false } = {}) {
    const m = this._menu;
    if (!m) return false;
    this._menu = null;
    m.anchor.setAttribute("aria-expanded", "false");
    // Keyboard users land back on the button that opened the menu.
    if (refocus && m.pop.contains(document.activeElement) && m.anchor.isConnected) m.anchor.focus({ preventScroll: true });
    this.root.removeEventListener("pointerdown", m.onDown, true);
    this.setOverlay("menu", false);
    conceal(m.pop, m.anchor.isConnected ? m.anchor.getBoundingClientRect() : null).then(() => m.pop.remove());
    return true;
  }

  /** Track open overlays (data-overlay on the root while any is open). */
  setOverlay(kind, on) {
    if (on) this.overlays.add(kind);
    else this.overlays.delete(kind);
    this.idle?.wake();
    if (this.overlays.size) this.root.dataset.overlay = [...this.overlays].join(" ");
    else delete this.root.dataset.overlay;
  }

  settingsMenu(anchor) {
    const v = this.viewer;
    const g = v.state.global;
    const body = h("div", { class: "ov-pop-body" },
      miniSeg({ label: "Mode", value: this.mode, options: MODES.map((m) => ({ value: m.id, label: m.name })), onChange: (m) => { this.closeMenu(); this.setMode(m); } }),
      miniSeg({ label: "Theme", value: this.theme, options: [{ value: "dark", label: "Dark" }, { value: "light", label: "Light" }], onChange: (t) => this.applyTheme(t) }),
      this.anchorField(),
      slider({ label: "Point size", min: 0.25, max: 4, value: g.pointSize, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setGlobal({ pointSize: x }) }),
      slider({ label: "Point brightness", min: 0, max: 2, value: g.pointOpacity, format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setGlobal({ pointOpacity: x }) }),
      slider({ label: "Star glow", min: 0, max: 4, value: g.glow, format: (x) => (x <= 0.02 ? "Markers" : x.toFixed(2)), onInput: (x) => v.setGlobal({ glow: x }) }),
      slider({ label: "Birth fade", min: 0, max: 30, step: 0.5, value: g.fadeTime, format: (x) => `${x.toFixed(1)} Myr`, onInput: (x) => v.setGlobal({ fadeTime: x }) }),
      toggle({ label: "Fade by opacity", hint: "Otherwise newborn clusters grow in size", checked: g.fadeByOpacity, onChange: (x) => v.setGlobal({ fadeByOpacity: x }) }),
      toggle({ label: "Fade out after birth", checked: g.fadeInOut, onChange: (x) => v.setGlobal({ fadeInOut: x }) }),
      toggle({ label: "Size by member count", checked: g.sizeByStars, onChange: (x) => v.setGlobal({ sizeByStars: x }) }),
      v.timeline.count > 1 ? slider({ label: "Motion trails", min: 0, max: this.maxTrailMyr(), step: this.maxTrailMyr() / 60, value: g.trails || 0, format: (x) => (x <= 0 ? "Off" : `${x.toFixed(x >= 10 ? 0 : 1)} Myr`), onInput: (x) => v.setGlobal({ trails: x }) }) : null,
      (this._orbitToggle = toggle({ label: "Auto-orbit camera", hint: "Any drag, scroll or key stops it", checked: !!v.controls.autoOrbit, onChange: (x) => this.setAutoOrbit(x) })),
      slider({ label: "Orbit speed", min: 0.1, max: 4, value: g.autoOrbitSpeed ?? 1, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => { g.autoOrbitSpeed = x; v.controls.autoOrbit = autoOrbitRate(g); } }),
      miniSeg({ label: "Orbit direction", value: g.autoOrbitDirection < 0 ? -1 : 1, options: [{ value: 1, label: "Clockwise" }, { value: -1, label: "Anticlockwise" }], onChange: (d) => { g.autoOrbitDirection = d; v.controls.autoOrbit = autoOrbitRate(g); } }),
      slider({ label: "Field of view", min: 1, max: 120, value: v.pose.fov, scale: "log", format: (x) => `${x < 10 ? x.toFixed(1) : Math.round(x)}°`, onInput: (x) => { v.pose.fov = Math.min(x, v.controls.maxFovFor()); v.renderer.invalidate(); } }),
      slider({ label: "Scroll & key speed", min: 0.2, max: 4, value: g.navSpeed ?? 1, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => { g.navSpeed = x; v.controls.speed = x; } }),
      toggle({ label: "Performance overlay", checked: !!this.perfEl, onChange: (x) => this.setPerfOverlay(x) }),
      ...this.plugins.flatMap((p) => p.settingsRows?.() || []),
      h("button", { class: "ov-btn ov-btn--sm ov-btn--block", type: "button", onclick: () => { this.closeMenu(); this.resetDisplay(); } }, "Reset display settings"),
    );
    this.menu(anchor, [{ body }], { title: "Display" });
  }

  /** Back to the figure's own display settings (classic "Reset controls"). */
  resetDisplay() {
    const v = this.viewer;
    const g0 = initialViewerState(v.manifest).global;
    const keep = { grid: v.state.global.grid, gridOpacity: v.state.global.gridOpacity, labels: v.state.global.labels };
    v.setGlobal({ ...g0, ...keep });
    v.controls.speed = g0.navSpeed ?? 1;
    this.setAutoOrbit(!!g0.autoOrbit);
    const fov = v.homePose?.fov;
    if (Number.isFinite(fov) && v.state.view.mode !== "sky") { v.pose.fov = fov; v.renderer.invalidate(); }
    this.toast("Display settings reset", { icon: icon("sliders"), ms: 1200 });
  }

  /** The Widgets menu: every tool panel this figure offers, ticked when open. */
  widgetsMenu(anchor) {
    this.menu(anchor, this.widgets.menuItems(), { title: "Widgets" });
  }

  /** Plugins add a button to the toolbar (before the separator). */
  addToolbarButton(btn) {
    this.toolbar.insertBefore(btn, this.toolbar.querySelector(".ov-toolbar-sep"));
    return btn;
  }

  // ------------------------------------------------------------- modes

  /** Switch between Focus and Detailed mode. */
  setMode(id) {
    const prev = this.mode;
    this.mode = applyMode(id);
    if (prev === this.mode) return;
    this.viewer.emit("mode", { mode: this.mode });
    this.applyLayout();
    if (!this.narrow) this.setLayersOpen(this.modeConfig.layersOpen);
    if (this.selection) {
      this.inspector.show(this.modeConfig.selection === "panel" ? this.selection : null);
      this.reticle.onSelect();
    }
    this.idle.wake();
    const settle = () => {
      this.syncViewSeg();
      this.reticle.update();
      this.anchorLayers();
      this.inspector.track?.draw();
      this.updateScale();
    };
    requestAnimationFrame(settle);
    setTimeout(settle, 320);
  }

  toggleMode() {
    this.setMode(this.mode === "focus" ? "detailed" : "focus");
  }

  /** Longest useful trail: a third of the timeline. */
  maxTrailMyr() {
    const tl = this.viewer.timeline;
    return Math.max(0.5, Math.abs(tl.times[tl.count - 1] - tl.times[0]) / 3);
  }

  /** Toggle motion trails at a length that reads well for this timeline. */
  toggleTrails() {
    const v = this.viewer;
    if (v.timeline.count < 2) return;
    const on = !(v.state.global.trails > 0);
    const len = on ? +(this.maxTrailMyr() * 0.3).toPrecision(2) : 0;
    v.setGlobal({ trails: len });
    this.toast(on ? `Motion trails · ${len} Myr` : "Motion trails off", { ms: 1200 });
  }

  setAutoOrbit(on) {
    const g = this.viewer.state.global;
    g.autoOrbit = !!on;
    this.viewer.controls.autoOrbit = autoOrbitRate(g);
    this._orbitToggle?.set(!!on);
    this.viewer.renderer.invalidate();
  }

  captureMenu(anchor) {
    const items = [
      { label: "Save PNG", icon: "camera", shortcut: "I", run: () => this.screenshot(1) },
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
    const story = this.plugins.find((p) => p.name === "states");
    const views = story?.project.items.length || 0;
    const items = [
      { label: views ? "Copy link to this view (with the saved views)" : "Copy link to this view", icon: "share", run: () => this.copyViewLink() },
      ...(views ? [{ label: "Copy presentation link", icon: "present", run: () => this.copyPresentationLink() }] : []),
      ...(this.canQuickLook ? this.arMenuItems() : []),
      { label: this.arHasTime ? "Save AR time-lapse (USDZ)" : "Save 3D model for AR (USDZ)", icon: "ar", run: () => this.saveArModel() },
      ...(this.arHasTime ? [{ label: "Save AR model of this moment (USDZ)", icon: "ar", run: () => this.saveArModel({ moment: true }) }] : []),
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
        background: sky?.captureBackground ? (ctx, w, h) => sky.captureBackground(ctx, w, h) : "#000",
      });
      if (clipboard && navigator.clipboard?.write) {
        await navigator.clipboard.write([new ClipboardItem({ "image/png": blob })]);
        this.toast("Image copied", { icon: icon("copy") });
      } else {
        downloadBlob(blob, `${slug(this.manifest.title || "oviz-figure")}-${stamp()}.png`);
        const size = v.renderer.lastCapture;
        this.toast(size ? `Saved PNG · ${size.width} × ${size.height}` : "Saved PNG", { icon: icon("camera") });
      }
    } catch (err) {
      console.error(err);
      this.toast("Screenshot failed");
    }
  }

  copyViewLink() {
    const n = this.plugins.find((p) => p.name === "states")?.project.items.length || 0;
    copyTextLater(this.viewLinkUrl())
      .then(() => this.toast(n ? `Link copied · it opens this view and carries the ${n} saved view${n === 1 ? "" : "s"}` : "Link to this view copied", { icon: icon("share") }))
      .catch(() => this.toast("Could not copy the link"));
  }

  /** A link that opens presenting the saved views from the first one. */
  copyPresentationLink() {
    const n = this.plugins.find((p) => p.name === "states")?.project.items.length || 0;
    if (!n) { this.toast("Save a view with N first, then copy a presentation link"); return; }
    copyTextLater(this.viewLinkUrl({ present: true }))
      .then(() => this.toast(`Presentation link copied · it opens presenting ${n} view${n === 1 ? "" : "s"}`, { icon: icon("present") }))
      .catch(() => this.toast("Could not copy the link"));
  }

  /**
   * A link to this view: this page with the view in its hash, plus the saved
   * views (`views: false` leaves them out); `present: true` makes it open
   * presenting them from the first.
   */
  async viewLinkUrl(opts = {}) {
    const url = new URL(location.href);
    url.hash = await this.viewLinkHash(opts);
    return url.toString();
  }

  /**
   * The view link's hash: the camera and time in readable form, plus the
   * whole State (what differs from how the figure opens), so every property
   * of the view (layer styles, volumes, the Sky background, filters, the
   * lasso, notes, widgets and panels) arrives as it was sent.
   */
  async viewLinkHash({ views = true, present = false } = {}) {
    const story = this.plugins.find((p) => p.name === "states");
    let state = "", shared = "";
    if (story?.homeState) {
      try { state = await encodeStatePart(captureState(this.viewer, this.sky), story.homeState); } catch (_) { state = ""; }
      if (views && story.project.items.length) {
        try { shared = await encodeViewsPart(story.project.items, story.homeState, story.project.defaultTransition); } catch (_) { shared = ""; }
      }
    }
    const extra = state ? { state } : this.viewLinkExtras();
    if (shared) {
      extra.views = shared;
      if (present) extra.present = 1;
    }
    return encodeViewHash(this.viewer, extra);
  }

  /**
   * What a view link keeps beyond the camera and time: the layers shown
   * (volumes too), the lasso selection with the outline that clips the
   * dust, and the selected object. Traces go by their index in the figure.
   */
  viewLinkExtras() {
    const v = this.viewer;
    const index = new Map(v.traces.map((t, i) => [t.key, i]));
    const extra = { layers: this.layers.toggleItems("all").map((it) => this.layers.itemVisible(it)) };
    const lasso = this.plugins.find((p) => p.name === "lasso");
    if (lasso?.selection) {
      extra.lasso = {
        selection: [...lasso.selection].filter(([key]) => index.has(key)).map(([key, idx]) => [index.get(key), Array.from(idx)]),
        isolate: lasso.isolate,
        filterOff: !!lasso.filterOff,
        mask: lasso.mask,
      };
    }
    const sel = this.selection;
    if (sel && sel.kind !== "member" && index.has(sel.trace) && Number.isInteger(sel.index)) extra.object = { trace: index.get(sel.trace), index: sel.index };
    return extra;
  }

  // ------------------------------------------------------------- palette sources

  commands() {
    const v = this.viewer;
    const tl = v.timeline;
    const list = [
      { title: tl.playing ? "Pause" : "Play time", icon: tl.playing ? "pause" : "play", shortcut: "Space", pinned: true, run: () => this.dock.togglePlay(tl.playing ? 0 : 1) },
      ...(tl.count > 1 ? [{ title: "Play backward", icon: "stepBack", shortcut: "⇧ Space", keywords: "reverse rewind time past", run: () => tl.play(-1) }] : []),
      { title: "Reset view", icon: "home", shortcut: "Home", pinned: true, run: () => v.resetView() },
      { title: "Go to present day (t = 0)", icon: "target", keywords: "now zero", run: () => tl.setTime(0) },
      { title: this.layersOpen ? "Hide layers panel" : "Show layers panel", icon: "layers", shortcut: "⇧ L", run: () => this.setLayersOpen(!this.layersOpen) },
      { title: v.state.global.grid ? "Hide Galactic grid" : "Show Galactic grid", icon: "grid", shortcut: "G", run: () => v.setGlobal({ grid: !v.state.global.grid }) },
      { title: "Save screenshot", icon: "camera", shortcut: "I", run: () => this.screenshot(1) },
      { title: "Copy link to this view", icon: "share", run: () => this.copyViewLink() },
      { title: this.theme === "dark" ? "Switch to light theme" : "Switch to dark theme", icon: this.theme === "dark" ? "sun" : "moon", run: () => this.applyTheme(this.theme === "dark" ? "light" : "dark") },
      { title: v.controls.autoOrbit ? "Stop auto-orbit" : "Auto-orbit camera", icon: "orbit", shortcut: "O", run: () => this.setAutoOrbit(!v.controls.autoOrbit) },
      ...["lsr", "sun", "free"].filter((k) => k !== v.anchor.kind && (k === "free" || v.anchorAvailable({ kind: k }))).map((k) => ({
        title: k === "free" ? "Free the camera (no anchor)" : `Anchor camera to ${anchorLabel({ kind: k })}`,
        icon: k === "free" ? "orbit" : "target",
        keywords: "camera anchor lock follow centre center pivot frame reference lsr sun local standard of rest",
        run: () => this.setCameraAnchor(k),
      })),
      ...(tl.count > 1 ? [{ title: v.state.global.trails > 0 ? "Hide motion trails" : "Show motion trails", icon: "trail", shortcut: "J", keywords: "streaks motion history comet", run: () => this.toggleTrails() }] : []),
      { title: "Hide interface (zen)", icon: "expand", shortcut: "Z", run: () => this.setZen(this.root.dataset.zen !== "true") },
      { title: this.mode === "focus" ? "Detailed mode" : "Focus mode", sub: this.mode === "focus" ? MODES[1].desc : MODES[0].desc, icon: "sidebar", shortcut: "U", keywords: "mode layout panels simple minimal focus detailed", run: () => this.toggleMode() },
      { title: "Keyboard shortcuts", icon: "keyboard", shortcut: "?", run: () => this.showHelp() },
      ...(this.canQuickLook ? this.arMenuItems().map((x) => ({ title: x.label, icon: "ar", keywords: "augmented reality iphone quick look table 3d time-lapse", run: x.run })) : []),
      { title: this.arHasTime ? "Save AR time-lapse (USDZ)" : "Save 3D model for AR (USDZ)", icon: "ar", keywords: "augmented reality iphone quick look usdz export 3d", run: () => this.saveArModel() },
      ...(this.arHasTime ? [{ title: "Save AR model of this moment (USDZ)", icon: "ar", keywords: "augmented reality iphone quick look usdz export 3d time", run: () => this.saveArModel({ moment: true }) }] : []),
      ...(this.canFullscreen ? [{ title: "Fullscreen", icon: "expand", shortcut: "M", run: () => this.toggleFullscreen() }] : []),
    ];
    if (this.manifest.sky?.enabled) {
      list.splice(2, 0, { title: v.state.view.mode === "sky" ? "Switch to Galactic 3D" : "Switch to Sky view", icon: v.state.view.mode === "sky" ? "cube" : "globe", shortcut: "V", pinned: true, run: () => this.setViewMode(v.state.view.mode === "sky" ? "3d" : "sky") });
    }
    list.push(...this.widgets.commands());
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
    const all = this.layers.toggleItems("all");
    if (all.length > 1) {
      const anyOn = all.some((it) => this.layers.itemVisible(it));
      items.push({ title: anyOn ? "Hide all layers" : "Show all layers", sub: this.manifest.volumes?.length ? "Data and volumes; showing again restores what was on" : "", icon: anyOn ? "eyeOff" : "eye", shortcut: "T", keywords: "toggle everything traces volumes", run: () => this.layers.toggleAll("all") });
    }
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
    reveal(el, null, { spring: SPRINGS.soft, items: ".ov-sheet-body > * > *" });
    this.setOverlay("sheet", true);
    return el;
  }

  closeSheet() {
    if (!this._sheet) return false;
    const { scrim, el } = this._sheet;
    this._sheet = null;
    this.setOverlay("sheet", false);
    scrim.style.transition = "opacity 200ms";
    scrim.style.opacity = "0";
    conceal(el, null).then(() => { el.remove(); scrim.remove(); });
    this.focusCanvas();
    return true;
  }

  showHelp() {
    const groups = [
      ["Camera", [["Orbit", "Drag"], ["Pan", "Right-drag / ⌃ Drag"], ["Zoom about the centre", "Scroll / Middle-drag"], ["Orbit / tilt", "W A S D"], ["Fly (4× faster)", "⇧ W A S D"], ["Zoom out / in", "Q E"], ["Move up / down", "R F"], ["Reset view and re-anchor", "Home"], ["Auto-orbit", "O"], ["3D ⇄ Sky", "V"]]],
      ["Time", [["Play / pause", "Space"], ["Play backward", "⇧ Space"], ["Step frame", "← →"], ["Step 5 frames", "⇧ ← →"], ["Slower / faster", "< >"], ["Present day", "0"]]],
      ["Layers & display", [["Toggle layer 1–9", "1–9"], ["Solo layer", "⇧ 1–9"], ["Show / hide all (volumes too)", "T"], ["Point size − / +", "[ ]"], ["Motion trails", "J"], ["Galactic grid", "G"], ["Sky background", "B"], ["Layers panel", "⇧ L"]]],
      ["Select", [["Search anything", `${MOD} K`], ["Select object", "Click"], ["Follow a cluster (its frame)", "Double-click"], ["Back to the LSR frame", "Double-click elsewhere"], ["Measure separation", "⇧ Click"], ["Lasso", "⇧ Drag / L"], ["Lasso filter on / off", "C"], ["Undo selection", `${MOD} Z`], ["Clear selection", "Esc"]]],
      ["Views & presenting", [["Save current view", "N"], ["Present", "P"], ["Next / previous view", "→ ←"], ["Views panel", "Y"], ["Save figure", `${MOD} S`]]],
      ["Capture", [["Screenshot", "I"], ["Hide interface", "Z"], ["Fullscreen", "M"], ["Focus / detailed mode", "U"], ["This help", "?"]]],
      ...this.plugins.flatMap((p) => p.helpGroups?.() || []),
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
    // A mouse click should not leave a button focused (Space would then
    // re-press it instead of playing time); keyboard activation keeps focus.
    this.root.addEventListener("click", (e) => {
      const b = e.target.closest?.("button");
      if (b && e.detail > 0 && document.activeElement === b) b.blur();
    });
    // Held movement keys, applied every frame (classic Oviz semantics).
    this.keysDown = new Set();
    this.shiftDown = false;
    v.renderer.beforeRender.push((now) => this.applyKeyMotion(now));
    const clearKeys = () => { this.keysDown.clear(); this.shiftDown = false; v.renderer.release("keys"); };
    window.addEventListener("blur", clearKeys);
    window.addEventListener("keyup", (e) => {
      // macOS sends no keyup for keys released while ⌘ is held: releasing ⌘
      // lets go of every held movement key, so none stays stuck.
      if (e.key === "Meta") { clearKeys(); return; }
      if (e.key === "Shift") this.shiftDown = false;
      this.keysDown.delete(e.key.toLowerCase());
      if (!this.keysDown.size) v.renderer.release("keys");
    });
    window.addEventListener("keydown", (e) => {
      if (isEditable(e.target)) return;
      // Let focused sliders, buttons and switches handle their own keys.
      if (controlConsumesKey(e)) return;
      const mod = e.metaKey || e.ctrlKey;
      const k = e.key;
      const lower = k.length === 1 ? k.toLowerCase() : k;
      if (mod && lower === "k") { this.palette.toggle(); e.preventDefault(); return; }
      if (mod && !e.altKey && lower === "z" && !e.shiftKey) {
        if (this.undoSelection()) e.preventDefault();
        return;
      }
      if (mod) {
        if (e.repeat) return;
        for (const p of this.plugins) if (p.onKey?.(e)) { e.preventDefault(); return; }
        return;
      }
      if (e.altKey) return;
      // The palette handles keys in its field; if focus left the field,
      // Escape still closes it.
      if (this.palette.open) {
        if (k === "Escape") { this.palette.close(); e.preventDefault(); }
        return;
      }
      if (k === "Shift") { this.shiftDown = true; return; }
      if ("wasdqerf".includes(lower) && lower.length === 1) {
        // Movement keys are held, not pressed: see applyKeyMotion.
        this.shiftDown = e.shiftKey;
        if (!this.keysDown.has(lower)) {
          this.keysDown.add(lower);
          this._lastKeyMotion = 0;
          v._cancelTween();
          v.controls.stop();
          v.emit("camera-takeover", {});
          v.emit("user-camera", {});
          v.renderer.hold("keys");
        }
        e.preventDefault();
        return;
      }
      // Held keys repeat only where repeating makes sense (stepping, sizes);
      // toggles such as N, L, B, Y or P fire once per press.
      if (e.repeat && !["ArrowLeft", "ArrowRight", "[", "]", "{", "}", "<", ">", ",", "."].includes(k)) {
        e.preventDefault();
        return;
      }
      for (const p of this.plugins) if (p.onKey?.(e)) { e.preventDefault(); return; }
      let handled = true;
      switch (lower) {
        // ---- classic Oviz keys
        case " ": this.dock.togglePlay(e.shiftKey ? -1 : 1); break;
        case "ArrowLeft": tl.pause(); tl.step(e.shiftKey ? -5 : -1); break;
        case "ArrowRight": tl.pause(); tl.step(e.shiftKey ? 5 : 1); break;
        case "?": this.showHelp(); break;
        case "v": this.setViewMode(v.state.view.mode === "sky" ? "3d" : "sky"); break;
        case "g": v.setGlobal({ grid: !v.state.global.grid }); break;
        case "o": this.setAutoOrbit(!v.controls.autoOrbit); break;
        case "j": this.toggleTrails(); break;
        case "u": this.toggleMode(); break;
        case "z": this.setZen(this.root.dataset.zen !== "true"); break;
        case "t": this.layers.toggleAll("all"); break;
        case "[": case "{": case "]": case "}": {
          const dir = k === "]" || k === "}" ? 1 : -1;
          const step = e.shiftKey ? 0.25 : 0.1;
          v.setGlobal({ pointSize: Math.min(4, Math.max(0.1, +(v.state.global.pointSize + dir * step).toFixed(2))) });
          this.toast(`Point size ${v.state.global.pointSize.toFixed(2)}×`, { ms: 900 });
          break;
        }
        case "Escape":
          if (this.closeMenu({ refocus: true }) || this.closeSheet()) break;
          // An open layer editor closes first (classic), then the panel.
          if (this.layersOpen && this.layers.expanded) { this.layers.expand(this.layers.expanded); break; }
          if (this.layersOpen && this.layersFloating) { this.setLayersOpen(false); break; }
          if (this.root.dataset.story === "true") { this.plugins.find((p) => p.name === "states")?.toggle(false); break; }
          if (this.root.dataset.zen === "true") { this.setZen(false); break; }
          if (this.selection || this.plugins.find((p) => p.name === "lasso")?.selection) {
            // One undo step covers both the click and the lasso selection.
            this.pushSelectionUndo();
            this.select(null, { recordUndo: false });
            this.plugins.find((p) => p.name === "lasso")?.clear();
            break;
          }
          handled = false;
          break;
        // ---- additions (keys the classic viewer did not use)
        case "/": this.palette.show(); break;
        case "0": tl.setTime(0); break;
        case "<": case ",": this.dock.cycleSpeed(-1); break;
        case ">": case ".": this.dock.cycleSpeed(1); break;
        case "i": this.screenshot(1); break;
        case "m": this.toggleFullscreen(); break;
        case "Home": v.resetView(); break;
        case "l": if (e.shiftKey) this.setLayersOpen(!this.layersOpen); else handled = false; break;
        default:
          if (/^Digit[1-9]$/.test(e.code)) {
            this.layers.toggleItemByIndex(Number(e.code.slice(5)) - 1, { solo: e.shiftKey });
          } else handled = false;
      }
      if (handled) e.preventDefault();
    });
  }

  /**
   * Classic Oviz keyboard flight, applied per frame while keys are held:
   * A/D orbit, W/S tilt (W toward overhead), Shift+W/A/S/D fly the camera
   * and target, Q/E zoom out/in, R/F move up/down; Shift is 4× faster.
   * In Sky view W/A/S/D turn the gaze and Q/E change the field of view.
   */
  applyKeyMotion(now) {
    const keys = this.keysDown;
    if (!keys || !keys.size) { this._lastKeyMotion = 0; return; }
    const v = this.viewer;
    const pose = v.renderer.camera.pose;
    const dt = this._lastKeyMotion ? Math.min(Math.max((now - this._lastKeyMotion) / 1000, 0), 0.1) : 1 / 60;
    this._lastKeyMotion = now;
    const c = v.controls;
    const sky = v.state.view.mode === "sky";
    // Flying (Shift+W/A/S/D) and R/F move within the anchor's frame, like panning.
    keyMotion(pose, keys, {
      dt,
      fast: this.shiftDown,
      sky,
      maxSpan: v.world.maxSpan,
      minDistance: c.minDistance, maxDistance: c.maxDistance, minFov: c.minFov, maxFov: c.maxFovFor(sky ? "sky" : "galactic"),
      speed: c.speed,
    });
    v.renderer.markInteraction();
    v._onCameraChange();
  }

  // ------------------------------------------------------------- selection undo

  pushSelectionUndo() {
    const lasso = this.plugins.find((p) => p.name === "lasso");
    (this._undo ||= []).push({ selection: this.selection, lasso: lasso?.capture() ?? null });
    if (this._undo.length > 32) this._undo.shift();
  }

  undoSelection() {
    const entry = this._undo?.pop();
    if (!entry) return false;
    this.select(entry.selection, { recordUndo: false });
    this.plugins.find((p) => p.name === "lasso")?.restore(entry.lasso);
    this.toast("Selection restored", { ms: 1000 });
    return true;
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
