// "Views & story": saved States as a filmstrip, presentation mode, drafts
// and self-export.

import { h, icon, iconButton, clear, kbd, MOD } from "./dom.js";
import { slider, toggle, select } from "./controls.js";
import { captureState, importStates, exportStatesBlock, makeTransition, easing, uid, readTransition } from "../app/states.js";
import { cloneJson } from "../app/state.js";
import { buildExportHtml, saveHtml } from "../app/export.js";
import { readDraft, writeDraft, deleteDraft } from "../app/drafts.js";
import { clamp } from "../core/math.js";
import { clonePose } from "../engine/camera.js";

export class StoryPlugin {
  constructor() {
    this.name = "states";
    this.active = -1;
    this.running = null;
    this.presenting = false;
    this.dirty = false;
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.sky = ui.plugins.find((p) => p.name === "sky") || null;
    this.homeState = null;
    this.project = importStates(v.manifest, v.state);
    if (!this.project.projectId) this.project.projectId = figureId(v.manifest);
    this.readOnly = this.project.presentOnly;
    this.buildDrawer();
    this.buildPresenter();
    // The user grabbing the camera (drag, wheel, keys, a flight, 3D/Sky)
    // ends a running transition; the camera stays where it is.
    v.on("camera-takeover", () => this.takeover());
    // Keep drafts when the tab goes away, and warn before leaving only when
    // unsaved views are not safely in this browser.
    window.addEventListener("pagehide", () => this.flushDraft());
    window.addEventListener("beforeunload", (e) => {
      if (!this.dirty || this.readOnly || !this.project.items.length) return;
      if (this.project.autosave !== false && this.draftSaved) return;
      e.preventDefault();
      e.returnValue = "";
    });
  }

  afterBoot() {
    const v = this.viewer;
    // "Current view only" exports carry one hidden State applied at load.
    if (v.manifest.states?.apply_first_on_load && this.project.items.length) {
      const first = this.project.items[0];
      const t = makeTransition(v, this.sky, first.state);
      t.finish();
      this.afterArrive(t);
      this.project.items = [];
      // Home returns to the saved view, as the figure opened (classic).
      const home = v.state.view.mode === "sky" ? v.returnPose : v.pose;
      if (home) v.homePose = clonePose(home);
      if (v.state.view.mode !== "sky" && v.state.view.anchor) v.defaultAnchor = { ...v.state.view.anchor };
    }
    this.homeState = captureState(v, this.sky);
    this.render();
    this.restoreDraft().then(() => this.render());
    // Present-only exports open straight into the first view. (The legacy
    // `default_mode` flag only chose the drawer tab, so it does not.)
    if (this.readOnly && this.project.items.length) {
      setTimeout(() => this.present(true, { startAt: 0, instant: true }), 60);
    }
    this.ui.statesBtn.hidden = this.readOnly && !this.project.items.length;
  }

  // ------------------------------------------------------------ drawer UI

  buildDrawer() {
    const ui = this.ui;
    this.el = h("section", { class: "ov-story ov-glass ov-chrome", "aria-label": "Views and story", "data-open": "false" });
    this.count = h("span", { class: "ov-story-count" });
    this.addBtn = h("button", { class: "ov-btn ov-btn--primary ov-btn--sm", type: "button", "aria-label": "Save view", onclick: () => this.add() }, icon("plus"), h("span", { class: "ov-btn-label" }, "Save view"));
    this.presentBtn = h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-label": "Present", onclick: () => this.present(true) }, icon("present"), h("span", { class: "ov-btn-label" }, "Present"));
    this.homeBtn = iconButton("home", "Original view (the figure as it opened)", () => this.goHome());
    this.exportBtn = h("button", { class: "ov-btn ov-btn--sm ov-btn--ghost", type: "button", "aria-label": "Export", onclick: (e) => ui.menu(e.currentTarget, [...this.linkMenuItems(), ...(this.project.items.length ? ["-"] : []), ...this.exportMenuItems()], { title: "Share & export", above: true }) }, icon("download"), h("span", { class: "ov-btn-label" }, "Export"));
    this.settingsBtn = iconButton("sliders", "Transition timing", (e) => this.timingMenu(e.currentTarget));
    const head = h("div", { class: "ov-story-head" },
      h("span", { class: "ov-story-title" }, "Story"), this.count,
      h("span", { class: "ov-grow" }),
      this.readOnly ? null : this.addBtn, this.presentBtn, this.homeBtn, this.readOnly ? null : this.exportBtn, this.readOnly ? null : this.settingsBtn,
      iconButton("close", "Close", () => this.toggle(false), { shortcut: "Y" }));
    this.strip = h("div", { class: "ov-story-strip", role: "list" });
    this.status = h("div", { class: "ov-story-status" });
    this.el.append(head, this.strip, this.status);
    ui.ui.append(this.el);
  }

  toggle(open = this.el.dataset.open !== "true") {
    // The drawer and a floating layers panel share the space above the bar.
    if (open && this.ui.layersOpen && this.ui.layersFloating) this.ui.setLayersOpen(false);
    this.el.dataset.open = String(open);
    this.ui.statesBtn.setAttribute("aria-pressed", String(open));
    this.ui.root.dataset.story = String(open);
    if (open) this.render();
  }

  render() {
    const items = this.project.items;
    this.count.textContent = items.length ? `${items.length} view${items.length === 1 ? "" : "s"}` : "";
    this.presentBtn.disabled = !items.length;
    clear(this.strip);
    if (!items.length) {
      this.strip.append(h("div", { class: "ov-story-empty" },
        h("div", { class: "ov-story-empty-title" }, "Build a story from saved views"),
        h("div", null, "Frame a view — camera, time, layers, Sky — then save it with ", kbd("N"), ". Present them in order with smooth transitions, or export a present-only file to share.")));
      this.status.textContent = "";
      this.ui.dock.setMarkers?.([]);
      return;
    }
    items.forEach((it, i) => this.strip.append(this.card(it, i)));
    if (!this.readOnly) {
      this.strip.append(h("button", { class: "ov-story-add", type: "button", onclick: () => this.add(), "data-tip": "Save current view  N" }, icon("plus")));
    }
    this.status.textContent = this.dirty ? "Unsaved changes · autosaved in this browser" : "";
    this.ui.dock.setMarkers?.(items.map((it, index) => ({ time: this.viewer.timeline.frameToTime(it.state.time?.frame ?? 0), name: it.name, index })));
  }

  card(it, i) {
    const thumb = h("div", { class: "ov-card-thumb" });
    if (it.thumb) thumb.style.backgroundImage = `url("${it.thumb}")`;
    else thumb.append(h("span", { class: "ov-card-mode" }, it.state.view?.mode === "sky" ? "Sky" : "3D"));
    const name = h("div", { class: "ov-card-name", title: it.name }, it.name);
    const card = h("div", { class: "ov-card", role: "listitem", tabindex: "0", "data-active": String(i === this.active), draggable: this.readOnly ? null : "true" },
      thumb,
      h("div", { class: "ov-card-meta" }, h("span", { class: "ov-card-index" }, String(i + 1)), name),
    );
    if (it.camera === "keep") card.append(h("span", { class: "ov-card-flag", title: "Keeps the current camera" }, icon("target")));
    card.addEventListener("click", () => this.goTo(i));
    card.addEventListener("keydown", (e) => {
      // Only keys the card itself handles; the rest (P, Space, arrows…)
      // still reach the figure after a card was clicked.
      if (e.target !== card) return;
      if (e.key === "Enter") { this.goTo(i); e.preventDefault(); e.stopPropagation(); }
      else if ((e.key === "Delete" || e.key === "Backspace") && !this.readOnly) { this.remove(i); e.preventDefault(); e.stopPropagation(); }
    });
    if (!this.readOnly) {
      const more = iconButton("sliders", `Edit ${it.name}`, (e) => { e.stopPropagation(); this.cardMenu(e.currentTarget, i); }, { cls: "ov-card-more" });
      card.append(more);
      name.addEventListener("dblclick", (e) => { e.stopPropagation(); this.rename(i, name); });
      card.addEventListener("dragstart", (e) => { e.dataTransfer.setData("text/oviz-state", String(i)); card.dataset.dragging = "true"; });
      card.addEventListener("dragend", () => { card.dataset.dragging = "false"; });
      card.addEventListener("dragover", (e) => { e.preventDefault(); card.dataset.over = "true"; });
      card.addEventListener("dragleave", () => { card.dataset.over = "false"; });
      card.addEventListener("drop", (e) => {
        e.preventDefault();
        const raw = e.dataTransfer.getData("text/oviz-state");
        const from = raw === "" ? NaN : Number(raw);
        card.dataset.over = "false";
        if (Number.isInteger(from) && from >= 0 && from !== i) this.move(from, i);
      });
    }
    return card;
  }

  cardMenu(anchor, i) {
    const it = this.project.items[i];
    this.ui.menu(anchor, [
      { label: "Go to view", icon: "target", run: () => this.goTo(i) },
      { label: "Update to current view", icon: "camera", run: () => this.update(i) },
      { label: "Rename…", icon: "edit", run: () => this.rename(i) },
      { label: "Edit caption…", icon: "edit", run: () => this.editCaption(i) },
      { label: it.camera === "keep" ? "Camera: keep current" : "Camera: fly to saved", icon: "follow", run: () => { it.camera = it.camera === "keep" ? "follow" : "keep"; this.changed(); } },
      { label: "Transition…", icon: "play", run: () => this.editTransition(i) },
      { label: "Duplicate", icon: "copy", run: () => this.duplicate(i) },
      ...(i > 0 ? [{ label: "Move earlier", icon: "stepBack", run: () => this.move(i, i - 1) }] : []),
      ...(i < this.project.items.length - 1 ? [{ label: "Move later", icon: "stepFwd", run: () => this.move(i, i + 1) }] : []),
      "-",
      { label: "Delete", icon: "trash", run: () => this.remove(i) },
    ], { title: it.name, above: true, align: "left" });
  }

  /** How the story flies into this view: its own length and easing, or the story's. */
  editTransition(i) {
    const it = this.project.items[i];
    const def = this.project.defaultTransition;
    const own = !!it.transition;
    const cur = { duration_ms: it.transition?.duration_ms ?? def.duration_ms, easing: it.transition?.easing ?? def.easing };
    const apply = () => { it.transition = { ...cur }; this.changed(); };
    const body = h("div", { class: "ov-dialog-grid" },
      h("div", { class: "ov-note" }, own ? "This view has its own transition." : "This view uses the story's transition until you change it here."),
      slider({ label: "Length", min: 0, max: 20, step: 0.1, value: cur.duration_ms / 1000, format: (x) => `${x.toFixed(1)} s`, onChange: (x) => { cur.duration_ms = Math.round(x * 1000); apply(); } }),
      select({ label: "Easing", value: cur.easing, options: [
        { value: "easeInOutCubic", label: "Smooth (ease in and out)" }, { value: "easeInOutQuad", label: "Gentle" },
        { value: "easeOutCubic", label: "Ease out" }, { value: "linear", label: "Linear" },
      ], onChange: (x) => { cur.easing = x; apply(); } }),
      h("div", null, h("button", { class: "ov-btn ov-btn--ghost", type: "button", onclick: () => { it.transition = null; this.changed(); this.ui.closeSheet(); } }, "Use the story's transition")));
    this.ui.sheet(`Transition · ${it.name}`, body);
  }

  timingMenu(anchor) {
    const dt = this.project.defaultTransition;
    const body = h("div", { class: "ov-pop-body" },
      slider({ label: "Transition length", min: 0, max: 20, step: 0.1, value: dt.duration_ms / 1000, format: (x) => `${x.toFixed(1)} s`, onChange: (x) => { dt.duration_ms = Math.round(x * 1000); this.changed(); } }),
      toggle({ label: "Autosave drafts in this browser", checked: this.project.autosave !== false, onChange: (x) => { this.project.autosave = x; this.changed(); } }),
    );
    this.ui.menu(anchor, [{ body }], { title: "Story settings", above: true });
  }

  // ------------------------------------------------------------ editing

  /** Save the current view. The N key saves quietly (classic); the button opens the story. */
  async add({ open = true } = {}) {
    if (this.readOnly) return;
    const v = this.viewer;
    // A transition in flight would be captured half-blended: land it first.
    if (this.running) this.cancel();
    const state = captureState(v, this.sky);
    const n = this.project.items.length + 1;
    const it = { id: uid(), name: this.autoName(state, n), caption: "", thumb: "", camera: "follow", transition: null, state };
    this.project.items.push(it);
    this.active = this.project.items.length - 1;
    it.thumb = await this.thumbnail();
    this.changed();
    if (open && this.el.dataset.open !== "true") this.toggle(true);
    this.ui.toast(`Saved view ${n}`, { icon: icon("bookmark"), action: { label: "Undo", run: () => this.remove(this.project.items.indexOf(it), { silent: true }) } });
  }

  autoName(state, n) {
    const t = this.viewer.timeline.frameToTime(state.time.frame);
    const when = Math.abs(t) < 1e-6 ? "today" : this.ui.formatTime(t);
    const sel = this.ui.selection && this.ui.selection.kind !== "member" ? this.viewer.describeObject(this.ui.selection.trace, this.ui.selection.index)?.name : "";
    const where = state.view.mode === "sky" ? "Sky" : sel || "Galaxy";
    return `${where}, ${when}`.replace(/^./, (c) => c.toUpperCase()) || `View ${n}`;
  }

  async update(i) {
    const it = this.project.items[i];
    if (this.running) this.cancel();
    it.state = captureState(this.viewer, this.sky);
    it.thumb = await this.thumbnail();
    this.active = i;
    this.changed();
    this.ui.toast(`Updated “${it.name}”`, { icon: icon("camera") });
  }

  rename(i, nameEl) {
    const it = this.project.items[i];
    const input = h("input", { class: "ov-input", value: it.name, "aria-label": "View name" });
    const target = nameEl || this.strip.children[i]?.querySelector(".ov-card-name");
    if (!target) return;
    target.replaceChildren(input);
    input.focus();
    input.select();
    let finished = false;
    const done = (commit) => {
      // Re-rendering removes the input, which fires blur: finish only once.
      if (finished) return;
      finished = true;
      if (commit && input.value.trim()) { it.name = input.value.trim(); this.changed(); }
      else this.render();
    };
    input.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Enter") done(true);
      if (e.key === "Escape") done(false);
    });
    input.addEventListener("blur", () => done(true));
  }

  editCaption(i) {
    const it = this.project.items[i];
    const ta = h("textarea", { class: "ov-input", rows: "4", style: { height: "auto", padding: "10px", lineHeight: "1.5" } }, it.caption || "");
    ta.addEventListener("keydown", (e) => e.stopPropagation());
    const save = h("button", { class: "ov-btn ov-btn--primary", type: "button", onclick: () => { it.caption = ta.value.trim(); this.changed(); this.ui.closeSheet(); } }, "Save caption");
    this.ui.sheet(`Caption · ${it.name}`, h("div", { class: "ov-dialog-grid" },
      h("div", { class: "ov-note" }, "Shown under the view title while presenting."), ta, h("div", null, save)));
    setTimeout(() => ta.focus(), 30);
  }

  /** Edit the list while the same view (if any) stays the active one. */
  keepActive(fn) {
    const cur = this.project.items[this.active];
    fn();
    this.active = cur ? this.project.items.indexOf(cur) : -1;
  }

  duplicate(i) {
    const it = cloneJson(this.project.items[i]);
    it.id = uid();
    it.name = `${it.name} (copy)`;
    this.keepActive(() => this.project.items.splice(i + 1, 0, it));
    this.changed();
  }

  move(from, to) {
    const items = this.project.items;
    this.keepActive(() => {
      const [it] = items.splice(from, 1);
      items.splice(to, 0, it);
    });
    this.changed();
  }

  remove(i, { silent = false } = {}) {
    if (i < 0 || i >= this.project.items.length) return;
    let it = null;
    this.keepActive(() => { [it] = this.project.items.splice(i, 1); });
    this.changed();
    if (!silent) {
      this.ui.toast(`Deleted “${it.name}”`, { action: { label: "Undo", run: () => {
        this.keepActive(() => this.project.items.splice(Math.min(i, this.project.items.length), 0, it));
        this.changed();
      } } });
    }
  }

  changed() {
    this.dirty = true;
    this.render();
    this.saveDraft();
    this.renderPresenter();
    this.emitDom("oviz:states-changed", { count: this.project.items.length });
  }

  /** DOM events for pages that embed the figure (classic Oviz fired these too). */
  emitDom(type, detail) {
    try { this.ui.root.dispatchEvent(new CustomEvent(type, { detail, bubbles: true })); } catch (_) { /* old browsers */ }
  }

  async thumbnail() {
    const v = this.viewer;
    const r = v.renderer;
    try {
      const W = 320, H = 180;
      const c = document.createElement("canvas");
      c.width = W; c.height = H;
      const ctx = c.getContext("2d");
      ctx.fillStyle = "#000";
      ctx.fillRect(0, 0, W, H);
      const src = r.canvas;
      const scale = Math.max(W / src.width, H / src.height);
      const sw = W / scale, sh = H / scale;
      const sx = (src.width - sw) / 2, sy = (src.height - sh) / 2;
      if (v.state.view.mode === "sky" && this.sky?.sky.ready && v.state.sky.backgroundVisible) {
        const bg = document.createElement("canvas");
        bg.width = src.width; bg.height = src.height;
        await this.sky.sky.drawInto(bg.getContext("2d"), bg.width, bg.height);
        ctx.drawImage(bg, sx, sy, sw, sh, 0, 0, W, H);
      }
      r.render();
      ctx.drawImage(src, sx, sy, sw, sh, 0, 0, W, H);
      return c.toDataURL("image/jpeg", 0.78);
    } catch (err) {
      return "";
    }
  }

  // ------------------------------------------------------------ navigation

  /**
   * The user took the camera mid-transition: land every other property now
   * and leave the camera where it is. (A 3D↔Sky transition lands exactly,
   * camera included, so Sky registration holds; "keep" views never held
   * the camera, so they carry on.)
   */
  takeover() {
    const r = this.running;
    if (!r || r.transition.keepCamera) return;
    this.cancel({ keepCamera: !r.transition.modeChange });
  }

  cancel({ keepCamera = false } = {}) {
    if (!this.running) return;
    const r = this.running;
    this.running = null;
    r.off();
    const pose = keepCamera ? clonePose(this.viewer.renderer.camera.pose) : null;
    r.transition.finish();
    if (pose) {
      Object.assign(this.viewer.renderer.camera.pose, pose);
      this.viewer.state.view.pose = this.viewer.renderer.camera.pose;
      // Stopped mid-flight: the view's anchor holds only if the camera is on it.
      this.viewer.reconcileAnchor?.();
    }
    this.viewer.renderer.release("state-transition");
    // An interrupted transition still arrives: fire the same events so the
    // view mode, panels and Sky layers reflect the applied State.
    this.afterArrive(r.transition);
    r.resolve();
  }

  goTo(i, { instant = false } = {}) {
    const items = this.project.items;
    if (!items.length) return Promise.resolve();
    i = ((i % items.length) + items.length) % items.length;
    const it = items[i];
    const v = this.viewer;
    // Rapid next/previous: the new flight starts from wherever the camera
    // is now, instead of snapping to the interrupted view first.
    if (this.running) this.cancel({ keepCamera: !this.running.transition.modeChange });
    this.active = i;
    this.markActive();
    this.renderPresenter();
    v.timeline.pause();
    const keepCamera = it.camera === "keep";
    const transition = makeTransition(v, this.sky, it.state, { keepCamera });
    const spec = it.transition || this.project.defaultTransition;
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    const dur = instant || reduce ? 0 : Math.max(0, spec.duration_ms) * (transition.modeChange ? 1.35 : 1);
    const ease = easing(spec.easing);
    if (this.presenting) this.renderProgress(dur);
    v.emit("state-transition-start", { index: i, id: it.id });
    this._target = { index: i + 1, id: it.id };
    return this.runTransition(transition, dur, ease);
  }

  /**
   * Apply any captured State (the public `applyState`): a smooth transition
   * by default, or `{ instant: true }`. No saved view becomes active.
   */
  applyState(state, { instant = false, duration_ms: ms, easing: easeName, keepCamera = false } = {}) {
    if (!state || typeof state !== "object" || !state.view) return Promise.resolve();
    const v = this.viewer;
    if (this.running) this.cancel({ keepCamera: !this.running.transition.modeChange });
    this.active = -1;
    this.markActive();
    v.timeline.pause();
    const transition = makeTransition(v, this.sky, state, { keepCamera });
    const spec = this.project.defaultTransition;
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    const dur = instant || reduce ? 0 : Math.max(0, Number(ms ?? spec.duration_ms) || 0) * (transition.modeChange ? 1.35 : 1);
    return this.runTransition(transition, dur, easing(easeName || spec.easing));
  }

  /** Animate a transition over `dur` ms; resolves once it has arrived exactly. */
  runTransition(transition, dur, ease) {
    const v = this.viewer;
    const which = this._target || { index: null, id: null };
    this._target = null;
    this.emitDom("oviz:transition-start", which);
    return new Promise((resolve0) => {
      const resolve = () => { this.emitDom("oviz:transition-end", which); resolve0(); };
      if (dur <= 0) {
        transition.finish();
        this.afterArrive(transition);
        resolve();
        return;
      }
      const start = performance.now();
      const off = v.addAnimator((now) => {
        const raw = clamp((now - start) / dur, 0, 1);
        transition.step(ease(raw), raw);
        if (raw >= 1) {
          this.running = null;
          off();
          transition.finish();
          v.renderer.release("state-transition");
          this.afterArrive(transition);
          resolve();
        }
      });
      v.renderer.hold("state-transition");
      this.running = { off, transition, resolve };
    });
  }

  afterArrive(transition) {
    const v = this.viewer;
    // Views imported without a thumbnail get one the first time they are visited.
    const it = this.project.items[this.active];
    if (it && !it.thumb) {
      setTimeout(async () => {
        if (this.running || this.project.items[this.active] !== it) return;
        it.thumb = await this.thumbnail();
        this.render();
      }, 250);
    }
    const mode = transition.to.view.mode;
    v.emit("viewmode", { mode, from: transition.from.view.mode });
    v.emit("viewmode-settled", { mode });
    v.emit("state-applied", {});
    v.emit("time", v.timeline);
    this.ui.renderColorbars();
    this.ui.layers.render();
    this.ui.syncViewSeg();
  }

  markActive() {
    [...this.strip.querySelectorAll(".ov-card")].forEach((c, j) => { c.dataset.active = String(j === this.active); });
  }

  next() { return this.goTo(this.active + 1); }
  previous() { return this.goTo(this.active < 0 ? 0 : this.active - 1); }

  /** Back to the figure as it opened (classic "Original"), with the usual transition. */
  goHome({ instant = false } = {}) {
    if (!this.homeState) return Promise.resolve();
    const v = this.viewer;
    if (this.running) this.cancel({ keepCamera: !this.running.transition.modeChange });
    this.active = -1;
    this.markActive();
    this.renderPresenter();
    v.timeline.pause();
    const t = makeTransition(v, this.sky, this.homeState);
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    const spec = this.project.defaultTransition;
    const dur = instant || reduce ? 0 : Math.max(0, spec.duration_ms) * (t.modeChange ? 1.35 : 1);
    v.emit("state-transition-start", { index: -1, id: "original" });
    this._target = { index: 0, id: "original" };
    return this.runTransition(t, dur, easeIO);
  }

  // ------------------------------------------------------------ presentation

  buildPresenter() {
    const ui = this.ui;
    this.pTitle = h("div", { class: "ov-present-title" });
    this.pCaption = h("div", { class: "ov-present-caption" });
    this.pDots = h("div", { class: "ov-present-dots" });
    this.pPrev = iconButton("stepBack", "Previous view", () => this.previous(), { shortcut: "←" });
    this.pNext = iconButton("stepFwd", "Next view", () => this.next(), { shortcut: "→" });
    this.pExit = h("button", { class: "ov-btn ov-btn--ghost ov-btn--sm", type: "button", onclick: () => this.present(false) }, this.readOnly ? "Explore" : "Exit", kbd("Esc"));
    // The time of the view (and of playing time), as the classic time bar showed.
    this.pTime = h("div", { class: "ov-present-time", "aria-live": "off" });
    this.pText = h("div", { class: "ov-present-text" }, this.pTitle, this.pCaption, this.pTime, this.pDots);
    this.presenter = h("div", { class: "ov-present ov-glass", "data-open": "false", role: "region", "aria-label": "Presentation" },
      this.pPrev,
      this.pText,
      this.pNext,
      h("div", { class: "ov-present-exit" }, this.pExit));
    // Story-style progress along the top; the current segment fills while
    // the camera flies to the view.
    this.pProgress = h("div", { class: "ov-present-progress", "data-open": "false", role: "group", "aria-label": "Views" });
    // A 3D reel of the views above the caption bar, shown while the pointer
    // moves and tucked away when the presenter is idle.
    this.pReel = h("div", { class: "ov-present-reel", "data-open": "false", "data-idle": "true" });
    ui.ui.append(this.pReel, this.presenter, this.pProgress);
    let idleTimer = 0;
    const wake = () => {
      if (!this.presenting) return;
      this.pReel.dataset.idle = "false";
      clearTimeout(idleTimer);
      idleTimer = setTimeout(() => { this.pReel.dataset.idle = "true"; }, 2600);
    };
    ui.root.addEventListener("pointermove", wake, { passive: true });
    this.pReel.addEventListener("pointerenter", () => clearTimeout(idleTimer));
    this.pReel.addEventListener("pointerleave", wake);
    const tl = this.viewer.timeline;
    this.pTime.hidden = tl.count <= 1;
    tl.on(() => { if (this.presenting) this.syncPresenterTime(); });
    window.addEventListener("resize", () => { if (this.presenting) requestAnimationFrame(() => this.liftKey()); });
  }

  syncPresenterTime() {
    const tl = this.viewer.timeline;
    if (tl.count <= 1) return;
    const t = tl.time;
    const text = Math.abs(t) < 1e-6 ? "Present day" : this.ui.formatTime(t);
    const playing = tl.playing ? (tl.direction < 0 ? " · playing backward" : " · playing") : "";
    this.pTime.textContent = `${text}${playing}`;
    // A view named after its time needs no second readout until time moves.
    const title = this.pTitle.textContent || "";
    this.pTime.hidden = !tl.playing && title.includes(this.ui.formatTime(t));
  }

  present(on, { startAt = null, instant = false } = {}) {
    const items = this.project.items;
    if (on && !items.length) return;
    this.presenting = on;
    this.ui.root.dataset.presenting = String(on);
    this.presenter.dataset.open = String(on);
    this.pProgress.dataset.open = String(on);
    this.pReel.dataset.open = String(on);
    this.pReel.dataset.idle = "true";
    this.ui.setZen(on);
    // The key shows only the layers in view while presenting.
    this.ui.legend?.render();
    if (!on) this.liftKey();
    if (on) {
      this.toggle(false);
      const i = startAt ?? (this.active >= 0 ? this.active : 0);
      this.goTo(i, { instant });
      this.renderPresenter();
      this.syncPresenterTime();
    }
  }

  renderPresenter() {
    const items = this.project.items;
    const it = items[this.active];
    const changed = this._shown !== it;
    this._shown = it;
    this.pTitle.textContent = it ? it.name : "";
    this.pCaption.textContent = it?.caption || "";
    this.pCaption.hidden = !it?.caption;
    clear(this.pDots);
    items.forEach((x, j) => {
      const d = h("button", { class: "ov-present-dot", type: "button", "aria-label": `Go to ${x.name}`, "aria-current": j === this.active ? "step" : null, onclick: () => this.goTo(j) });
      this.pDots.append(d);
    });
    if (changed && this.presenting) {
      // Replay the caption entrance for each new view.
      this.pText.classList.remove("ov-anim");
      void this.pText.offsetWidth;
      this.pText.classList.add("ov-anim");
    }
    this.renderProgress();
    this.renderReel();
    requestAnimationFrame(() => {
      this.ui.root.style.setProperty("--ov-present-h", `${this.presenter.offsetHeight}px`);
      this.liftKey();
    });
  }

  /**
   * Where the presenter would cover the key and scale bar (phones, narrow
   * windows), they step up just above it; elsewhere they stay put.
   */
  liftKey() {
    const root = this.ui.root;
    const bl = this.ui.bottomLeft;
    const cur = parseFloat(root.style.getPropertyValue("--ov-key-lift")) || 0;
    let lift = 0;
    const host = bl?.offsetParent;
    if (this.presenting && host) {
      // Layout offsets ignore the lift's own (animating) transform.
      const a = this.presenter.getBoundingClientRect(), hr = host.getBoundingClientRect();
      const left = hr.left + bl.offsetLeft, bottom = hr.top + bl.offsetTop + bl.offsetHeight;
      if (a.left < left + bl.offsetWidth && left < a.right) lift = Math.max(0, Math.round(bottom - a.top + 10));
    }
    if (lift !== cur) root.style.setProperty("--ov-key-lift", `${lift}px`);
  }

  renderProgress(durationMs = 0) {
    const items = this.project.items;
    if (this.pProgress.childElementCount !== items.length) {
      this.pProgress.replaceChildren(...items.map((x, j) => h("button", { class: "ov-present-seg", type: "button", "aria-label": `Go to ${x.name}`, onclick: () => this.goTo(j) }, h("span", null, h("i")))));
    }
    [...this.pProgress.children].forEach((seg, j) => {
      const fill = seg.querySelector("i");
      seg.setAttribute("aria-current", j === this.active ? "step" : "false");
      if (j === this.active && durationMs > 0) {
        fill.style.transition = "none";
        fill.style.width = "0%";
        void fill.offsetWidth;
        fill.style.transition = `width ${durationMs}ms linear`;
        fill.style.width = "100%";
      } else if (j !== this.active || !fill.style.width) {
        fill.style.transition = "none";
        fill.style.width = j <= this.active ? "100%" : "0%";
      }
    });
  }

  renderReel() {
    const items = this.project.items;
    const reel = this.pReel;
    if (reel.childElementCount !== items.length || [...reel.children].some((c, j) => c._item !== items[j] || c._thumb !== (items[j].thumb || ""))) {
      reel.replaceChildren(...items.map((it, j) => {
        const card = h("button", { class: "ov-reel-card", type: "button", "aria-label": `Go to ${it.name}`, onclick: () => this.goTo(j) },
          h("b", null, String(j + 1)), h("span", null, it.name));
        if (it.thumb) card.style.backgroundImage = `url("${it.thumb}")`;
        else card.append(h("em", null, it.state.view?.mode === "sky" ? "Sky" : "3D"));
        card._item = it;
        card._thumb = it.thumb || "";
        return card;
      }));
    }
    const n = items.length;
    [...reel.children].forEach((card, j) => {
      const o = j - Math.max(this.active, 0);
      const a = Math.abs(o);
      card.style.transform = `translateX(${o * 150}px) translateZ(${-a * 90}px) rotateY(${Math.max(-50, Math.min(50, -o * 26))}deg)`;
      card.style.opacity = a > 3 ? "0" : String(1 - a * 0.2);
      card.style.zIndex = String(n - a);
      card.style.pointerEvents = a > 3 ? "none" : "";
      card.setAttribute("aria-current", o === 0 ? "step" : "false");
    });
  }

  onKey(e) {
    const mod = e.metaKey || e.ctrlKey;
    if (mod && (e.key === "s" || e.key === "S")) { this.saveDocument(); return true; }
    if (mod) return false;
    const k = e.key;
    if (this.presenting) {
      // Classic presentation keys: arrows move between views; Space still
      // plays time; Esc or P leaves.
      if (k === "ArrowRight" || k === "PageDown") { if (!e.repeat) this.next(); return true; }
      if (k === "ArrowLeft" || k === "PageUp") { if (!e.repeat) this.previous(); return true; }
      if (k === "Home") { this.goTo(0); return true; }
      if (k === "End") { this.goTo(this.project.items.length - 1); return true; }
      if (k === "Escape" || k === "p" || k === "P") { this.present(false); return true; }
      return false;
    }
    if (k === "p" || k === "P") {
      if (!this.project.items.length) this.ui.toast("Save a view with N first, then press P to present");
      else this.present(true);
      return true;
    }
    if ((k === "n" || k === "N") && !this.readOnly) { this.add({ open: false }); return true; }
    if (k === "y" || k === "Y") { this.toggle(); return true; }
    return false;
  }

  // ------------------------------------------------------------ palette & commands

  paletteItems() {
    return this.project.items.map((it, i) => ({
      title: it.name,
      sub: it.caption || `View ${i + 1}`,
      icon: "bookmark",
      hint: String(i + 1),
      run: () => this.goTo(i),
    }));
  }

  commands() {
    const list = [
      { title: "Present story", icon: "present", shortcut: "P", run: () => this.present(true) },
      { title: "Open views & story", icon: "bookmark", shortcut: "Y", run: () => this.toggle(true) },
      { title: "Go to the original view", sub: "The figure as it opened: camera, time, layers and Sky", icon: "home", keywords: "reset states story start initial first", run: () => this.goHome() },
    ];
    // ("Copy link to this view" carries the views already.)
    list.push(...this.linkMenuItems().slice(1).map((x) => ({ title: x.label, icon: x.icon, keywords: "share url link presentation story views", run: x.run })));
    if (!this.readOnly) {
      list.unshift({ title: "Save current view", icon: "plus", shortcut: "N", pinned: true, run: () => this.add() });
      list.push(...this.exportMenuItems().filter((x) => x !== "-").map((x) => ({ title: x.label, icon: x.icon, run: x.run })));
    }
    return list;
  }

  /**
   * Views that came with a shared link (see app/viewhash.js encodeViewsPart):
   * they replace the story for this visit, and autosave apart from the
   * views this browser keeps for the figure.
   */
  loadShared({ items, defaultTransition }) {
    if (!Array.isArray(items) || !items.length) return;
    this.sharedSession = true;
    this.project.items = items.map((it) => ({
      id: uid(), name: it.name, caption: it.caption || "", thumb: "",
      camera: it.camera === "keep" ? "keep" : "follow",
      transition: readTransition(it.transition), state: it.state,
    }));
    if (defaultTransition) this.project.defaultTransition = readTransition(defaultTransition, this.project.defaultTransition);
    this.active = -1;
    this.dirty = false;
    this.ui.statesBtn.hidden = false;
    this.render();
    this.renderPresenter();
  }

  /** Links that carry the views: to explore them, or to open presenting them. */
  linkMenuItems() {
    const n = this.project.items.length;
    if (!n) return [];
    return [
      { label: "Copy link with these views", icon: "share", run: () => this.ui.copyViewLink() },
      { label: "Copy presentation link", icon: "present", run: () => this.ui.copyPresentationLink() },
    ];
  }

  exportMenuItems() {
    if (this.readOnly) return [];
    return [
      { label: "Save figure (editable)", icon: "download", shortcut: `${MOD === "⌘" ? "⌘" : "Ctrl"}S`, run: () => this.saveDocument() },
      { label: "Export presentation (present-only)", icon: "present", run: () => this.exportHtml({ presentOnly: true }) },
      { label: "Export current view only", icon: "camera", run: () => this.exportHtml({ currentOnly: true }) },
    ];
  }

  manifestForExport({ presentOnly = false, currentOnly = false } = {}) {
    const m = cloneJson(this.viewer.manifest);
    // Saved files reopen the way they were being viewed (theme and mode).
    m.viewer = { ...(m.viewer || {}), mode: this.ui.mode };
    if (currentOnly) {
      const cur = captureState(this.viewer, this.sky);
      m.states = exportStatesBlock({ ...this.project, projectId: uid(), items: [{ id: uid(), name: "View", caption: "", thumb: "", camera: "follow", transition: null, state: cur }] }, { defaultMode: "edit" });
      m.states.apply_first_on_load = true;
      m.states.items_hidden = true;
      return m;
    }
    m.states = exportStatesBlock(this.project, { presentOnly, defaultMode: presentOnly ? "present" : "edit" });
    return m;
  }

  async exportHtml(opts = {}) {
    const m = this.manifestForExport(opts);
    const text = buildExportHtml(m, { theme: this.ui.theme, mode: this.ui.mode });
    const base = slugName(this.viewer.manifest.title || "oviz-figure");
    const suffix = opts.presentOnly ? "-presentation" : opts.currentOnly ? "-view" : "";
    const res = await saveHtml(text, `${base}${suffix}.html`);
    if (res.ok) this.ui.toast(`Exported ${res.name}`, { icon: icon("download") });
  }

  async saveDocument() {
    if (this.readOnly) return;
    const m = this.manifestForExport({});
    const text = buildExportHtml(m, { theme: this.ui.theme, mode: this.ui.mode });
    const res = await saveHtml(text, `${slugName(this.viewer.manifest.title || "oviz-figure")}.html`, { reuseHandle: true });
    if (res.ok) {
      this.dirty = false;
      this.clearDraft();
      this.render();
      this.ui.toast(`Saved ${res.name}`, { icon: icon("download") });
    }
  }

  // ------------------------------------------------------------ drafts

  draftKey() {
    // Views opened from a shared link autosave apart from this browser's own.
    return `oviz.draft.${this.project.projectId}${this.sharedSession ? ".shared" : ""}`;
  }

  saveDraft() {
    if (this.readOnly || this.project.autosave === false) return;
    this.draftSaved = false;
    clearTimeout(this._draftTimer);
    this._draftTimer = setTimeout(() => this.flushDraft(), 300);
  }

  /** Write the draft now (IndexedDB, else localStorage); say so once if it fails. */
  flushDraft() {
    clearTimeout(this._draftTimer);
    if (this.readOnly || this.project.autosave === false || !this.dirty) return Promise.resolve();
    const payload = { savedAt: Date.now(), items: this.project.items, defaultTransition: this.project.defaultTransition };
    return writeDraft(this.draftKey(), JSON.stringify(payload)).then((ok) => {
      this.draftSaved = ok;
      if (!ok && !this._draftWarned) {
        this._draftWarned = true;
        this.ui.toast(`Your views could not be autosaved in this browser: save the figure (${MOD === "⌘" ? "⌘" : "Ctrl"} S) to keep them`, { icon: icon("bookmark"), ms: 6000 });
      }
    });
  }

  clearDraft() {
    deleteDraft(this.draftKey());
  }

  async restoreDraft() {
    if (this.readOnly || this.project.autosave === false || this.sharedSession) return;
    const raw = await readDraft(this.draftKey());
    // A shared link's views arrived meanwhile: they win.
    if (!raw || this.sharedSession) return;
    try {
      const d = JSON.parse(raw);
      if (!Array.isArray(d.items)) return;
      // Any edit counts (captions, updated views, camera, timing), not only
      // the list of names.
      const same = JSON.stringify(d.items) === JSON.stringify(this.project.items)
        && (!d.defaultTransition || JSON.stringify(d.defaultTransition) === JSON.stringify(this.project.defaultTransition));
      if (same) return;
      const backup = this.project.items;
      this.project.items = d.items;
      if (d.defaultTransition) this.project.defaultTransition = d.defaultTransition;
      this.dirty = true;
      this.ui.toast("Restored your unsaved views from this browser", {
        icon: icon("bookmark"),
        ms: 6000,
        action: { label: "Discard", run: () => { this.project.items = backup; this.clearDraft(); this.dirty = false; this.render(); } },
      });
    } catch (_) { /* corrupt draft */ }
  }

  // ------------------------------------------------------------ public API

  /** A view by 1-based index, id or name; −1 when there is none. */
  indexOf(target) {
    const items = this.project.items;
    if (typeof target === "number") return target >= 1 && target <= items.length ? target - 1 : -1;
    return items.findIndex((x) => x.id === target || x.name === target);
  }

  api() {
    return {
      list: () => this.project.items.map((it, i) => ({ index: i + 1, id: it.id, name: it.name })),
      goTo: (target) => {
        if (target === "original") return this.goHome();
        const i = typeof target === "number" ? target - 1 : this.project.items.findIndex((x) => x.id === target || x.name === target);
        return i >= 0 && i < this.project.items.length ? this.goTo(i) : Promise.resolve();
      },
      original: (opts) => this.goHome(opts),
      next: () => this.next(),
      previous: () => this.previous(),
      add: () => this.add({ open: false }),
      // Editing by 1-based index, id or name (classic Oviz.get(id).states).
      update: (target) => { const i = this.indexOf(target); return i >= 0 ? this.update(i) : Promise.resolve(); },
      rename: (target, name) => { const i = this.indexOf(target); if (i >= 0 && String(name || "").trim()) { this.project.items[i].name = String(name).trim(); this.changed(); } },
      duplicate: (target) => { const i = this.indexOf(target); if (i >= 0) this.duplicate(i); },
      move: (target, to) => { const i = this.indexOf(target); const j = Number(to) - 1; if (i >= 0 && Number.isInteger(j) && j >= 0 && j < this.project.items.length) this.move(i, j); },
      remove: (target) => { const i = this.indexOf(target); if (i >= 0) this.remove(i, { silent: true }); },
      setCameraBehavior: (target, behavior) => { const i = this.indexOf(target); if (i >= 0) { this.project.items[i].camera = behavior === "keep" ? "keep" : "follow"; this.changed(); } },
      save: () => this.saveDocument(),
      activeIndex: () => (this.active >= 0 ? this.active + 1 : 0),
      capture: () => captureState(this.viewer, this.sky),
      present: (on = true) => this.present(on),
      exportHtml: (opts = {}) => (opts.download === false ? buildExportHtml(this.manifestForExport(opts), { theme: this.ui.theme }) : this.exportHtml(opts)),
    };
  }
}

function easeIO(t) {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
}

function slugName(s) {
  return String(s).trim().toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "oviz-figure";
}

/** Stable id for a figure (title + blob fingerprint), for per-figure storage. */
export function figureId(m) {
  let hash = 2166136261;
  const text = `${m.title}|${(m.blobs || []).map((b) => `${b.id}:${b.bytes}`).join(",")}`;
  for (let i = 0; i < text.length; i++) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }
  return (hash >>> 0).toString(36);
}
