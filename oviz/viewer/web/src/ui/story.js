// "Views & story": saved States as a filmstrip, presentation mode, drafts
// and self-export.

import { h, icon, iconButton, clear, kbd, MOD, localStorageGet, localStorageSet } from "./dom.js";
import { slider, toggle } from "./controls.js";
import { captureState, importStates, exportStatesBlock, makeTransition, easing, uid } from "../app/states.js";
import { cloneJson } from "../app/state.js";
import { buildExportHtml, saveHtml } from "../app/export.js";
import { clamp } from "../core/math.js";

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
    v.on("transition", () => {});
    // Cancel a running transition when the user grabs the camera.
    v.canvas.addEventListener("pointerdown", () => this.cancel(), true);
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
    }
    this.homeState = captureState(v, this.sky);
    this.restoreDraft();
    this.render();
    this.ui.dock.setMarkers?.([]);
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
    this.addBtn = h("button", { class: "ov-btn ov-btn--primary ov-btn--sm", type: "button", onclick: () => this.add() }, icon("plus"), "Save view");
    this.presentBtn = h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.present(true) }, icon("present"), "Present");
    this.exportBtn = h("button", { class: "ov-btn ov-btn--sm ov-btn--ghost", type: "button", onclick: (e) => ui.menu(e.currentTarget, this.exportMenuItems(), { title: "Export", above: true }) }, icon("download"), "Export");
    this.settingsBtn = iconButton("sliders", "Transition timing", (e) => this.timingMenu(e.currentTarget));
    const head = h("div", { class: "ov-story-head" },
      h("span", { class: "ov-story-title" }, "Story"), this.count,
      h("span", { class: "ov-grow" }),
      this.readOnly ? null : this.addBtn, this.presentBtn, this.readOnly ? null : this.exportBtn, this.readOnly ? null : this.settingsBtn,
      iconButton("close", "Close", () => this.toggle(false), { shortcut: "Y" }));
    this.strip = h("div", { class: "ov-story-strip", role: "list" });
    this.status = h("div", { class: "ov-story-status" });
    this.el.append(head, this.strip, this.status);
    ui.ui.append(this.el);
  }

  toggle(open = this.el.dataset.open !== "true") {
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
      return;
    }
    items.forEach((it, i) => this.strip.append(this.card(it, i)));
    if (!this.readOnly) {
      this.strip.append(h("button", { class: "ov-story-add", type: "button", onclick: () => this.add(), "data-tip": "Save current view  N" }, icon("plus")));
    }
    this.status.textContent = this.dirty ? "Unsaved changes · autosaved in this browser" : "";
    this.ui.dock.setMarkers?.(items.map((it) => this.viewer.timeline.frameToTime(it.state.time?.frame ?? 0)));
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
      if (e.key === "Enter") { this.goTo(i); e.preventDefault(); }
      if ((e.key === "Delete" || e.key === "Backspace") && !this.readOnly) { this.remove(i); e.preventDefault(); }
      e.stopPropagation();
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
        const from = Number(e.dataTransfer.getData("text/oviz-state"));
        card.dataset.over = "false";
        if (Number.isFinite(from) && from !== i) this.move(from, i);
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
      { label: "Duplicate", icon: "copy", run: () => this.duplicate(i) },
      "-",
      { label: "Delete", icon: "trash", run: () => this.remove(i) },
    ], { title: it.name, above: true, align: "left" });
  }

  timingMenu(anchor) {
    const dt = this.project.defaultTransition;
    const body = h("div", { class: "ov-pop-body" },
      slider({ label: "Transition length", min: 0, max: 6, step: 0.1, value: dt.duration_ms / 1000, format: (x) => `${x.toFixed(1)} s`, onChange: (x) => { dt.duration_ms = Math.round(x * 1000); this.changed(); } }),
      toggle({ label: "Autosave drafts in this browser", checked: this.project.autosave !== false, onChange: (x) => { this.project.autosave = x; this.changed(); } }),
    );
    this.ui.menu(anchor, [{ body }], { title: "Story settings", above: true });
  }

  // ------------------------------------------------------------ editing

  async add() {
    if (this.readOnly) return;
    const v = this.viewer;
    const state = captureState(v, this.sky);
    const n = this.project.items.length + 1;
    const it = { id: uid(), name: this.autoName(state, n), caption: "", thumb: "", camera: "follow", transition: null, state };
    this.project.items.push(it);
    this.active = this.project.items.length - 1;
    it.thumb = await this.thumbnail();
    this.changed();
    if (this.el.dataset.open !== "true") this.toggle(true);
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
    const done = (commit) => {
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

  duplicate(i) {
    const it = cloneJson(this.project.items[i]);
    it.id = uid();
    it.name = `${it.name} (copy)`;
    this.project.items.splice(i + 1, 0, it);
    this.changed();
  }

  move(from, to) {
    const items = this.project.items;
    const [it] = items.splice(from, 1);
    items.splice(to, 0, it);
    this.active = to;
    this.changed();
  }

  remove(i, { silent = false } = {}) {
    if (i < 0) return;
    const [it] = this.project.items.splice(i, 1);
    if (this.active >= this.project.items.length) this.active = this.project.items.length - 1;
    this.changed();
    if (!silent) this.ui.toast(`Deleted “${it.name}”`, { action: { label: "Undo", run: () => { this.project.items.splice(i, 0, it); this.changed(); } } });
  }

  changed() {
    this.dirty = true;
    this.render();
    this.saveDraft();
    this.renderPresenter();
  }

  async thumbnail() {
    const v = this.viewer;
    const r = v.renderer;
    try {
      const W = 320, H = 180;
      const c = document.createElement("canvas");
      c.width = W; c.height = H;
      const ctx = c.getContext("2d");
      ctx.fillStyle = "#05070b";
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

  cancel() {
    if (!this.running) return;
    const r = this.running;
    this.running = null;
    r.off();
    r.transition.finish();
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
    if (this.running) this.cancel();
    this.active = i;
    this.markActive();
    this.renderPresenter();
    v.timeline.pause();
    this.ui.following = null;
    const keepCamera = it.camera === "keep";
    const transition = makeTransition(v, this.sky, it.state, { keepCamera });
    const spec = it.transition || this.project.defaultTransition;
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    const dur = instant || reduce ? 0 : Math.max(0, spec.duration_ms) * (transition.modeChange ? 1.35 : 1);
    const ease = easing(spec.easing);
    v.emit("state-transition-start", { index: i, id: it.id });
    return new Promise((resolve) => {
      if (dur <= 0) {
        transition.finish();
        this.afterArrive(transition);
        resolve();
        return;
      }
      const start = performance.now();
      const off = v.addAnimator((now) => {
        const raw = clamp((now - start) / dur, 0, 1);
        transition.step(ease(raw));
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

  goHome() {
    if (!this.homeState) return;
    const v = this.viewer;
    this.active = -1;
    this.markActive();
    const t = makeTransition(v, this.sky, this.homeState);
    const start = performance.now();
    const off = v.addAnimator((now) => {
      const raw = clamp((now - start) / 1100, 0, 1);
      t.step(easeIO(raw));
      if (raw >= 1) { off(); t.finish(); this.afterArrive(t); v.renderer.release("state-transition"); }
    });
    v.renderer.hold("state-transition");
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
    this.presenter = h("div", { class: "ov-present ov-glass", "data-open": "false", role: "region", "aria-label": "Presentation" },
      this.pPrev,
      h("div", { class: "ov-present-text" }, this.pTitle, this.pCaption, this.pDots),
      this.pNext,
      h("div", { class: "ov-present-exit" }, this.pExit));
    ui.ui.append(this.presenter);
  }

  present(on, { startAt = null, instant = false } = {}) {
    const items = this.project.items;
    if (on && !items.length) return;
    this.presenting = on;
    this.ui.root.dataset.presenting = String(on);
    this.presenter.dataset.open = String(on);
    this.ui.setZen(on);
    if (on) {
      this.toggle(false);
      const i = startAt ?? (this.active >= 0 ? this.active : 0);
      this.goTo(i, { instant });
      this.renderPresenter();
    }
  }

  renderPresenter() {
    const items = this.project.items;
    const it = items[this.active];
    this.pTitle.textContent = it ? it.name : "";
    this.pCaption.textContent = it?.caption || "";
    this.pCaption.hidden = !it?.caption;
    clear(this.pDots);
    items.forEach((x, j) => {
      const d = h("button", { class: "ov-present-dot", type: "button", "aria-label": `Go to ${x.name}`, "aria-current": j === this.active ? "step" : null, onclick: () => this.goTo(j) });
      this.pDots.append(d);
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
    if ((k === "n" || k === "N") && !this.readOnly) { this.add(); return true; }
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
    ];
    if (!this.readOnly) {
      list.unshift({ title: "Save current view", icon: "plus", shortcut: "N", pinned: true, run: () => this.add() });
      list.push(...this.exportMenuItems().filter((x) => x !== "-").map((x) => ({ title: x.label, icon: x.icon, run: x.run })));
    }
    return list;
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
    if (currentOnly) {
      const cur = captureState(this.viewer, this.sky);
      m.states = exportStatesBlock({ ...this.project, items: [{ id: uid(), name: "View", caption: "", thumb: "", camera: "follow", transition: null, state: cur }] }, { defaultMode: "edit" });
      m.states.apply_first_on_load = true;
      m.states.items_hidden = true;
      return m;
    }
    m.states = exportStatesBlock(this.project, { presentOnly, defaultMode: presentOnly ? "present" : "edit" });
    return m;
  }

  async exportHtml(opts = {}) {
    const m = this.manifestForExport(opts);
    const text = buildExportHtml(m, { theme: this.ui.theme });
    const base = slugName(this.viewer.manifest.title || "oviz-figure");
    const suffix = opts.presentOnly ? "-presentation" : opts.currentOnly ? "-view" : "";
    const res = await saveHtml(text, `${base}${suffix}.html`);
    if (res.ok) this.ui.toast(`Exported ${res.name}`, { icon: icon("download") });
  }

  async saveDocument() {
    if (this.readOnly) return;
    const m = this.manifestForExport({});
    const text = buildExportHtml(m, { theme: this.ui.theme });
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
    return `oviz.draft.${this.project.projectId}`;
  }

  saveDraft() {
    if (this.readOnly || this.project.autosave === false) return;
    clearTimeout(this._draftTimer);
    this._draftTimer = setTimeout(() => {
      const payload = { savedAt: Date.now(), items: this.project.items, defaultTransition: this.project.defaultTransition };
      localStorageSet(this.draftKey(), JSON.stringify(payload));
    }, 300);
  }

  clearDraft() {
    try { localStorage.removeItem(this.draftKey()); } catch (_) { /* ignore */ }
  }

  restoreDraft() {
    if (this.readOnly || this.project.autosave === false) return;
    const raw = localStorageGet(this.draftKey());
    if (!raw) return;
    try {
      const d = JSON.parse(raw);
      if (!Array.isArray(d.items)) return;
      const embedded = JSON.stringify(this.project.items.map((x) => [x.id, x.name]));
      const draft = JSON.stringify(d.items.map((x) => [x.id, x.name]));
      if (embedded === draft && d.items.length === this.project.items.length) return;
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

  api() {
    return {
      list: () => this.project.items.map((it, i) => ({ index: i + 1, id: it.id, name: it.name })),
      goTo: (target) => this.goTo(typeof target === "number" ? target - 1 : this.project.items.findIndex((x) => x.id === target || x.name === target)),
      next: () => this.next(),
      previous: () => this.previous(),
      add: () => this.add(),
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

function figureId(m) {
  let hash = 2166136261;
  const text = `${m.title}|${(m.blobs || []).map((b) => `${b.id}:${b.bytes}`).join(",")}`;
  for (let i = 0; i < text.length; i++) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }
  return (hash >>> 0).toString(36);
}
