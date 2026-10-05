// Layers panel: group presets, the data layers (traces and volumes, one
// list), the distribution filter, the Sky background and reference guides.

import { h, icon, iconButton, clear, fmtInt, fmt } from "./dom.js";
import { slider, toggle, select, miniSeg, colorPicker, numberInput } from "./controls.js";
import { lutGradient } from "../core/color.js";
import { applyGroupDefaults, groupMode, volumeGroupMode, clampVolumeState, toggleAllPlan, soloPlan } from "../app/state.js";

export class LayersPanel {
  constructor(ui) {
    this.ui = ui;
    this.viewer = ui.viewer;
    this.manifest = ui.viewer.manifest;
    this.expanded = null;
    this.rows = new Map();
    // Section show/hide buttons, and what each hid (so showing restores it).
    this.sectionEyes = [];
    this.snapshots = new Map();
    this.el = h("section", { class: "ov-panel ov-layers ov-glass ov-chrome", "aria-label": "Layers", "data-open": "true" });
    this.head = h("div", { class: "ov-panel-head" },
      h("span", { class: "ov-panel-title" }, "Layers"),
      iconButton("close", "Hide layers", () => ui.setLayersOpen(false), { shortcut: "⇧L" }),
    );
    this.body = h("div", { class: "ov-panel-body" });
    this.el.append(this.head, this.body);
    this.render();
    const v = this.viewer;
    v.on("style", () => this.sync());
    v.on("volume", () => this.sync());
    v.on("global", () => this.sync());
    v.on("volume-ready", () => this.sync());
    // Volumes shown only at the present day say so while time is elsewhere.
    let timeSync = 0;
    v.on("time", () => {
      if (timeSync) return;
      timeSync = requestAnimationFrame(() => { timeSync = 0; v._resolveParams?.(); this.sync(); });
    });
    v.on("viewmode", () => this.render());
    v.on("state-applied", () => {
      // A remembered mix belongs to the view it was hidden from.
      this.snapshots.clear();
      this.render();
    });
  }

  get lutFor() {
    return this.ui.luts;
  }

  // ---------------------------------------------------------------- render

  render() {
    clear(this.body);
    this.rows.clear();
    this.sectionEyes = [];
    // Rows are rebuilt closed, so no editor is open any more.
    this.expanded = null;
    const v = this.viewer;
    const m = this.manifest;
    // "All" leads the list, as in the classic group menu.
    const groups = [...(m.groups?.order || [])].sort((a, b) => (String(b).toLowerCase() === "all") - (String(a).toLowerCase() === "all"));
    if (groups.length > 1) {
      const sel = h("select", { "aria-label": "Layer group" }, groups.map((g) => h("option", { value: g, selected: g === v.state.group }, g)));
      sel.addEventListener("change", () => this.setGroup(sel.value));
      this.body.append(h("div", { class: "ov-group" }, h("label", null, "View"), sel));
    }
    const traces = this.groupTraces();
    const volumes = this.groupVolumes();
    if (traces.length || volumes.length) {
      // Volumes are data layers like the traces: one list, numbered 1–9 in
      // order, and one eye (the T key) that shows or hides all of it.
      const sec = h("div", { class: "ov-section" },
        h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Data"), this.sectionEye("all")),
      );
      traces.forEach((t, i) => sec.append(...this.traceRow(t, i)));
      volumes.forEach((spec, j) => sec.append(...this.volumeRow(spec, traces.length + j)));
      this.body.append(sec);
    }
    const filterSection = this.ui.filter?.section();
    if (filterSection) this.body.append(filterSection);
    if (v.state.view.mode === "sky" && this.ui.sky) {
      const sec = this.ui.sky.layersSection?.();
      if (sec) this.body.append(sec);
    }
    this.body.append(this.referenceSection());
    this.sync();
  }

  setGroup(group) {
    const v = this.viewer;
    applyGroupDefaults(v.state, this.manifest, group);
    v.volumes.invalidate();
    v.renderer.invalidate();
    v.invalidatePick();
    v.emit("style", { key: "*" });
    v.emit("volume", { key: "*" });
    this.render();
  }

  // ---------------------------------------------------------------- show / hide all

  /**
   * The data traces the current group lists (the panel's Data rows). In Sky
   * view the eye is at the Sun, so the Sun is not a layer there. Guides
   * (annotations the figure marks as such) are listed under Guides instead.
   */
  groupTraces() {
    const v = this.viewer;
    const m = this.manifest;
    const sky = v.state.view.mode === "sky";
    const sun = sky ? v.world.sunTrace : null;
    // 3D-only context (skyHidden) leaves the list in Sky view, like the Sun.
    return m.traces.filter((t) => t.showInLegend && t.role !== "guide" && t.key !== sun && !(sky && t.skyHidden) && groupMode(m, v.state.group, t.key) !== false);
  }

  /** The figure's own guides (labels, scale bars, model curves): one switch each. */
  guideTraces() {
    const v = this.viewer;
    const sky = v.state.view.mode === "sky";
    return this.manifest.traces.filter((t) => t.role === "guide" && !(sky && t.skyHidden) && groupMode(this.manifest, v.state.group, t.key) !== false);
  }

  /** The volumes the current group lists (every volume when groups list none). */
  groupVolumes() {
    const m = this.manifest;
    return uniqueVolumes(m.volumes || []).filter((s) => volumeGroupMode(m, this.viewer.state.group, s) !== false);
  }

  /**
   * The layers a show/hide-all covers, in legend order: data traces, then
   * volumes. `scope` is "all" (the T key), "traces" or "volumes".
   */
  toggleItems(scope = "all") {
    const items = [];
    if (scope !== "volumes") for (const t of this.groupTraces()) items.push({ id: `trace:${t.key}`, kind: "trace", key: t.key, name: t.name });
    if (scope !== "traces") {
      for (const s of this.groupVolumes()) items.push({ id: `volume:${s.stateKey}`, kind: "volume", key: s.stateKey, name: s.name });
    }
    return items;
  }

  itemVisible(it) {
    const v = this.viewer;
    return it.kind === "volume" ? !!v.state.volumes[it.key]?.visible : !!v.state.traces[it.key]?.visible;
  }

  setItemVisible(it, on) {
    const v = this.viewer;
    if (this.itemVisible(it) === on) return;
    if (it.kind === "volume") v.setVolumeState(it.key, { visible: on });
    else v.setTraceStyle(it.key, { visible: on });
  }

  /**
   * Show or hide every layer in `scope`, as the classic legend did: hiding
   * remembers what was on, so showing again restores that mix (a volume
   * that was on comes back with the traces) rather than turning on
   * everything. A section with no memory of its own restores its part of
   * the last hide-all. Returns false when there is nothing to toggle.
   */
  toggleAll(scope = "all") {
    const items = this.toggleItems(scope);
    if (!items.length) return false;
    const group = this.viewer.state.group;
    const key = `${group}::${scope}`;
    const visible = new Set(items.filter((it) => this.itemVisible(it)).map((it) => it.id));
    const remembered = this.snapshots.get(key) || this.snapshots.get(`${group}::all`) || null;
    const plan = toggleAllPlan(items.map((it) => it.id), visible, remembered);
    if (plan.snapshot) this.snapshots.set(key, plan.snapshot);
    else this.snapshots.delete(key);
    for (const it of items) this.setItemVisible(it, plan.visible.get(it.id));
    this.sync();
    return true;
  }

  /** A section's show/hide-all button; `sync` keeps its icon and label current. */
  sectionEye(scope) {
    const b = iconButton("eye", "", () => this.toggleAll(scope), { cls: "ov-section-eye", shortcut: scope === "all" ? "T" : "" });
    b.dataset.scope = scope;
    this.sectionEyes.push(b);
    return b;
  }

  syncSectionEyes() {
    for (const b of this.sectionEyes) {
      const scope = b.dataset.scope;
      const on = this.toggleItems(scope).some((it) => this.itemVisible(it));
      const noun = scope === "volumes" ? "all volumes" : scope === "traces" ? "all data layers" : this.manifest.volumes?.length ? "all layers, volumes included" : "all layers";
      const label = `${on ? "Hide" : "Show"} ${noun}`;
      b.replaceChildren(icon(on ? "eye" : "eyeOff"));
      b.setAttribute("aria-label", label);
      b.setAttribute("aria-pressed", String(on));
      b.dataset.tip = scope === "all" ? `${label}  T` : label;
    }
  }

  /**
   * Show one layer alone: a trace key, or `{ volume: stateKey }`. Every
   * other layer, volumes included, goes off (the classic double-click);
   * soloing it again restores the mix that was on before.
   */
  solo(key) {
    const id = typeof key === "object" && key ? `volume:${key.volume}` : `trace:${key}`;
    const items = this.toggleItems("all");
    const group = this.viewer.state.group;
    const snapKey = `${group}::solo`;
    const visible = new Set(items.filter((it) => this.itemVisible(it)).map((it) => it.id));
    const plan = soloPlan(items.map((it) => it.id), visible, id, this.snapshots.get(snapKey));
    if (plan.snapshot) this.snapshots.set(snapKey, plan.snapshot);
    else this.snapshots.delete(snapKey);
    for (const it of items) this.setItemVisible(it, plan.visible.get(it.id));
    this.sync();
  }

  /**
   * The number keys: 1–9 toggle the first nine layers (data, then volumes),
   * Shift shows that one alone, volumes included, like the classic viewer.
   */
  toggleItemByIndex(index, { solo = false } = {}) {
    const items = this.toggleItems("all");
    const target = items[index];
    if (!target) return false;
    if (solo) this.solo(target.kind === "volume" ? { volume: target.key } : target.key);
    else this.setItemVisible(target, !this.itemVisible(target));
    return true;
  }

  // ---------------------------------------------------------------- traces

  traceRow(trace, i) {
    const v = this.viewer;
    const st = () => v.state.traces[trace.key] || {};
    const swatch = h("span", { class: "ov-swatch" });
    const swBtn = h("button", { class: "ov-swatch-btn", type: "button", "aria-label": `Toggle ${trace.name}` }, swatch);
    const name = h("button", { class: "ov-row-name", type: "button", title: `${trace.name}\nClick to show/hide · double-click to solo` }, trace.name);
    const count = h("span", { class: "ov-row-count" }, trace.points ? fmtInt(trace.points.count) : "");
    const eye = iconButton("eye", `Show or hide ${trace.name}`, () => toggle(), { cls: "ov-eye", shortcut: i < 9 ? String(i + 1) : "" });
    const chev = iconButton("chevron", `Edit ${trace.name}`, () => this.expand(trace.key), { cls: "ov-chev" });
    const row = h("div", { class: "ov-row", role: "listitem" }, swBtn, name, count, h("div", { class: "ov-row-actions" }, eye, chev));
    const toggle = () => v.setTraceStyle(trace.key, { visible: !st().visible });
    swBtn.addEventListener("click", toggle);
    let clickTimer = 0;
    name.addEventListener("click", () => {
      clearTimeout(clickTimer);
      clickTimer = setTimeout(toggle, 200);
    });
    name.addEventListener("dblclick", () => {
      clearTimeout(clickTimer);
      this.solo(trace.key);
    });
    const editorHost = h("div");
    this.rows.set(trace.key, { row, swatch, eye, trace, editorHost, kind: "trace" });
    return [row, editorHost];
  }

  expand(key) {
    this.expanded = this.expanded === key ? null : key;
    for (const [k, r] of this.rows) {
      const open = k === this.expanded;
      r.row.dataset.expanded = String(open);
      clear(r.editorHost);
      if (open) r.editorHost.append(r.kind === "volume" ? this.volumeEditor(r.spec) : this.traceEditor(r.trace));
    }
  }

  traceEditor(trace) {
    const v = this.viewer;
    const st = v.state.traces[trace.key] || {};
    const wrap = h("div", { class: "ov-editor" });
    const cb = trace.colorBy;
    const mode = st.colorMode || cb?.defaultMode || "fixed";
    const colorField = h("div", { class: "ov-field" },
      h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, "Colour")),
      colorPicker({ value: st.color || trace.color, onChange: (c) => v.setTraceStyle(trace.key, { color: c, colorMode: "fixed" }) }),
    );
    const cmapHost = h("div", { class: "ov-field" });
    const renderCmap = () => {
      clear(cmapHost);
      const s = v.state.traces[trace.key] || {};
      if ((s.colorMode || mode) !== "by_value" || !cb) { cmapHost.hidden = true; return; }
      cmapHost.hidden = false;
      const names = cb.colormaps?.length ? cb.colormaps : [cb.colormap];
      const cur = s.colormap || cb.colormap;
      const ramp = h("div", { class: "ov-ramp" });
      const lut = this.lutFor.get(cur);
      if (lut) ramp.style.background = lutGradient(lut);
      cmapHost.append(
        select({ label: `${cb.label || "Value"} colormap`, value: cur, options: names.map((n) => ({ value: n, label: this.ui.cmapLabel(n) })), onChange: (n) => { v.setTraceStyle(trace.key, { colormap: n }); renderCmap(); } }),
        ramp,
        h("div", { class: "ov-field-head" }, h("span", { class: "ov-faint ov-num" }, fmt(s.cmin ?? cb.cmin)), h("span", { class: "ov-faint ov-num" }, fmt(s.cmax ?? cb.cmax))),
      );
    };
    if (cb) {
      wrap.append(miniSeg({
        label: "Colour by",
        value: mode,
        options: [{ value: "fixed", label: "Solid" }, { value: "by_value", label: cb.label || "Value" }],
        onChange: (m) => {
          v.setTraceStyle(trace.key, { colorMode: m, ...(m === "by_value" ? {} : {}) });
          colorField.hidden = m === "by_value";
          renderCmap();
        },
      }));
    }
    colorField.hidden = mode === "by_value";
    wrap.append(colorField, cmapHost);
    renderCmap();
    if (trace.points || trace.lines || trace.labels) {
      wrap.append(slider({
        label: "Opacity", min: 0, max: 1, step: 0.01, value: st.opacity ?? trace.opacity,
        format: (x) => `${Math.round(x * 100)}%`,
        onInput: (x) => v.setTraceStyle(trace.key, { opacity: x }),
      }));
      wrap.append(slider({
        label: trace.labels ? "Text size" : trace.lines ? "Line width" : "Point size",
        min: 0.1, max: 5, value: st.sizeScale ?? 1, scale: "log",
        format: (x) => `${x.toFixed(2)}×`,
        onInput: (x) => v.setTraceStyle(trace.key, { sizeScale: x }),
      }));
    }
    // Member-star motion arrows in Sky view (classic per-trace setting).
    if (trace.points?.members && this.manifest.sky?.members) {
      const on = !!st.memberArrows;
      const box = h("div", { style: { display: on ? "grid" : "none", gap: "10px" } },
        slider({ label: "Arrow length", min: 0.1, max: 5, value: st.memberArrowLength ?? 1, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setTraceStyle(trace.key, { memberArrowLength: x }) }),
        slider({ label: "Arrow width", min: 0.2, max: 4, value: st.memberArrowWidth ?? 1, scale: "log", format: (x) => `${x.toFixed(2)}×`, onInput: (x) => v.setTraceStyle(trace.key, { memberArrowWidth: x }) }));
      wrap.append(toggle({
        label: "Motion arrows in Sky", hint: "Each member star's motion relative to its cluster, as 3 Myr of travel",
        checked: on,
        onChange: (x) => { v.setTraceStyle(trace.key, { memberArrows: x }); box.style.display = x ? "grid" : "none"; },
      }), box);
    }
    const actions = h("div", { class: "ov-field-row" });
    if (trace.points) {
      actions.append(h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.ui.frameTrace(trace.key) }, icon("target"), "Frame"));
      // Follow the layer's median through time (classic "focus group").
      if (trace.points.count > 1) {
        const following = v.anchor.kind === "layer" && v.anchor.trace === trace.key;
        actions.append(h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-pressed": String(following), title: "The camera moves with this group's median through time", onclick: () => { this.ui.followLayer(trace.key); this.expand(trace.key); this.expand(trace.key); } }, icon("follow"), "Follow"));
      }
    }
    actions.append(h("button", { class: "ov-btn ov-btn--sm ov-btn--ghost", type: "button", onclick: () => {
      const vis = v.state.traces[trace.key]?.visible;
      v.state.traces[trace.key] = { visible: vis, inGroup: true };
      v.setTraceStyle(trace.key, {});
      this.expand(trace.key); this.expand(trace.key);
    } }, "Reset"));
    wrap.append(actions);
    return wrap;
  }

  // ---------------------------------------------------------------- volumes

  volumeRow(spec, index = 99) {
    const v = this.viewer;
    const key = spec.stateKey;
    const swatch = h("span", { class: "ov-swatch ov-swatch--volume" });
    const swBtn = h("button", { class: "ov-swatch-btn", type: "button", "aria-label": `Toggle ${spec.name}` }, swatch);
    const name = h("button", { class: "ov-row-name", type: "button", title: `${spec.name}${spec.unit ? ` (${spec.unit})` : ""}\nClick to show/hide · double-click to solo` }, spec.name);
    const status = h("span", { class: "ov-row-count" });
    const eye = iconButton("eye", `Show or hide ${spec.name}`, () => toggleVol(), { cls: "ov-eye", shortcut: index < 9 ? String(index + 1) : "" });
    const chev = iconButton("chevron", `Edit ${spec.name}`, () => this.expand(key), { cls: "ov-chev" });
    const row = h("div", { class: "ov-row" }, swBtn, name, status, h("div", { class: "ov-row-actions" }, eye, chev));
    const toggleVol = () => v.setVolumeState(key, { visible: !v.state.volumes[key]?.visible });
    swBtn.addEventListener("click", toggleVol);
    let clickTimer = 0;
    name.addEventListener("click", () => {
      clearTimeout(clickTimer);
      clickTimer = setTimeout(toggleVol, 200);
    });
    name.addEventListener("dblclick", () => {
      clearTimeout(clickTimer);
      this.solo({ volume: key });
    });
    const editorHost = h("div");
    this.rows.set(key, { row, swatch, eye, spec, editorHost, kind: "volume", status });
    return [row, editorHost];
  }

  volumeEditor(spec) {
    const v = this.viewer;
    const key = spec.stateKey;
    const vs = clampVolumeState(spec, v.state.volumes[key]);
    const [dmin, dmax] = spec.dataRange;
    const wrap = h("div", { class: "ov-editor" });
    const ramp = h("div", { class: "ov-ramp" });
    const paint = (name) => {
      const lut = this.lutFor.get(name);
      if (lut) ramp.style.background = lutGradient(lut);
    };
    paint(vs.colormap);
    const names = spec.colormaps.length ? spec.colormaps : [vs.colormap];
    const proc = spec.procedural;
    if (proc?.windowOptions?.length) {
      // Event densities: how many Myr of events each moment counts.
      const cur = Number(vs.kdeTimeWindowMyr) || proc.windowMyr;
      const opts = [...new Set([...proc.windowOptions, cur])].sort((a, b) => a - b);
      wrap.append(select({
        label: "Time window", value: String(cur),
        options: opts.map((w) => ({ value: String(w), label: `${w} Myr trailing` })),
        onChange: (w) => v.setVolumeState(key, { kdeTimeWindowMyr: Number(w) }),
      }));
    }
    if (spec.note) wrap.append(h("div", { class: "ov-note" }, spec.note));
    wrap.append(
      select({ label: "Colormap", value: vs.colormap, options: names.map((n) => ({ value: n, label: this.ui.cmapLabel(n) })), onChange: (n) => { v.setVolumeState(key, { colormap: n }); paint(n); } }),
      ramp,
    );
    // A fixed display encoding (pre-scaled data, a colour field) has no
    // meaningful window or stretch to edit.
    if (!spec.fixedDisplay) wrap.append(
      miniSeg({
        label: "Stretch", value: vs.stretch === "log10" ? "log10" : vs.stretch,
        options: [{ value: "linear", label: "Linear" }, { value: "log10", label: "Log" }, { value: "asinh", label: "Asinh" }],
        onChange: (s) => v.setVolumeState(key, { stretch: s }),
      }),
    );
    const unit = spec.unit ? ` ${spec.unit}` : "";
    const fmtV = (x) => `${formatSci(x)}${unit}`;
    const span = dmax - dmin || 1;
    const low = slider({ label: "Window low", min: dmin, max: dmax, step: span / 1000, value: vs.vmin, format: fmtV, onInput: (x) => { v.setVolumeState(key, { vmin: x }); lowNum.set(+x.toPrecision(4)); } });
    const high = slider({ label: "Window high", min: dmin, max: dmax, step: span / 1000, value: vs.vmax, format: fmtV, onInput: (x) => { v.setVolumeState(key, { vmax: x }); highNum.set(+x.toPrecision(4)); } });
    // Exact values, as the classic number fields allowed (e.g. vmax 0.07).
    const exact = (field, label, apply) => numberInput({ label, value: +Number(vs[field]).toPrecision(4), step: niceStep(span), onChange: (x) => { const c = clampVolumeState(spec, { ...v.state.volumes[key], [field]: x }); v.setVolumeState(key, { [field]: c[field] }); apply(c[field]); } });
    const lowNum = exact("vmin", `Low${unit}`, (x) => low.set(x));
    const highNum = exact("vmax", `High${unit}`, (x) => high.set(x));
    if (!spec.fixedDisplay) wrap.append(low, high, h("div", { class: "ov-field-row" }, lowNum, highNum));
    wrap.append(
      slider({ label: "Opacity", min: 0, max: 1, value: vs.opacity, format: (x) => `${Math.round(x * 100)}%`, onInput: (x) => v.setVolumeState(key, { opacity: x }) }),
      slider({ label: "Density gain", min: 1, max: 200, value: vs.alphaCoef, scale: "log", format: (x) => x.toFixed(0), onInput: (x) => v.setVolumeState(key, { alphaCoef: x }) }),
      slider({ label: "Samples", min: 24, max: 768, step: 8, value: vs.steps, format: (x) => String(Math.round(x)), onInput: (x) => v.setVolumeState(key, { steps: Math.round(x) }) }),
    );
    if (spec.supportsShowAllTimes && spec.timeMyr == null) {
      wrap.append(toggle({
        label: "Show at every time", hint: spec.coRotate ? "Co-rotates with the Galactic frame" : "",
        checked: vs.showAllTimes, onChange: (x) => v.setVolumeState(key, { showAllTimes: x }),
      }));
    }
    wrap.append(h("div", { class: "ov-note" }, `Data range ${formatSci(dmin)}–${formatSci(dmax)}${unit} · ${spec.dims.join("×")} voxels`));
    return wrap;
  }

  /** "t = 0 only" when a visible volume is not drawn at this time (classic render-mode note). */
  volumeTimeNote(spec, vs) {
    const v = this.viewer;
    if (!vs.visible) return "";
    const drawn = (v.volumes.params?.volumeDraws || []).some((d) => d.key === spec.key || v.volumeSpecs?.find((s) => s.key === d.key)?.stateKey === spec.stateKey);
    if (drawn) return "";
    const sky = v.state.view.mode === "sky";
    return h("span", { class: "ov-faint", title: sky ? "Sky view shows volumes at the present day only (press 0)" : "Shown at the present day only (press 0), or turn on Show at every time" }, "t = 0 only");
  }

  // ---------------------------------------------------------------- reference

  referenceSection() {
    const v = this.viewer;
    const g = v.state.global;
    const sec = h("div", { class: "ov-section" }, h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Guides")));
    const inner = h("div", { style: { padding: "2px 14px 4px", display: "grid", gap: "4px" } });
    this.guideToggles = new Map();
    for (const t of this.guideTraces()) {
      const el = toggle({ label: t.name, checked: v.state.traces[t.key]?.visible !== false, onChange: (x) => v.setTraceStyle(t.key, { visible: x }) });
      this.guideToggles.set(t.key, el);
      inner.append(el);
    }
    const hasRef = this.manifest.traces.some((t) => !t.showInLegend);
    if (hasRef) inner.append(this.gridToggle = toggle({ label: "Galactic grid", hint: "Radius circles, quadrants and labels", checked: g.grid, onChange: (x) => v.setGlobal({ grid: x }) }));
    for (const img of this.manifest.images || []) {
      inner.append(toggle({
        label: "Milky Way image", checked: v.state.images[img.key]?.visible !== false,
        onChange: (x) => { v.state.images[img.key].visible = x; v.renderer.invalidate(); },
      }));
    }
    inner.append(this.labelsToggle = toggle({ label: "Text labels", checked: g.labels !== false, onChange: (x) => v.setGlobal({ labels: x }) }));
    sec.append(inner);
    return sec;
  }

  // ---------------------------------------------------------------- sync

  sync() {
    const v = this.viewer;
    for (const [key, r] of this.rows) {
      if (r.kind === "trace") {
        const st = v.state.traces[key] || {};
        const visible = !!st.visible;
        r.row.dataset.visible = String(visible);
        r.eye.replaceChildren(icon(visible ? "eye" : "eyeOff"));
        const cb = r.trace.colorBy;
        const mode = st.colorMode || cb?.defaultMode || "fixed";
        if (mode === "by_value" && cb) {
          const lut = this.lutFor.get(st.colormap || cb.colormap);
          r.swatch.className = "ov-swatch ov-swatch--ramp";
          r.swatch.style.background = lut ? lutGradient(lut, "to right", 5) : "";
          r.swatch.style.color = "";
        } else {
          r.swatch.className = "ov-swatch";
          r.swatch.style.background = "";
          r.swatch.style.color = st.color || r.trace.legendColor || r.trace.color;
        }
      } else if (r.kind === "volume") {
        const vs = v.state.volumes[key] || {};
        const visible = !!vs.visible;
        r.row.dataset.visible = String(visible);
        r.eye.replaceChildren(icon(visible ? "eye" : "eyeOff"));
        const lut = this.lutFor.get(vs.colormap);
        r.swatch.style.background = lut ? lutGradient(lut, "135deg", 4) : r.spec.legendColor;
        const loaded = v.volumes.volumes.has(r.spec.key);
        r.status.replaceChildren(loaded ? this.volumeTimeNote(r.spec, vs) : h("span", { class: "ov-spinner", title: "Loading…" }));
      }
    }
    this.syncSectionEyes();
    for (const [key, el] of this.guideToggles || []) el.set(v.state.traces[key]?.visible !== false);
    this.gridToggle?.set(v.state.global.grid);
    this.labelsToggle?.set(v.state.global.labels !== false);
  }
}

function niceStep(span) {
  const s = Math.abs(span) / 200;
  if (!(s > 0)) return "any";
  const p = Math.pow(10, Math.floor(Math.log10(s)));
  return String(+(Math.ceil(s / p) * p).toPrecision(2));
}

function uniqueVolumes(specs) {
  const seen = new Set();
  const out = [];
  for (const s of specs) {
    if (seen.has(s.stateKey)) continue;
    seen.add(s.stateKey);
    out.push(s);
  }
  return out;
}

function formatSci(x) {
  if (!Number.isFinite(x)) return "—";
  const a = Math.abs(x);
  if (a === 0) return "0";
  if (a >= 1e4 || a < 1e-3) return x.toExponential(2).replace("e", "×10^").replace("-", "−");
  if (a >= 100) return x.toFixed(0);
  if (a >= 10) return x.toFixed(1);
  if (a >= 1) return x.toFixed(2);
  return x.toPrecision(3);
}

