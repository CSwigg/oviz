// Layers panel: group presets, data traces, volumes and reference guides.

import { h, icon, iconButton, clear, fmtInt, fmt } from "./dom.js";
import { slider, toggle, select, miniSeg, colorPicker } from "./controls.js";
import { lutGradient } from "../core/color.js";
import { applyGroupDefaults, groupMode, clampVolumeState } from "../app/state.js";

export class LayersPanel {
  constructor(ui) {
    this.ui = ui;
    this.viewer = ui.viewer;
    this.manifest = ui.viewer.manifest;
    this.expanded = null;
    this.rows = new Map();
    this.el = h("section", { class: "ov-panel ov-layers ov-glass ov-chrome", "aria-label": "Layers", "data-open": "true" });
    this.head = h("div", { class: "ov-panel-head" },
      h("span", { class: "ov-panel-title" }, "Layers"),
      iconButton("close", "Hide layers", () => ui.setLayersOpen(false), { shortcut: "L" }),
    );
    this.body = h("div", { class: "ov-panel-body" });
    this.el.append(this.head, this.body);
    this.render();
    const v = this.viewer;
    v.on("style", () => this.sync());
    v.on("volume", () => this.sync());
    v.on("global", () => this.sync());
    v.on("volume-ready", () => this.sync());
    v.on("viewmode", () => this.render());
    v.on("state-applied", () => this.render());
  }

  get lutFor() {
    return this.ui.luts;
  }

  // ---------------------------------------------------------------- render

  render() {
    clear(this.body);
    this.rows.clear();
    const v = this.viewer;
    const m = this.manifest;
    const groups = m.groups?.order || [];
    if (groups.length > 1) {
      const sel = h("select", { "aria-label": "Layer group" }, groups.map((g) => h("option", { value: g, selected: g === v.state.group }, g)));
      sel.addEventListener("change", () => this.setGroup(sel.value));
      this.body.append(h("div", { class: "ov-group" }, h("label", null, "View"), sel));
    }
    const traces = m.traces.filter((t) => t.showInLegend && groupMode(m, v.state.group, t.key) !== false);
    if (traces.length) {
      const sec = h("div", { class: "ov-section" },
        h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Data"),
          iconButton("eye", "Show or hide all data layers", () => this.toggleAll(traces), { shortcut: "T" })),
      );
      traces.forEach((t, i) => sec.append(...this.traceRow(t, i)));
      this.body.append(sec);
    }
    const filterSection = this.ui.filter?.section();
    if (filterSection) this.body.append(filterSection);
    const volKeys = uniqueVolumes(m.volumes || []);
    if (volKeys.length) {
      const sec = h("div", { class: "ov-section" }, h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Volumes")));
      for (const spec of volKeys) sec.append(...this.volumeRow(spec));
      this.body.append(sec);
    }
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
    v.renderer.invalidate();
    v.invalidatePick();
    v.emit("style", { key: "*" });
    this.render();
  }

  toggleAll(traces) {
    const v = this.viewer;
    const anyOn = traces.some((t) => v.state.traces[t.key]?.visible);
    for (const t of traces) v.setTraceStyle(t.key, { visible: !anyOn });
  }

  solo(key) {
    const v = this.viewer;
    const traces = this.manifest.traces.filter((t) => t.showInLegend && groupMode(this.manifest, v.state.group, t.key) !== false);
    const onlyThis = traces.every((t) => (t.key === key) === !!v.state.traces[t.key]?.visible);
    for (const t of traces) v.setTraceStyle(t.key, { visible: onlyThis ? true : t.key === key });
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
    const actions = h("div", { class: "ov-field-row" });
    if (trace.points) {
      actions.append(h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.ui.frameTrace(trace.key) }, icon("target"), "Frame"));
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

  volumeRow(spec) {
    const v = this.viewer;
    const key = spec.stateKey;
    const swatch = h("span", { class: "ov-swatch ov-swatch--volume" });
    const swBtn = h("button", { class: "ov-swatch-btn", type: "button", "aria-label": `Toggle ${spec.name}` }, swatch);
    const name = h("button", { class: "ov-row-name", type: "button", title: `${spec.name}${spec.unit ? ` (${spec.unit})` : ""}` }, spec.name);
    const status = h("span", { class: "ov-row-count" });
    const eye = iconButton("eye", `Show or hide ${spec.name}`, () => toggleVol(), { cls: "ov-eye" });
    const chev = iconButton("chevron", `Edit ${spec.name}`, () => this.expand(key), { cls: "ov-chev" });
    const row = h("div", { class: "ov-row" }, swBtn, name, status, h("div", { class: "ov-row-actions" }, eye, chev));
    const toggleVol = () => v.setVolumeState(key, { visible: !v.state.volumes[key]?.visible });
    swBtn.addEventListener("click", toggleVol);
    name.addEventListener("click", toggleVol);
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
    wrap.append(
      select({ label: "Colormap", value: vs.colormap, options: names.map((n) => ({ value: n, label: this.ui.cmapLabel(n) })), onChange: (n) => { v.setVolumeState(key, { colormap: n }); paint(n); } }),
      ramp,
      miniSeg({
        label: "Stretch", value: vs.stretch === "log10" ? "log10" : vs.stretch,
        options: [{ value: "linear", label: "Linear" }, { value: "log10", label: "Log" }, { value: "asinh", label: "Asinh" }],
        onChange: (s) => v.setVolumeState(key, { stretch: s }),
      }),
    );
    const unit = spec.unit ? ` ${spec.unit}` : "";
    const fmtV = (x) => `${formatSci(x)}${unit}`;
    const span = dmax - dmin || 1;
    wrap.append(
      slider({ label: "Window low", min: dmin, max: dmax, step: span / 1000, value: vs.vmin, format: fmtV, onInput: (x) => v.setVolumeState(key, { vmin: x }) }),
      slider({ label: "Window high", min: dmin, max: dmax, step: span / 1000, value: vs.vmax, format: fmtV, onInput: (x) => v.setVolumeState(key, { vmax: x }) }),
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

  // ---------------------------------------------------------------- reference

  referenceSection() {
    const v = this.viewer;
    const g = v.state.global;
    const sec = h("div", { class: "ov-section" }, h("div", { class: "ov-section-head" }, h("span", { class: "ov-grow" }, "Guides")));
    const inner = h("div", { style: { padding: "2px 14px 4px", display: "grid", gap: "4px" } });
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
        r.status.replaceChildren(loaded ? "" : h("span", { class: "ov-spinner", title: "Loading…" }));
      }
    }
    this.gridToggle?.set(v.state.global.grid);
    this.labelsToggle?.set(v.state.global.labels !== false);
  }
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

export function formatSci(x) {
  if (!Number.isFinite(x)) return "—";
  const a = Math.abs(x);
  if (a === 0) return "0";
  if (a >= 1e4 || a < 1e-3) return x.toExponential(2).replace("e", "×10^").replace("-", "−");
  if (a >= 100) return x.toFixed(0);
  if (a >= 10) return x.toFixed(1);
  if (a >= 1) return x.toFixed(2);
  return x.toPrecision(3);
}

