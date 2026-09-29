// Compact legend: the figure's data layers as coloured dots and names, the
// way a printed figure has a key. Click a name to show or hide the layer,
// double-click to show it alone. The full layers panel stays one click away.

import { h } from "./dom.js";
import { groupMode } from "../app/state.js";
import { lutGradient } from "../core/color.js";

const MAX_ITEMS = 8;

export class CompactLegend {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.rows = [];
    this.group = null;
    this.el = h("div", { class: "ov-legend", role: "group", "aria-label": "Legend" });
    v.on("style", () => this.sync());
    v.on("state-applied", () => this.render());
    this.render();
  }

  traces() {
    const v = this.viewer;
    const m = v.manifest;
    return m.traces.filter((t) => t.showInLegend && groupMode(m, v.state.group, t.key) !== false);
  }

  render() {
    const v = this.viewer;
    this.group = v.state.group;
    const all = this.traces();
    this.rows = all.slice(0, MAX_ITEMS).map((t) => {
      const dot = h("span", { class: "ov-legend-dot" });
      const b = h("button", {
        class: "ov-legend-item", type: "button",
        "data-tip": `${t.name}\nClick to show or hide · double-click to show only this`,
      }, dot, h("span", { class: "ov-legend-name" }, t.name));
      let timer = 0;
      b.addEventListener("click", () => {
        clearTimeout(timer);
        timer = setTimeout(() => v.setTraceStyle(t.key, { visible: !v.state.traces[t.key]?.visible }), 200);
      });
      b.addEventListener("dblclick", () => {
        clearTimeout(timer);
        this.ui.layers.solo(t.key);
      });
      return { t, b, dot };
    });
    const extra = all.length - this.rows.length;
    const more = extra > 0
      ? h("button", { class: "ov-legend-more", type: "button", onclick: () => this.ui.setLayersOpen(true) }, `+${extra} more`)
      : null;
    this.el.replaceChildren(...this.rows.map((r) => r.b), ...(more ? [more] : []));
    this.el.hidden = !all.length;
    this.sync();
  }

  sync() {
    const v = this.viewer;
    if (this.group !== v.state.group) { this.render(); return; }
    for (const r of this.rows) {
      const st = v.state.traces[r.t.key] || {};
      r.b.dataset.visible = String(!!st.visible);
      r.b.setAttribute("aria-pressed", String(!!st.visible));
      const cb = r.t.colorBy;
      const mode = st.colorMode || cb?.defaultMode || "fixed";
      const lut = mode === "by_value" && cb ? this.ui.luts.get(st.colormap || cb.colormap) : null;
      if (lut) {
        r.dot.style.background = lutGradient(lut, "to right", 5);
        r.dot.dataset.ramp = "true";
      } else {
        r.dot.style.background = st.color || r.t.legendColor || r.t.color;
        delete r.dot.dataset.ramp;
      }
    }
  }
}
