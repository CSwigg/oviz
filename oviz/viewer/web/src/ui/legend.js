// Compact legend: the figure's data layers as coloured dots and names, the
// way a printed figure has a key. Click a name to show or hide the layer,
// double-click to show it alone. Volumes are listed with the data, as in
// the layers panel. On touch screens the rows grow to finger size, a tap
// toggles at once, a long press shows the layer alone, and a first row
// hides or shows everything (the T key). The full layers panel stays one
// click away.

import { h, icon } from "./dom.js";
import { lutGradient } from "../core/color.js";

const MAX_ITEMS = 8;
const LONG_PRESS_MS = 480;

export class CompactLegend {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.rows = [];
    this.group = null;
    this.el = h("div", { class: "ov-legend", role: "group", "aria-label": "Legend" });
    v.on("style", () => this.sync());
    v.on("volume", () => this.sync());
    v.on("state-applied", () => this.render());
    // Present-day maps fade out away from t = 0: the presenting key drops them.
    v.on("time", () => this.syncDrawn());
    this.render();
  }

  /** Mark volume rows not drawn at this time (the presenting key hides them). */
  syncDrawn() {
    const v = this.viewer;
    const rows = this.rows.filter((r) => r.it.kind === "volume");
    if (!rows.length) return;
    const drawn = new Set();
    try {
      const stateKey = new Map(v.volumeSpecs.map((s) => [s.key, s.stateKey]));
      for (const d of v._volumeParams(v.timeline.frame, v.timeline.time).volumeDraws || []) {
        if (d.fade > 0.002) drawn.add(stateKey.get(d.key));
      }
    } catch (_) { return; }
    for (const r of rows) r.b.dataset.drawn = String(drawn.has(r.it.key));
  }

  /** Data traces, then volumes (the layers panel's order), at most MAX_ITEMS. */
  items() {
    let all = this.ui.layers.toggleItems("all");
    // Presenting, the key lists only what is in view.
    if (this.ui.root.dataset.presenting === "true") all = all.filter((it) => this.ui.layers.itemVisible(it));
    const volumes = all.filter((it) => it.kind === "volume").slice(0, MAX_ITEMS - 1);
    const traces = all.filter((it) => it.kind === "trace");
    const room = Math.max(1, MAX_ITEMS - volumes.length);
    return { shown: [...traces.slice(0, room), ...volumes], extra: Math.max(0, traces.length - room) };
  }

  render() {
    const v = this.viewer;
    this.group = v.state.group;
    const { shown, extra } = this.items();
    const layers = this.ui.layers;
    this.rows = shown.map((it) => {
      const dot = h("span", { class: "ov-legend-dot", "data-kind": it.kind });
      const b = h("button", {
        class: "ov-legend-item", type: "button",
        "data-tip": `${it.name}\nClick to show or hide · double-click to show only this`,
      }, dot, h("span", { class: "ov-legend-name" }, it.name));
      const flip = () => layers.setItemVisible(it, !layers.itemVisible(it));
      const solo = () => layers.solo(it.kind === "volume" ? { volume: it.key } : it.key);
      let timer = 0;
      const touched = bindTouch(b, {
        tap: flip,
        press: () => {
          solo();
          const alone = layers.itemVisible(it) && layers.toggleItems("all").every((x) => x.id === it.id || !layers.itemVisible(x));
          this.ui.toast?.(alone ? `Only ${it.name} · press and hold again to bring the rest back` : "The other layers are back", { ms: 2200 });
        },
      });
      b.addEventListener("click", () => {
        if (touched()) return; // the tap already toggled it
        clearTimeout(timer);
        timer = setTimeout(flip, 200);
      });
      b.addEventListener("dblclick", () => {
        if (touched()) return;
        clearTimeout(timer);
        solo();
      });
      return { it, b, dot };
    });
    // Touch screens have no T key: this row hides or shows every layer.
    this.allIcon = h("span", { class: "ov-legend-all-icon" });
    this.allLabel = h("span", { class: "ov-legend-name" });
    this.allBtn = h("button", { class: "ov-legend-item ov-legend-all", type: "button" }, this.allIcon, this.allLabel);
    const allTouched = bindTouch(this.allBtn, { tap: () => layers.toggleAll("all") });
    this.allBtn.addEventListener("click", () => { if (!allTouched()) layers.toggleAll("all"); });
    const more = extra > 0
      ? h("button", { class: "ov-legend-more", type: "button", onclick: () => this.ui.setLayersOpen(true) }, `+${extra} more`)
      : null;
    this.el.replaceChildren(...(shown.length > 1 ? [this.allBtn] : []), ...this.rows.map((r) => r.b), ...(more ? [more] : []));
    this.el.hidden = !shown.length;
    this.sync();
  }

  sync() {
    const v = this.viewer;
    if (this.group !== v.state.group) { this.render(); return; }
    if (this.allBtn) {
      const anyOn = this.ui.layers.toggleItems("all").some((it) => this.ui.layers.itemVisible(it));
      this.allIcon.replaceChildren(icon(anyOn ? "eyeOff" : "eye"));
      this.allLabel.textContent = anyOn ? "Hide all" : "Show all";
      this.allBtn.setAttribute("aria-label", anyOn ? "Hide all layers" : "Show all layers");
    }
    for (const r of this.rows) {
      const visible = this.ui.layers.itemVisible(r.it);
      r.b.dataset.visible = String(visible);
      r.b.setAttribute("aria-pressed", String(visible));
      if (r.it.kind === "volume") {
        const vs = v.state.volumes[r.it.key] || {};
        const lut = this.ui.luts.get(vs.colormap);
        r.dot.style.background = lut ? lutGradient(lut, "135deg", 4) : "";
        continue;
      }
      const t = v.traceByKey.get(r.it.key);
      const st = v.state.traces[r.it.key] || {};
      const cb = t?.colorBy;
      const mode = st.colorMode || cb?.defaultMode || "fixed";
      const lut = mode === "by_value" && cb ? this.ui.luts.get(st.colormap || cb.colormap) : null;
      if (lut) {
        r.dot.style.background = lutGradient(lut, "to right", 5);
        r.dot.dataset.ramp = "true";
      } else {
        r.dot.style.background = st.color || t?.legendColor || t?.color;
        delete r.dot.dataset.ramp;
      }
    }
    this.syncDrawn();
  }
}

/**
 * Touch taps and long presses on a button. A tap runs `tap` at once (no
 * wait for a double-click); holding still for LONG_PRESS_MS runs `press`.
 * Moving the finger cancels both. Returns a function that says whether
 * the button was just touched, so its click (and dblclick) handlers can
 * leave touch input to these.
 */
function bindTouch(el, { tap, press }) {
  let down = null;
  let touchedAt = -Infinity;
  const cancel = () => { if (down) clearTimeout(down.timer); down = null; };
  el.addEventListener("pointerdown", (e) => {
    if (e.pointerType !== "touch") return;
    cancel();
    down = { x: e.clientX, y: e.clientY, long: false, timer: 0 };
    if (press) down.timer = setTimeout(() => { if (down) { down.long = true; press(); } }, LONG_PRESS_MS);
  });
  el.addEventListener("pointermove", (e) => {
    if (down && Math.hypot(e.clientX - down.x, e.clientY - down.y) > 12) cancel();
  });
  el.addEventListener("pointerup", (e) => {
    if (e.pointerType !== "touch") return;
    touchedAt = performance.now();
    if (!down) return;
    clearTimeout(down.timer);
    const long = down.long;
    down = null;
    if (!long) tap();
  });
  el.addEventListener("pointercancel", cancel);
  // No context menu (or text callout) after a long press.
  el.addEventListener("contextmenu", (e) => { if (down || performance.now() - touchedAt < 800) e.preventDefault(); });
  return () => performance.now() - touchedAt < 800;
}
