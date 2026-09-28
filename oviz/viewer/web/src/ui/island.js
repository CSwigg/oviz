// Island skin: the search capsule at the top becomes a live "island" that
// stretches to show what is happening (the time while playing, the object
// you picked, short notifications) and flows open into search. In other
// skins it stays a plain search button.

import { h } from "./dom.js";
import { formatDistance } from "../core/math.js";

export class Island {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.btn = ui.search;
    this.live = h("span", { class: "ov-island-live", "aria-live": "polite" });
    this.btn.append(this.live);
    this.flashTimer = 0;
    this.flashing = null;
    this.state = "";
    v.on("time", () => this.refresh());
    v.on("select", (hit) => {
      this.refresh();
      // Names arrive with the object metadata; redraw once they do.
      if (hit && hit.kind !== "member") v.objectMeta(hit.trace).then(() => { this.key = ""; this.refresh(); });
    });
    v.on("skin", () => this.refresh());
    this.ui.dock?.tl?.on?.(() => this.refresh());
  }

  get active() {
    return document.documentElement.dataset.ovizSkin === "island";
  }

  /** Show a short message in the island (returns false outside the skin). */
  flash(text, { iconName = null, ms = 1800 } = {}) {
    if (!this.active) return false;
    this.flashing = { text, iconName };
    clearTimeout(this.flashTimer);
    this.flashTimer = setTimeout(() => { this.flashing = null; this.refresh(); }, ms);
    this.refresh();
    return true;
  }

  refresh() {
    if (!this.active) {
      if (this.state) { this.state = ""; this.key = ""; this.btn.dataset.island = ""; this.live.replaceChildren(); this.btn.style.width = ""; }
      return;
    }
    const v = this.viewer;
    const tl = v.timeline;
    const sel = this.ui.selection;
    const f = tl.count > 1 ? tl.frame / (tl.count - 1) : 1;
    const state = this.flashing ? "flash" : tl.playing ? (sel ? "both" : "play") : sel ? "select" : "idle";
    const key = `${state}|${sel ? `${sel.trace}:${sel.index}` : ""}|${this.flashing?.text || ""}`;
    if (key === this.key) {
      // Same content: only the clock moves.
      if (this._time) this._time.textContent = this.ui.formatTime(tl.time);
      if (this._fill) this._fill.style.width = `${(f * 100).toFixed(1)}%`;
      return;
    }
    this.key = key;
    this._time = this._fill = null;
    this.btn.dataset.island = state;
    if (state === "idle") {
      this.state = state;
      this.live.replaceChildren();
      this.btn.style.width = "";
      return;
    }
    const parts = [];
    if (state === "flash") {
      parts.push(h("span", { class: "ov-island-dot" }), h("span", { class: "ov-island-text" }, this.flashing.text));
    } else {
      if (tl.playing) {
        this._time = h("span", { class: "ov-island-time" }, this.ui.formatTime(tl.time));
        this._fill = h("i", { style: { width: `${(f * 100).toFixed(1)}%` } });
        parts.push(h("span", { class: "ov-island-bars", "aria-hidden": "true" }, h("i"), h("i"), h("i"), h("i")),
          this._time, h("span", { class: "ov-island-track" }, this._fill));
      }
      if (sel) {
        const d = sel.kind === "member" ? null : v.describeObject(sel.trace, sel.index);
        const name = sel.kind === "member" ? `${sel.member.cluster} member` : d?.name || "Selection";
        const color = sel.kind === "member" ? sel.member.color : d?.color;
        const dist = sel.kind === "member" ? sel.member.dist : d?.distanceFromSun;
        parts.push(h("span", { class: "ov-island-swatch", style: { background: color || "#fff" } }),
          h("span", { class: "ov-island-text" }, name),
          dist != null ? h("span", { class: "ov-island-sub" }, formatDistance(dist)) : null);
      }
    }
    this.live.replaceChildren(...parts.filter(Boolean));
    // Measure the content and let CSS spring the capsule to that width.
    const want = Math.min(Math.max(this.live.scrollWidth + 44, 170), Math.min(560, this.ui.root.clientWidth - 32));
    this.btn.style.width = `${Math.round(want)}px`;
    this.state = state;
  }
}
