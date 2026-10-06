// Selection reticle: a ring around the selected object, recomputed after
// every rendered frame so it tracks the object through camera moves and
// time playback. In Focus mode a small callout card sits beside it instead
// of the inspector opening (details stay one click away).

import { h, iconButton, fmt, fmtInt } from "./dom.js";
import { formatDistance, clamp } from "../core/math.js";

const SVGNS = "http://www.w3.org/2000/svg";
const R = 13;

function el(tag, attrs = {}) {
  const e = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  return e;
}

export class SelectionReticle {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.svg = el("svg", { class: "ov-reticle", "aria-hidden": "true" });
    this.mark = el("g", { class: "ov-reticle-mark", opacity: 0 });
    this.svg.append(this.mark);
    ui.ui.prepend(this.svg);
    this.callout = h("div", { class: "ov-callout ov-glass", role: "dialog", "aria-label": "Selected object", hidden: true });
    ui.ui.append(this.callout);
    this.scr = [0, 0, 0];
    this.visible = false;
    v.renderer.afterRender.push(() => this.update());
    v.on("select", () => this.onSelect());
    v.on("viewmode", () => this.update());
    window.addEventListener("resize", () => this.update());
  }

  get calloutMode() {
    return this.ui.modeConfig?.selection === "callout";
  }

  onSelect() {
    const sel = this.ui.selection;
    this.mark.replaceChildren();
    this.callout.hidden = true;
    if (!sel) { this.hide(); return; }
    if (this.calloutMode) {
      const ready = sel.kind === "member" ? Promise.resolve() : this.viewer.objectMeta(sel.trace);
      ready.then(() => {
        if (this.ui.selection !== sel) return;
        this.buildCallout(sel);
        this.update();
      });
    }
    // Rebuilt per selection so the lock-on animation replays.
    const lock = el("g", { class: "ov-reticle-lock" });
    lock.append(el("circle", { class: "ov-reticle-halo", r: R }));
    this.mark.append(lock);
    this.update();
  }

  buildCallout(sel) {
    const v = this.viewer;
    let kicker, name, color;
    const facts = [];
    if (sel.kind === "member") {
      kicker = `Member of ${sel.member.cluster}`;
      name = `Star ${sel.index + 1}`;
      color = sel.member.color;
      if (sel.member.dist != null) facts.push(formatDistance(sel.member.dist));
    } else {
      const d = v.describeObject(sel.trace, sel.index);
      if (!d) return;
      kicker = d.traceName;
      name = d.name;
      color = d.color;
      if (d.ageNow != null) facts.push(`${fmt(d.ageNow, d.ageNow >= 100 ? 0 : 1)} Myr old`);
      if (d.dist != null) facts.push(formatDistance(d.dist));
      if (d.nStars != null) facts.push(`${fmtInt(d.nStars)} stars`);
    }
    // Native replaceChildren would print a null child as "null": filter them.
    this.callout.replaceChildren(...[
      h("div", { class: "ov-callout-head" },
        h("span", { class: "ov-callout-dot", style: { background: color } }),
        h("span", { class: "ov-callout-kicker" }, kicker),
        iconButton("close", "Clear selection", () => this.ui.select(null), { cls: "ov-callout-close", shortcut: "Esc" })),
      h("div", { class: "ov-callout-name" }, name),
      facts.length ? h("div", { class: "ov-callout-facts" }, facts.join(" · ")) : null,
      h("div", { class: "ov-callout-actions" },
        // Member stars are seen from the Sun: there is nowhere to fly to.
        sel.kind === "member" ? null : h("button", { class: "ov-callout-btn", type: "button", onclick: () => v.flyToObject(sel.trace, sel.index) }, "Fly to"),
        h("button", { class: "ov-callout-btn ov-callout-btn--primary", type: "button", onclick: () => this.ui.openDetails() }, "Details")),
    ].filter(Boolean));
    this.callout.hidden = false;
  }

  placeCallout(x, y) {
    const c = this.callout;
    if (c.hidden) return;
    const W = this.ui.root.clientWidth, H = this.ui.root.clientHeight;
    const cw = c.offsetWidth, ch = c.offsetHeight;
    let cx = x + 26;
    if (cx + cw > W - 12) cx = x - 26 - cw;
    cx = clamp(cx, 12, Math.max(12, W - cw - 12));
    const cy = clamp(y - ch / 2, 64, Math.max(64, H - ch - 100));
    c.style.transform = `translate(${cx.toFixed(0)}px, ${cy.toFixed(0)}px)`;
    c.dataset.side = cx < x ? "left" : "right";
  }

  hide() {
    // Always hide the card: a selection that starts off-screen must not
    // show a freshly built card at its last (or initial 0, 0) position.
    this.callout.style.visibility = "hidden";
    delete this.ui.root.dataset.callout;
    if (!this.visible) return;
    this.visible = false;
    this.mark.setAttribute("opacity", "0");
    delete this.ui.root.dataset.reticle;
  }

  update() {
    const sel = this.ui.selection;
    const v = this.viewer;
    if (!sel) return this.hide();
    const p = this.ui.selectionPosition(sel);
    const cam = v.renderer.camera;
    const s = p && cam.project(p, this.scr);
    const W = cam.width, H = cam.height;
    if (!s || s[0] < -40 || s[1] < -40 || s[0] > W + 40 || s[1] > H + 40) return this.hide();
    // Work in page space (the canvas need not start at the page origin).
    const cr = v.canvas.getBoundingClientRect();
    const rr = this.ui.root.getBoundingClientRect();
    const x = s[0] + cr.left - rr.left;
    const y = s[1] + cr.top - rr.top;
    this.mark.setAttribute("transform", `translate(${x.toFixed(1)} ${y.toFixed(1)})`);
    this.mark.setAttribute("opacity", "1");
    this.visible = true;
    this.ui.root.dataset.reticle = "true";
    // The label card steps aside while the details panel shows the object.
    const details = this.ui.inspector?.el.dataset.open === "true";
    const showCallout = !this.callout.hidden && !details;
    this.callout.style.visibility = showCallout ? "" : "hidden";
    if (showCallout) this.ui.root.dataset.callout = "true";
    else delete this.ui.root.dataset.callout;
    this.placeCallout(x, y);
  }
}
