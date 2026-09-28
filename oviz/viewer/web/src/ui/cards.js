// Cards skin: panels behave like physical cards. Picking a new object deals
// the previous inspector card back into a stack behind the new one (tap a
// card in the stack to return to that object), and cards can be flicked
// away: the inspector to the right, the layers card to the left.

import { flickToDismiss, springEasing, SPRINGS, reducedMotion } from "./motion.js";

const MAX_STACK = 3;

export class CardStack {
  constructor(ui) {
    this.ui = ui;
    this.ghosts = [];
    const insp = ui.inspector.el;
    flickToDismiss(insp, insp.querySelector(".ov-panel-head"), {
      axis: "x", sign: 1,
      enabled: () => this.active && !ui.narrow,
      onDismiss: () => ui.select(null),
    });
    const layers = ui.layers.el;
    const head = layers.querySelector(".ov-panel-head");
    if (head) {
      flickToDismiss(layers, head, {
        axis: "x", sign: -1,
        enabled: () => this.active && !ui.narrow,
        onDismiss: () => ui.setLayersOpen(false),
      });
    }
    ui.viewer.on("skin", () => { if (!this.active) this.clear(false); });
  }

  get active() {
    return document.documentElement.dataset.ovizSkin === "cards";
  }

  /** The inspector is about to switch from `prev` to `next`. */
  push(prev, next) {
    if (!this.active || this.ui.narrow) return;
    if (!next) { this.clear(true); return; }
    const same = (a, b) => a && b && a.trace === b.trace && a.index === b.index && a.kind === b.kind;
    // Returning to a card in the stack lifts it out of the stack.
    const back = this.ghosts.find((g) => same(g._hit, next));
    if (back) this.drop(back);
    if (!prev || same(prev, next) || this.ui.inspector.el.dataset.open !== "true") return;
    const src = this.ui.inspector.el;
    const g = src.cloneNode(true);
    g.classList.add("ov-card-ghost");
    g.removeAttribute("data-tilting");
    g.removeAttribute("style");
    g.setAttribute("aria-hidden", "true");
    for (const el of g.querySelectorAll("[id]")) el.removeAttribute("id");
    // Canvases clone blank; copy the chart pixels across.
    const from = src.querySelectorAll("canvas");
    g.querySelectorAll("canvas").forEach((c, i) => {
      try { c.getContext("2d").drawImage(from[i], 0, 0); } catch (_) { /* size mismatch */ }
    });
    g._hit = prev;
    g.addEventListener("click", () => this.ui.select(g._hit));
    src.before(g);
    this.ghosts.unshift(g);
    while (this.ghosts.length > MAX_STACK) this.drop(this.ghosts[this.ghosts.length - 1]);
    this.layout();
    // Deal the new card in on top.
    if (!reducedMotion()) {
      const s = springEasing(SPRINGS.wobble);
      src.animate([{ translate: "46px 10px", rotate: "3deg", opacity: 0.4 }, { translate: "0 0", rotate: "0deg", opacity: 1 }],
        { duration: s.duration, easing: s.easing });
      g.animate([{ translate: "0 12px" }, { translate: "0 0" }], { duration: s.duration, easing: s.easing });
    }
  }

  layout() {
    this.ghosts.forEach((g, i) => {
      g.style.setProperty("--d", String(i + 1));
      g.dataset.depth = String(i + 1);
    });
  }

  drop(g) {
    this.ghosts = this.ghosts.filter((x) => x !== g);
    const done = () => g.remove();
    if (reducedMotion() || !g.animate) done();
    else g.animate([{ opacity: 1 }, { opacity: 0, translate: "0 -24px" }], { duration: 220, easing: "ease-in" }).finished.then(done, done);
    this.layout();
  }

  clear(animate) {
    for (const g of [...this.ghosts]) {
      if (animate) this.drop(g);
      else g.remove();
    }
    this.ghosts = [];
  }
}
