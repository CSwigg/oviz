// Layouts: each mode arranges the shared components its own way.
// Components are moved, never rebuilt, so every feature, binding and
// keyboard shortcut keeps working; each move is recorded and undone exactly
// when the mode (or the phone/desktop breakpoint) changes.

import { h } from "./dom.js";

export class LayoutManager {
  constructor(ui) {
    this.ui = ui;
    this.moves = [];
    this.created = [];
    this.current = null;
  }

  /** Move `el` into `parent` (before `before`), remembering where it was. */
  move(el, parent, before = null) {
    if (!el || !parent) return;
    this.moves.push({ el, parent: el.parentNode, next: el.nextSibling });
    parent.insertBefore(el, before);
  }

  /** Insert a layout-only element; it is removed on reset. */
  add(el, parent, before = null) {
    parent.insertBefore(el, before);
    this.created.push(el);
    return el;
  }

  reset() {
    for (const m of this.moves.splice(0).reverse()) {
      if (!m.parent) continue;
      const next = m.next && m.next.parentNode === m.parent ? m.next : null;
      m.parent.insertBefore(m.el, next);
    }
    for (const el of this.created.splice(0)) el.remove();
    this.current = null;
  }

  apply(mode, narrow) {
    this.reset();
    const build = LAYOUTS[mode.layout];
    if (build) build(this, this.ui, narrow);
    this.current = { mode: mode.id, layout: mode.layout || "default", narrow };
  }
}

const sep = () => h("span", { class: "ov-bar-sep", "aria-hidden": "true" });

/** The time bar (a bare bar when the figure has no timeline). */
function timeBar(lm, ui) {
  return ui.dock.el.hidden
    ? lm.add(h("div", { class: "ov-dock ov-glass ov-bar-only" }), ui.bottom, ui.bottomRight)
    : ui.dock.el;
}

const LAYOUTS = {
  /**
   * Focus: the figure and one quiet bar. Time, layers, views, 3D/Sky, the
   * mode switch and a "more" menu share the bottom bar; a plain key sits in
   * the corner. The top keeps only the title and a search button.
   */
  focus(lm, ui) {
    const bar = timeBar(lm, ui);
    if (!ui.dock.el.hidden) lm.add(sep(), bar);
    lm.move(ui.layersBtn, bar);
    if (!ui.statesBtn.hidden) lm.move(ui.statesBtn, bar);
    if (ui.hasSky) lm.move(ui.viewToggle, bar);
    lm.move(ui.modeBtn, bar);
    lm.move(ui.moreBtn, bar);
    lm.move(ui.legend.el, ui.bottomLeft, ui.scaleRow);
  },

  /**
   * Detailed: every control in view. The toolbar stays at the top, the
   * layers panel is open, the full transport and 3D/Sky switch are on the
   * bottom; the key shows whenever the layers panel is closed.
   */
  detailed(lm, ui) {
    const bar = timeBar(lm, ui);
    if (!ui.dock.el.hidden) lm.add(sep(), bar);
    lm.move(ui.modeBtn, bar);
    lm.move(ui.legend.el, ui.bottomLeft, ui.scaleRow);
  },
};
