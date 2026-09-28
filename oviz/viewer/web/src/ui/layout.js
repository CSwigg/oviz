// Layouts: an interface style can rearrange the shared components into its
// own structure. Components are moved, never rebuilt, so every feature,
// binding and keyboard shortcut keeps working; each move is recorded and
// undone exactly when the style (or the phone/desktop breakpoint) changes.

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

  apply(skin, narrow) {
    this.reset();
    const build = LAYOUTS[skin.layout];
    if (build) build(this, this.ui, narrow);
    this.current = { skin: skin.id, layout: skin.layout || "default", narrow };
  }
}

const sep = () => h("span", { class: "ov-bar-sep", "aria-hidden": "true" });

const LAYOUTS = {
  /**
   * Focus: the figure and one quiet bar. Time, layers, views, 3D/Sky and a
   * "more" menu share the bottom bar; a plain legend sits in the corner.
   */
  focus(lm, ui) {
    const bar = ui.dock.el.hidden
      ? lm.add(h("div", { class: "ov-dock ov-glass ov-bar-only" }), ui.bottom, ui.bottomRight)
      : ui.dock.el;
    if (!ui.dock.el.hidden) lm.add(sep(), bar);
    lm.move(ui.layersBtn, bar);
    if (!ui.statesBtn.hidden) lm.move(ui.statesBtn, bar);
    if (ui.hasSky) lm.move(ui.viewToggle, bar);
    lm.move(ui.moreBtn, bar);
    lm.move(ui.legend.el, ui.bottomLeft, ui.scale);
  },

  /**
   * Maps: a search box at the top left (layers open under it, and so do an
   * object's details), map buttons stacked on the right, time at the bottom.
   */
  maps(lm, ui, narrow) {
    const box = lm.add(h("div", { class: "ov-searchbox ov-glass" }), ui.top, ui.top.firstChild);
    lm.move(ui.layersBtn, box);
    lm.move(ui.search, box);
    // Phones drop the toolbar: its tools live behind "more" in the box.
    if (narrow) lm.move(ui.moreBtn, box);
    const stack = lm.add(h("div", { class: "ov-mapctl", role: "toolbar", "aria-label": "Map controls" }), ui.ui);
    if (ui.hasSky) {
      const g = h("div", { class: "ov-mapctl-group ov-glass" });
      stack.append(g);
      lm.move(ui.viewToggle, g);
    }
    const zoom = h("div", { class: "ov-mapctl-group ov-glass" });
    const nav = h("div", { class: "ov-mapctl-group ov-glass" });
    stack.append(zoom, nav);
    lm.move(ui.zoomInBtn, zoom);
    lm.move(ui.zoomOutBtn, zoom);
    lm.move(ui.homeBtn, nav);
    lm.move(ui.fsBtn, nav);
    lm.move(ui.legend.el, ui.bottomLeft, ui.scale);
  },

  /**
   * Studio: a sidebar holds the title, search, the selected object and the
   * layers; a timeline runs along the bottom edge. The stage is inset so
   * nothing covers the figure. Phones keep the standard sheets.
   */
  studio(lm, ui, narrow) {
    if (narrow) return;
    const side = lm.add(h("aside", { class: "ov-sidebar ov-chrome", "aria-label": "Figure sidebar" }), ui.ui, ui.ui.firstChild);
    const head = h("div", { class: "ov-side-head" });
    const scroll = h("div", { class: "ov-side-scroll" });
    const foot = h("div", { class: "ov-side-foot" });
    side.append(head, scroll, foot);
    lm.move(ui.brand, head);
    lm.move(ui.sidebarBtn, head);
    lm.move(ui.search, side, scroll);
    lm.move(ui.inspector.el, scroll);
    lm.move(ui.layers.el, scroll);
    lm.move(ui.toolbar, foot);
    const bar = ui.dock.el.hidden
      ? lm.add(h("div", { class: "ov-dock ov-bar-only" }), ui.bottom, ui.bottomRight)
      : ui.dock.el;
    lm.add(sep(), bar);
    if (ui.hasSky) lm.move(ui.viewSeg, bar);
    lm.move(ui.homeBtn, bar);
    lm.move(ui.fsBtn, bar);
  },
};
