// Widgets: floating tool panels (the classic viewer's Widgets menu), such as
// the Relative SFH, the Birth tree and the Sky view.
//
// A plugin registers a widget with `ui.widgets.register(def)`:
//   { id, title, icon, description, size: {w, h}, minSize: {w, h},
//     available: () => bool,
//     mount: (panel) => ({ onOpen, onClose, onResize(w, h), capture, apply, destroy }) }
// `panel` is { id, el, head, body, host, setTitle(text), addHeadButton(node), close() }.
// Panels open lazily, can be dragged by their title and resized from the
// corner, come to the front when touched, and on phones become bottom sheets
// (one at a time). Which widgets are open, where, and each widget's own
// `capture()` go into States.

import { h, iconButton, clear } from "./dom.js";
import { clamp } from "../core/math.js";

const GAP = 12;

export class WidgetHost {
  constructor(ui) {
    this.ui = ui;
    this.defs = new Map();
    this.panels = new Map(); // id → {def, el, ctl, rect, open}
    this.listeners = new Set();
    this.z = 5;
    this.layer = h("div", { class: "ov-widgets" });
    ui.ui.append(this.layer);
    window.addEventListener("resize", () => this.fitAll());
  }

  /** Register a widget (see the module comment). Re-registering replaces it. */
  register(def) {
    if (!def?.id || typeof def.mount !== "function") throw new Error("widget needs an id and mount()");
    this.defs.set(def.id, def);
    this.changed();
    return def;
  }

  onChange(fn) {
    this.listeners.add(fn);
    return () => this.listeners.delete(fn);
  }

  changed() {
    for (const fn of this.listeners) fn();
  }

  /** Widgets this figure can show. */
  available() {
    return [...this.defs.values()].filter((d) => !d.available || d.available());
  }

  isOpen(id) {
    return !!this.panels.get(id)?.open;
  }

  toggle(id) {
    return this.isOpen(id) ? this.close(id) : this.open(id);
  }

  open(id, { rect, state } = {}) {
    const def = this.defs.get(id);
    if (!def || (def.available && !def.available())) return null;
    let p = this.panels.get(id);
    if (!p) p = this.build(def);
    // Phones show one sheet at a time.
    if (this.ui.narrow) for (const other of this.panels.values()) if (other !== p && other.open) this.close(other.def.id);
    if (rect) p.rect = { ...p.rect, ...rect };
    this.place(p);
    p.open = true;
    p.el.dataset.open = "true";
    this.front(p);
    p.ctl.onOpen?.();
    if (state !== undefined) p.ctl.apply?.(state);
    this.resized(p);
    this.ui.idle?.wake?.();
    this.changed();
    return p;
  }

  close(id) {
    const p = this.panels.get(id);
    if (!p || !p.open) return;
    p.open = false;
    p.el.dataset.open = "false";
    p.ctl.onClose?.();
    this.changed();
  }

  closeAll() {
    for (const id of [...this.panels.keys()]) this.close(id);
  }

  build(def) {
    const title = h("span", { class: "ov-panel-title" }, def.title);
    const actions = h("div", { class: "ov-widget-actions" });
    const head = h("div", { class: "ov-panel-head ov-widget-head" }, title, actions,
      iconButton("close", `Close ${def.title}`, () => this.close(def.id)));
    const body = h("div", { class: "ov-widget-body" });
    const grip = h("span", { class: "ov-widget-grip", "aria-hidden": "true" });
    const el = h("section", { class: "ov-panel ov-widget ov-glass ov-chrome", "data-open": "false", "data-widget": def.id, "aria-label": def.title }, head, body, grip);
    this.layer.append(el);
    const size = def.size || { w: 380, h: 300 };
    const n = this.panels.size;
    const rect = { x: null, y: null, w: size.w, h: size.h, cascade: n };
    const p = { def, el, head, body, rect, open: false, ctl: null };
    const panel = {
      id: def.id, el, head, body, host: this,
      setTitle: (text) => { title.textContent = text; },
      addHeadButton: (node) => actions.append(node),
      close: () => this.close(def.id),
    };
    p.ctl = def.mount(panel) || {};
    this.panels.set(def.id, p);
    el.addEventListener("pointerdown", () => this.front(p), true);
    this.bindDrag(p);
    this.bindResize(p, grip);
    new ResizeObserver(() => { if (p.open) this.resized(p); }).observe(body);
    return p;
  }

  front(p) {
    p.el.style.zIndex = String(++this.z);
  }

  /** Stage bounds the panels stay inside (below the top bar, above the dock). */
  bounds() {
    const r = this.ui.root.getBoundingClientRect();
    return { w: r.width, h: r.height, top: 60, bottom: Math.max(120, r.height - 96) };
  }

  place(p) {
    const b = this.bounds();
    const minW = p.def.minSize?.w ?? 240, minH = p.def.minSize?.h ?? 160;
    const r = p.rect;
    r.w = clamp(r.w, minW, Math.max(minW, b.w - 2 * GAP));
    r.h = clamp(r.h, minH, Math.max(minH, b.bottom - b.top));
    if (r.x == null || r.y == null) {
      // New panels cascade down the right-hand side.
      r.x = b.w - r.w - GAP - (r.cascade % 4) * 28;
      r.y = b.top + (r.cascade % 4) * 28;
    }
    r.x = clamp(r.x, GAP, Math.max(GAP, b.w - r.w - GAP));
    r.y = clamp(r.y, b.top - 48, Math.max(b.top - 48, b.bottom - Math.min(r.h, 120)));
    Object.assign(p.el.style, { left: `${r.x}px`, top: `${r.y}px`, width: `${r.w}px`, height: `${r.h}px` });
  }

  fitAll() {
    for (const p of this.panels.values()) if (p.open) this.place(p);
  }

  resized(p) {
    const b = p.body;
    p.ctl.onResize?.(b.clientWidth, b.clientHeight);
  }

  bindDrag(p) {
    let start = null;
    p.head.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 || this.ui.narrow || e.target.closest("button, input, select, a")) return;
      start = { x: e.clientX, y: e.clientY, rx: p.rect.x, ry: p.rect.y };
      p.head.setPointerCapture(e.pointerId);
      p.el.dataset.dragging = "true";
    });
    p.head.addEventListener("pointermove", (e) => {
      if (!start) return;
      p.rect.x = start.rx + e.clientX - start.x;
      p.rect.y = start.ry + e.clientY - start.y;
      this.place(p);
    });
    const end = () => { start = null; delete p.el.dataset.dragging; };
    p.head.addEventListener("pointerup", end);
    p.head.addEventListener("pointercancel", end);
  }

  bindResize(p, grip) {
    let start = null;
    grip.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      e.stopPropagation();
      start = { x: e.clientX, y: e.clientY, w: p.rect.w, h: p.rect.h };
      grip.setPointerCapture(e.pointerId);
      p.el.dataset.dragging = "true";
    });
    grip.addEventListener("pointermove", (e) => {
      if (!start) return;
      p.rect.w = start.w + e.clientX - start.x;
      p.rect.h = start.h + e.clientY - start.y;
      this.place(p);
    });
    const end = () => { start = null; delete p.el.dataset.dragging; };
    grip.addEventListener("pointerup", end);
    grip.addEventListener("pointercancel", end);
  }

  // ------------------------------------------------------------ menus

  /** Items for the Widgets menu: one per available widget, ticked when open. */
  menuItems() {
    return this.available().map((d) => ({
      label: d.title,
      icon: d.icon || "widgets",
      checked: this.isOpen(d.id),
      run: () => this.toggle(d.id),
    }));
  }

  commands() {
    return this.available().map((d) => ({
      title: `${this.isOpen(d.id) ? "Close" : "Open"} ${d.title}`,
      sub: d.description || "",
      icon: d.icon || "widgets",
      keywords: `widget panel ${d.keywords || ""}`,
      run: () => this.toggle(d.id),
    }));
  }

  // ------------------------------------------------------------ States

  capture() {
    const open = {};
    for (const [id, p] of this.panels) {
      if (!p.open) continue;
      const { x, y, w, h: hh } = p.rect;
      open[id] = { rect: { x, y, w, h: hh }, state: p.ctl.capture ? p.ctl.capture() : null };
    }
    return Object.keys(open).length ? { open } : null;
  }

  apply(s) {
    const want = s && typeof s === "object" && s.open && typeof s.open === "object" ? s.open : {};
    for (const [id, p] of this.panels) if (p.open && !want[id]) this.close(id);
    for (const [id, w] of Object.entries(want)) {
      const rect = w?.rect && typeof w.rect === "object" ? w.rect : undefined;
      this.open(id, { rect, state: w?.state ?? null });
    }
  }
}

/** Clear and fill a widget body. */
export function widgetContent(body, ...children) {
  clear(body);
  body.append(...children.filter(Boolean));
  return body;
}
