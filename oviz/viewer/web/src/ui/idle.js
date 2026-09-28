// Idle fade: in styles that ask for it, the controls fade out while the
// pointer rests over the figure (like a video player) and return the moment
// it moves. Never while a menu, panel or story is open, while the pointer is
// over the controls, on touch screens, or while presenting.

const CHROME = ".ov-top, .ov-bottom, .ov-panel, .ov-pop, .ov-callout, .ov-story, .ov-legend, .ov-mapctl, .ov-sidebar";

export class IdleFade {
  constructor(ui, { delay = 2800 } = {}) {
    this.ui = ui;
    this.delay = delay;
    this.timer = 0;
    this.over = false;
    this.fine = window.matchMedia?.("(hover: hover) and (pointer: fine)");
    const root = ui.root;
    const wake = () => this.wake();
    for (const type of ["pointermove", "pointerdown", "wheel", "keydown"]) root.addEventListener(type, wake, { passive: true });
    root.addEventListener("pointerover", (e) => { this.over = !!e.target.closest?.(CHROME); }, { passive: true });
    root.addEventListener("pointerleave", () => { this.over = false; });
  }

  get enabled() {
    return !!this.ui.skinConfig?.autoHide && !!this.fine?.matches;
  }

  busy() {
    const ui = this.ui;
    const ds = ui.root.dataset;
    return this.over || ui.layersOpen || ui.overlays?.size || ui.palette?.open || ui._menu || ui._sheet
      || ds.story === "true" || ds.presenting === "true" || ds.lasso === "true" || ui.viewer.controls.pointers.size > 0;
  }

  wake() {
    if (this.ui.root.dataset.idle) delete this.ui.root.dataset.idle;
    clearTimeout(this.timer);
    if (this.enabled) this.timer = setTimeout(() => this.sleep(), this.delay);
  }

  sleep() {
    if (!this.enabled) return;
    if (this.busy()) {
      this.timer = setTimeout(() => this.sleep(), this.delay);
      return;
    }
    this.ui.root.dataset.idle = "true";
  }
}
