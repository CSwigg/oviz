// Notes: text annotations pinned to an object (they follow its orbit through
// time) or to a fixed point, like the classic viewer's text labels. Each has
// a text size and can be hidden; fixed notes have X/Y/Z in pc. Captured in
// saved views, editable in place, listed in the Notes widget.

import { h, icon, iconButton, clear } from "./dom.js";
import { slider, numberInput } from "./controls.js";
import { cloneJson } from "../app/state.js";
import { uid } from "../app/states.js";

const NOTE_SIZE = 12.5; // px: the note chip's own size

export class NotesPlugin {
  constructor() {
    this.name = "notes";
    this.notes = [];
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    v.stateExtensions.set("notes", { capture: () => cloneJson(this.notes), apply: (n) => this.setNotes(n || []) });
    ui.notes = this;
    // Classic manual labels in the initial state become notes.
    const legacy = v.manifest.initialState?.manual_labels;
    if (Array.isArray(legacy)) this.setNotes(legacy.map(fromLegacy).filter(Boolean));
    v.labels.el.addEventListener("click", (e) => {
      const el = e.target.closest?.(".ov-label--note");
      if (el) this.edit(el.dataset.noteId, el);
    });
    // Place a note where you click (classic "Place in scene").
    v.canvas.addEventListener("click", (e) => {
      if (!this.placing) return;
      this.placing = false;
      ui.root.dataset.placing = "false";
      const r = v.canvas.getBoundingClientRect();
      const x = e.clientX - r.left, y = e.clientY - r.top;
      const hit = v.pick(x, y);
      if (hit && hit.kind !== "member") { this.addForObject(hit); return; }
      const pos = v.pickSurface(x, y) || this.onFocusPlane(x, y);
      if (pos) this.addAt(pos);
    }, true);
    ui.widgets?.register({
      id: "notes", title: "Notes", icon: "edit", description: "Text labels in the scene",
      size: { w: 330, h: 300 }, minSize: { w: 260, h: 180 },
      mount: (panel) => this.mountWidget(panel),
    });
  }

  /** The point under the pointer on the plane through the orbit centre, facing the camera. */
  onFocusPlane(x, y) {
    const cam = this.viewer.renderer.camera;
    const { origin: o, dir: d } = cam.ray(x, y);
    const t0 = this.viewer.pose.target, f = cam.forward;
    const den = d[0] * f[0] + d[1] * f[1] + d[2] * f[2];
    if (Math.abs(den) < 1e-9) return null;
    const t = ((t0[0] - o[0]) * f[0] + (t0[1] - o[1]) * f[1] + (t0[2] - o[2]) * f[2]) / den;
    return t > 0 ? [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t] : null;
  }

  place() {
    this.placing = true;
    this.ui.root.dataset.placing = "true";
    this.ui.toast("Click in the scene to place a note (on an object, it follows the object)", { icon: icon("edit"), ms: 2600 });
  }

  commands() {
    const hidden = this.notes.filter((n) => n.visible === false).length;
    return [
      { title: "Place a note (click in the scene)", icon: "edit", keywords: "annotate label text add", run: () => this.place() },
      { title: "Add note at view centre", icon: "edit", keywords: "annotate label text", run: () => this.addAt(this.viewer.pose.target.slice()) },
      ...(this.notes.length ? [{ title: "Notes…", sub: `${this.notes.length} note${this.notes.length === 1 ? "" : "s"}${hidden ? `, ${hidden} hidden` : ""}`, icon: "edit", keywords: "labels text list", run: () => this.ui.widgets.open("notes") }] : []),
      ...(hidden ? [{ title: "Show hidden notes", icon: "eye", run: () => this.patchAll({ visible: true }) }] : []),
      ...(this.notes.length ? [{ title: "Remove all notes", icon: "trash", run: () => this.setNotes([]) }] : []),
    ];
  }

  patch(id, fields) {
    const n = this.notes.find((x) => x.id === id);
    if (!n) return;
    Object.assign(n, fields);
    this.sync();
  }

  patchAll(fields) {
    for (const n of this.notes) Object.assign(n, fields);
    this.sync();
  }

  remove(id) {
    this.viewer.labels.removeDynamic(`note:${id}`);
    this.notes = this.notes.filter((x) => x.id !== id);
    this.sync();
  }

  addForObject(hit) {
    const d = this.viewer.describeObject(hit.trace, hit.index);
    this.add({ id: uid(), text: d?.name || "Note", anchor: { trace: hit.trace, index: hit.index }, color: "#ffffff" });
  }

  addAt(pos) {
    this.add({ id: uid(), text: "Note", anchor: { pos }, color: "#ffffff" }, true);
  }

  add(note, editNow = false) {
    this.notes.push(note);
    this.sync();
    this.ui.toast("Note added — click it to edit", { icon: icon("edit"), ms: 1800 });
    if (editNow) requestAnimationFrame(() => this.edit(note.id));
  }

  setNotes(list) {
    for (const n of this.notes) this.viewer.labels.removeDynamic(`note:${n.id}`);
    this.notes = cloneJson(list);
    this.sync();
  }

  sync() {
    const v = this.viewer;
    const labels = v.labels;
    for (const n of this.notes) {
      const id = `note:${n.id}`;
      // Hidden notes stay in the list (and in States) but leave the scene.
      if (n.visible === false) { labels.removeDynamic(id); continue; }
      const get = n.anchor.pos
        ? () => n.anchor.pos
        : () => v.objectPosition(n.anchor.trace, n.anchor.index);
      labels.setDynamic(id, n.text, get, "ov-label--note");
      const el = labels.dynamic.get(id).el;
      el.dataset.noteId = n.id;
      el.style.color = n.color || "#ffffff";
      el.style.fontSize = Number.isFinite(n.size) && n.size > 0 ? `${n.size}px` : "";
      el._w = 0;
    }
    v.renderer.invalidate();
    this.renderWidget?.();
  }

  edit(id, anchorEl) {
    const n = this.notes.find((x) => x.id === id);
    if (!n) return;
    const el = anchorEl || this.viewer.labels.dynamic.get(`note:${id}`)?.el;
    if (!el) return;
    const input = h("input", { class: "ov-input", value: n.text, "aria-label": "Note text" });
    let cancelled = false;
    input.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Enter") { commit(); this.ui.closeMenu(); }
      // Closing the menu blurs the field; the blur must not commit.
      if (e.key === "Escape") { cancelled = true; this.ui.closeMenu(); }
    });
    const commit = () => {
      if (cancelled) return;
      const t = input.value.trim();
      if (t && t !== n.text) { n.text = t; this.sync(); }
    };
    const colors = ["#ffffff", "#f4c46a", "#7cc4ff", "#ff7a70", "#6fdc9b"].map((c) =>
      h("button", {
        class: "ov-chip", type: "button", style: { background: c }, "aria-label": c, "aria-pressed": String(n.color === c),
        onclick: (e) => {
          n.color = c;
          this.sync();
          for (const b of colors) b.setAttribute("aria-pressed", String(b === e.currentTarget));
        },
      }));
    const size = slider({ label: "Text size", min: 6, max: 96, value: n.size || NOTE_SIZE, scale: "log", format: (x) => `${x.toFixed(x < 20 ? 1 : 0)} px`, onInput: (x) => { n.size = +x.toFixed(1); this.sync(); } });
    const fields = [input, h("div", { class: "ov-colors" }, colors), size];
    if (n.anchor.pos) {
      // A fixed note sits at X/Y/Z in pc (classic text label fields).
      const axis = (k, label) => numberInput({ label, value: +Number(n.anchor.pos[k]).toFixed(2), step: "any", onChange: (x) => { n.anchor.pos[k] = x; this.sync(); } });
      fields.push(h("div", { class: "ov-field-row ov-note-xyz" }, axis(0, "X pc"), axis(1, "Y pc"), axis(2, "Z pc")));
    }
    const body = h("div", { class: "ov-pop-body" }, ...fields);
    this.ui.menu(el, [
      { body },
      { label: "Hide note", icon: "eyeOff", run: () => this.patch(n.id, { visible: false }) },
      { label: "Delete note", icon: "trash", run: () => this.remove(n.id) },
    ], { title: "Note", align: "left" });
    input.addEventListener("blur", commit);
    requestAnimationFrame(() => { input.focus(); input.select(); });
  }

  // ---------------------------------------------------------------- widget

  mountWidget(panel) {
    this.widgetBody = panel.body;
    panel.addHeadButton(iconButton("plus", "Place a note (click in the scene)", () => this.place()));
    this.renderWidget = () => {
      const body = this.widgetBody;
      if (!body || !body.isConnected) return;
      clear(body);
      if (!this.notes.length) {
        body.append(h("div", { class: "ov-note", style: { padding: "12px" } }, "No notes yet. Use + to place one, or Add note in an object's details."));
        return;
      }
      for (const n of this.notes) {
        const vis = n.visible !== false;
        const name = h("button", { class: "ov-row-name", type: "button", title: "Edit" }, n.text);
        name.addEventListener("click", () => { if (!vis) this.patch(n.id, { visible: true }); requestAnimationFrame(() => this.edit(n.id)); });
        body.append(h("div", { class: "ov-row", "data-visible": String(vis) },
          h("span", { class: "ov-swatch-btn" }, h("span", { class: "ov-swatch", style: { color: n.color || "#fff" } })),
          name,
          h("span", { class: "ov-row-count" }, n.anchor.pos ? "Fixed" : "On object"),
          h("div", { class: "ov-row-actions" },
            iconButton(vis ? "eye" : "eyeOff", vis ? "Hide note" : "Show note", () => this.patch(n.id, { visible: !vis }), { cls: "ov-eye" }),
            iconButton("target", "Go to note", () => this.goTo(n)),
            iconButton("trash", "Delete note", () => this.remove(n.id)))));
      }
    };
    return { onOpen: () => this.renderWidget(), destroy: () => { this.widgetBody = null; } };
  }

  /** Centre the orbit on a note, keeping the zoom. */
  goTo(n) {
    const v = this.viewer;
    const pos = n.anchor.pos || v.objectPosition(n.anchor.trace, n.anchor.index);
    if (pos) v.recenter?.(pos);
  }
}

function fromLegacy(l) {
  if (!l || typeof l !== "object") return null;
  const x = Number(l.x ?? l.position?.x), y = Number(l.y ?? l.position?.y), z = Number(l.z ?? l.position?.z);
  if (![x, y, z].every(Number.isFinite)) return null;
  const size = Number(l.size ?? l.font_size);
  return {
    id: String(l.id || uid()), text: String(l.text || "Note"), anchor: { pos: [x, y, z] }, color: String(l.color || "#ffffff"),
    // Classic labels keep their size and hidden state.
    ...(Number.isFinite(size) && size > 0 ? { size } : {}),
    ...(l.visible === false ? { visible: false } : {}),
  };
}
