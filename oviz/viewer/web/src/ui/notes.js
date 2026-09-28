// Notes: text annotations pinned to an object (they follow its orbit through
// time) or to a fixed point. Captured in saved views; editable in place.

import { h, icon } from "./dom.js";
import { cloneJson } from "../app/state.js";
import { uid } from "../app/states.js";

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
  }

  commands() {
    return [
      { title: "Add note at view centre", icon: "edit", keywords: "annotate label text", run: () => this.addAt(this.viewer.pose.target.slice()) },
      ...(this.notes.length ? [{ title: "Remove all notes", icon: "trash", run: () => this.setNotes([]) }] : []),
    ];
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
      const get = n.anchor.pos
        ? () => n.anchor.pos
        : () => v.objectPosition(n.anchor.trace, n.anchor.index);
      labels.setDynamic(id, n.text, get, "ov-label--note");
      const el = labels.dynamic.get(id).el;
      el.dataset.noteId = n.id;
      el.style.color = n.color || "#ffffff";
      el._w = 0;
    }
    v.renderer.invalidate();
  }

  edit(id, anchorEl) {
    const n = this.notes.find((x) => x.id === id);
    if (!n) return;
    const el = anchorEl || this.viewer.labels.dynamic.get(`note:${id}`)?.el;
    if (!el) return;
    const input = h("input", { class: "ov-input", value: n.text, "aria-label": "Note text" });
    input.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Enter") { commit(); this.ui.closeMenu(); }
      if (e.key === "Escape") this.ui.closeMenu();
    });
    const commit = () => {
      const t = input.value.trim();
      if (t && t !== n.text) { n.text = t; this.sync(); }
    };
    const colors = ["#ffffff", "#f4c46a", "#7cc4ff", "#ff7a70", "#6fdc9b"].map((c) =>
      h("button", { class: "ov-chip", type: "button", style: { background: c }, "aria-label": c, "aria-pressed": String(n.color === c), onclick: () => { n.color = c; this.sync(); } }));
    const body = h("div", { class: "ov-pop-body" }, input, h("div", { class: "ov-colors" }, colors));
    this.ui.menu(el, [
      { body },
      { label: "Delete note", icon: "trash", run: () => { this.notes = this.notes.filter((x) => x !== n); this.viewer.labels.removeDynamic(`note:${n.id}`); this.sync(); } },
    ], { title: "Note", align: "left" });
    input.addEventListener("blur", commit);
    requestAnimationFrame(() => { input.focus(); input.select(); });
  }
}

function fromLegacy(l) {
  if (!l || typeof l !== "object") return null;
  const x = Number(l.x ?? l.position?.x), y = Number(l.y ?? l.position?.y), z = Number(l.z ?? l.position?.z);
  if (![x, y, z].every(Number.isFinite)) return null;
  return { id: String(l.id || uid()), text: String(l.text || "Note"), anchor: { pos: [x, y, z] }, color: String(l.color || "#ffffff") };
}
