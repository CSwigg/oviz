// Text labels as positioned DOM nodes: crisp at every zoom, selectable,
// accessible, and cheap for the tens-to-hundreds of labels a figure uses.
// Positions update only when the camera or time changes.

import { frameOffset, cpuFramePosition } from "../engine/frames.js";

export class LabelsOverlay {
  constructor(container) {
    this.el = document.createElement("div");
    this.el.className = "ov-labels";
    this.el.setAttribute("aria-hidden", "true");
    container.appendChild(this.el);
    this.groups = [];
    this.params = null;
    this._key = "";
    this.dynamic = new Map(); // id → {el, getPosition}
  }

  addTrace(trace, data) {
    const L = trace.labels;
    const nodes = [];
    const family = L.family || "";
    for (let i = 0; i < L.count; i++) {
      const n = document.createElement("div");
      n.className = "ov-label";
      n.textContent = L.text[i];
      n.style.color = L.color[i] || "#e2e8f0";
      const px = Math.max(9, Math.min(64, (L.size[i] || 14) * 0.5));
      n.style.fontSize = `${px}px`;
      if (family && !/helvetica|arial|inter/i.test(family)) n.style.fontFamily = family;
      this.el.appendChild(n);
      nodes.push(n);
    }
    this.groups.push({
      key: trace.key, trace, nodes, positions: data.position,
      posFrames: L.position.frames || 1, count: L.count, offsets: data.offset || null,
    });
  }

  /** Add or update a transient label (hover/selection/measurement). */
  setDynamic(id, text, getPosition, className = "") {
    let d = this.dynamic.get(id);
    if (!d) {
      const el = document.createElement("div");
      el.className = `ov-label ov-label--dyn ${className}`;
      this.el.appendChild(el);
      d = { el, getPosition };
      this.dynamic.set(id, d);
    }
    d.getPosition = getPosition;
    d.el.textContent = text;
    this._key = "";
  }

  removeDynamic(id) {
    const d = this.dynamic.get(id);
    if (!d) return;
    d.el.remove();
    this.dynamic.delete(id);
  }

  update(frame) {
    const p = this.params;
    if (!p) return;
    const cam = frame.camera;
    const key = `${cam.version}|${p.frame.toFixed(4)}|${p.labelSignature}|${this.dynamic.size}`;
    if (key === this._key) return;
    this._key = key;
    const tmp = [0, 0, 0];
    const scr = [0, 0, 0];
    const W = frame.cssWidth, H = frame.cssHeight;
    // Greedy de-overlap: transient labels first, then traces in order.
    const placed = [];
    const fits = (x, y, w, h) => {
      for (const b of placed) {
        if (x < b[0] + b[2] && x + w > b[0] && y < b[1] + b[3] && y + h > b[1]) return false;
      }
      placed.push([x, y, w, h]);
      return true;
    };
    for (const d of this.dynamic.values()) {
      const pos = d.getPosition();
      const s = pos && cam.project(pos, scr);
      if (!s) { d.el.style.opacity = "0"; continue; }
      d.el.style.opacity = "1";
      d.el.style.transform = `translate3d(${s[0].toFixed(1)}px, ${s[1].toFixed(1)}px, 0)`;
      const w = d.el._w || (d.el._w = d.el.offsetWidth || 80);
      placed.push([s[0], s[1] - 12, w, 22]);
    }
    for (const g of this.groups) {
      const style = p.styles.get(g.key);
      const visible = style && style.visible && style.presence > 0.01 && p.labelsVisible !== false;
      const off = frameOffset(g.offsets, p.frame);
      const opacityBase = visible ? Math.min(1, style.opacityScale * style.presence) : 0;
      const sizeScale = style ? style.sizeScale || 1 : 1;
      for (let i = 0; i < g.count; i++) {
        const n = g.nodes[i];
        if (!visible) { if (n.style.opacity !== "0") n.style.opacity = "0"; continue; }
        const pos = cpuFramePosition(g.positions, g.posFrames, g.count, i, p.frame, tmp);
        const s = pos && cam.project([pos[0] + off[0], pos[1] + off[1], pos[2] + off[2]], scr);
        const w = (n._w || (n._w = n.offsetWidth || 60)) * sizeScale;
        const h = (n._h || (n._h = n.offsetHeight || 16)) * sizeScale;
        if (!s || s[0] < -w || s[0] > W + w || s[1] < 0 || s[1] > H || !fits(s[0] - w / 2 - 4, s[1] - h / 2 - 2, w + 8, h + 4)) {
          n.style.opacity = "0";
          continue;
        }
        // Soften labels as they approach the viewport edge.
        const edge = Math.min(s[0], W - s[0], s[1], H - s[1]);
        n.style.opacity = String(opacityBase * Math.min(1, Math.max(0, edge / 48)));
        n.style.transform = `translate3d(${s[0].toFixed(1)}px, ${s[1].toFixed(1)}px, 0) translate(-50%, -50%) scale(${sizeScale})`;
      }
    }
  }
}
