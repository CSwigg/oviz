// Selection reticle: a lock-on ring around the selected object, tied to the
// inspector card by a leader line. Positions are recomputed after every
// rendered frame, so the mark tracks the object through camera moves and
// time playback.

const SVGNS = "http://www.w3.org/2000/svg";
const R = 17;

function el(tag, attrs = {}) {
  const e = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  return e;
}

function arc(r, a0, a1) {
  const p = (a) => [Math.cos((a * Math.PI) / 180) * r, Math.sin((a * Math.PI) / 180) * r];
  const [x0, y0] = p(a0);
  const [x1, y1] = p(a1);
  return `M${x0.toFixed(2)} ${y0.toFixed(2)}A${r} ${r} 0 0 1 ${x1.toFixed(2)} ${y1.toFixed(2)}`;
}

export class SelectionReticle {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.svg = el("svg", { class: "ov-reticle", "aria-hidden": "true" });
    this.mark = el("g", { class: "ov-reticle-mark", opacity: 0 });
    this.leaderGroup = el("g");
    this.svg.append(this.leaderGroup, this.mark);
    ui.ui.prepend(this.svg);
    this.scr = [0, 0, 0];
    this.visible = false;
    v.renderer.afterRender.push(() => this.update());
    v.on("select", () => this.onSelect());
    v.on("viewmode", () => this.update());
    window.addEventListener("resize", () => this.update());
  }

  onSelect() {
    const sel = this.ui.selection;
    this.mark.replaceChildren();
    this.leaderGroup.replaceChildren();
    this.leader = null;
    if (!sel) { this.hide(); return; }
    // Rebuilt per selection so the lock-on and draw-in animations replay.
    const lock = el("g", { class: "ov-reticle-lock" });
    const spin = el("g", { class: "ov-reticle-spin" });
    spin.append(
      el("circle", { r: R, stroke: "none" }),
      el("path", { class: "ov-reticle-ring", d: arc(R, -58, 18) }),
      el("path", { class: "ov-reticle-ring", d: arc(R, 122, 198) }),
    );
    lock.append(el("circle", { class: "ov-reticle-halo", r: R + 6 }), spin);
    for (const a of [45, 135, 225, 315]) {
      const c = Math.cos((a * Math.PI) / 180), s = Math.sin((a * Math.PI) / 180);
      lock.append(el("path", { class: "ov-reticle-tick", d: `M${(c * (R + 3)).toFixed(2)} ${(s * (R + 3)).toFixed(2)}L${(c * (R + 8)).toFixed(2)} ${(s * (R + 8)).toFixed(2)}` }));
    }
    this.mark.append(lock);
    this.leader = el("path", { class: "ov-reticle-leader", pathLength: 1 });
    this.end = el("circle", { class: "ov-reticle-end", r: 2.4 });
    this.leaderGroup.append(this.leader, this.end);
    // Follow the inspector while it swings open.
    const until = performance.now() + 700;
    const follow = () => {
      this.update();
      if (performance.now() < until && this.ui.selection === sel) requestAnimationFrame(follow);
    };
    requestAnimationFrame(follow);
  }

  hide() {
    if (!this.visible) return;
    this.visible = false;
    this.mark.setAttribute("opacity", "0");
    this.leaderGroup.setAttribute("opacity", "0");
    delete this.ui.root.dataset.reticle;
  }

  update() {
    const sel = this.ui.selection;
    const v = this.viewer;
    if (!sel) return this.hide();
    const p = sel.kind === "member" ? sel.member.position : v.objectPosition(sel.trace, sel.index);
    const cam = v.renderer.camera;
    const s = p && cam.project(p, this.scr);
    const W = cam.width, H = cam.height;
    if (!s || s[0] < -40 || s[1] < -40 || s[0] > W + 40 || s[1] > H + 40) return this.hide();
    // The canvas may be inset in the page (the Atlas frame): work in page space.
    const cr = v.canvas.getBoundingClientRect();
    const rr = this.ui.root.getBoundingClientRect();
    const x = s[0] + cr.left - rr.left;
    const y = s[1] + cr.top - rr.top;
    this.mark.setAttribute("transform", `translate(${x.toFixed(1)} ${y.toFixed(1)})`);
    this.mark.setAttribute("opacity", "1");
    this.visible = true;
    this.ui.root.dataset.reticle = "true";
    this.updateLeader(x, y);
  }

  updateLeader(x, y) {
    const ui = this.ui;
    const panel = ui.inspector.el;
    if (!this.leader || ui.narrow || panel.dataset.open !== "true" || ui.root.dataset.zen === "true") {
      this.leaderGroup.setAttribute("opacity", "0");
      return;
    }
    const root = ui.root.getBoundingClientRect();
    const r = panel.getBoundingClientRect();
    const ex = r.left - root.left - 7;
    const ey = Math.min(r.top - root.top + 40, r.bottom - root.top - 12);
    if (ex - x < R + 40) {
      this.leaderGroup.setAttribute("opacity", "0");
      return;
    }
    // Run level out of the ring, angle toward the card only near it (so the
    // line stays off the data), then enter the card level.
    const sx = x + R + 7, sy = y;
    const mx = ex - 22;
    const kx = Math.max(sx, mx - Math.abs(ey - y));
    this.leader.setAttribute("d", `M${sx.toFixed(1)} ${sy.toFixed(1)}L${kx.toFixed(1)} ${sy.toFixed(1)}L${mx.toFixed(1)} ${ey.toFixed(1)}L${ex.toFixed(1)} ${ey.toFixed(1)}`);
    this.end.setAttribute("cx", ex.toFixed(1));
    this.end.setAttribute("cy", ey.toFixed(1));
    this.leaderGroup.setAttribute("opacity", "1");
  }
}
