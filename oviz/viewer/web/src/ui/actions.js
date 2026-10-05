// Actions (classic Oviz, declared with make_plot(actions=...)): buttons in an
// action bar that each run a short sequence of steps: go to a saved State,
// crossfade to a legend group, move the camera onto a group, a trace or a
// point (optionally orbiting), or play time until a stop condition. Steps
// run one after another ("after_previous") or together ("with_previous"),
// each after its delay. Clicking the running action again returns to the
// view it started from; any camera input interrupts it.

import { h, icon } from "./dom.js";
import { captureState, easing } from "../app/states.js";
import { applyGroupDefaults, cloneState } from "../app/state.js";
import { clonePose, poseFromEyeTarget } from "../engine/camera.js";
import { clamp } from "../core/math.js";

/** Group steps into phases: a phase starts at every "after_previous" step. */
export function actionPhases(steps) {
  const phases = [];
  for (const [i, step] of (steps || []).entries()) {
    if (i === 0 || step.start !== "with_previous") phases.push([step]);
    else phases[phases.length - 1].push(step);
  }
  return phases;
}

export class ActionsPlugin {
  constructor() {
    this.name = "actions";
    this.running = null; // { key, before, cancelled, timers, stopOrbit }
  }

  attach(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.items = Array.isArray(v.manifest.actions) ? v.manifest.actions.filter((a) => a && a.key && Array.isArray(a.steps)) : [];
    if (!this.items.length) return;
    this.bar = h("div", { class: "ov-actions ov-glass ov-chrome", role: "toolbar", "aria-label": "Actions" });
    this.buttons = new Map();
    for (const a of this.items) {
      const b = h("button", { class: "ov-action", type: "button", "aria-pressed": "false", title: a.description || a.label, onclick: () => this.toggle(a.key) }, a.label);
      this.buttons.set(a.key, b);
      this.bar.append(b);
    }
    ui.ui.append(this.bar);
    ui.root.dataset.actions = "true";
    // Any camera input (drag, wheel, keys, Home) interrupts a running action.
    v.on("user-camera", () => this.interrupt());
  }

  commands() {
    return (this.items || []).map((a) => ({ title: `Action: ${a.label}`, sub: a.description || "", icon: "play", keywords: "action tour step", run: () => this.run(a.key) }));
  }

  toggle(key) {
    if (this.running?.key === key) this.restore();
    else this.run(key);
  }

  sync() {
    for (const [key, b] of this.buttons || []) b.setAttribute("aria-pressed", String(this.running?.key === key));
  }

  /** Stop the running action where it is (its state stays). */
  interrupt() {
    const r = this.running;
    if (!r) return;
    r.cancelled = true;
    r.timers.forEach(clearTimeout);
    r.onCancel.forEach((fn) => fn());
    if (r.orbitOwned && !r.orbitPersist) this.ui.setAutoOrbit(false);
    this.running = null;
    this.sync();
  }

  /** Back to the view the action started from (classic: click it again). */
  restore() {
    const r = this.running;
    if (!r) return;
    this.interrupt();
    this.story()?.applyState(r.before, { duration_ms: 900 });
  }

  story() {
    return this.ui.plugins.find((p) => p.name === "states") || null;
  }

  async run(key) {
    const a = this.items.find((x) => x.key === key);
    if (!a) return;
    const previous = this.running?.before;
    this.interrupt();
    const v = this.viewer;
    const r = { key, before: previous || captureState(v, this.ui.sky), cancelled: false, timers: [], onCancel: [], orbitOwned: false, orbitPersist: true };
    this.running = r;
    this.sync();
    for (const phase of actionPhases(a.steps)) {
      if (r.cancelled) return;
      await Promise.all(phase.map((step) => this.step(step, r)));
    }
    if (this.running === r && !r.orbitOwned) { this.running = null; this.sync(); }
  }

  wait(ms, r) {
    return new Promise((resolve) => {
      if (!(ms > 0)) { resolve(); return; }
      r.timers.push(setTimeout(resolve, ms));
    });
  }

  async step(step, r) {
    await this.wait(step.delay_ms, r);
    if (r.cancelled) return;
    const v = this.viewer;
    if (step.type === "state") {
      const s = this.story();
      if (!s) return;
      const target = step.state;
      if (target === "original" || target === 0) await s.goHome();
      else {
        const i = s.indexOf(typeof target === "number" ? target : String(target));
        if (i >= 0) await s.goTo(i);
      }
    } else if (step.type === "legend_group") {
      // A smooth crossfade to the group's legend preset, camera untouched.
      const st = cloneState(captureState(v, this.ui.sky));
      applyGroupDefaults(st, v.manifest, String(step.group));
      await this.story()?.applyState(st, { duration_ms: step.duration_ms ?? 500, easing: step.easing, keepCamera: true });
      if (v.state.group !== step.group) this.ui.layers.setGroup(String(step.group));
      this.ui.layers.render();
    } else if (step.type === "camera") {
      await this.camera(step, r);
    } else if (step.type === "time") {
      await this.time(step, r);
    }
  }

  /** Objects a camera target covers: a legend group's visible traces, one trace, or none (a point). */
  targetPoints(target) {
    const v = this.viewer;
    const keys = target.kind === "trace" ? [target.key || v.traces.find((t) => t.name === target.name)?.key]
      : target.kind === "group" ? v.traces.filter((t) => t.points && t.showInLegend && (v.manifest.groups?.visibility?.[target.name]?.[t.key] === true)).map((t) => t.key)
        : [];
    const pts = [];
    for (const key of keys.filter(Boolean)) {
      const t = v.traceByKey.get(key);
      for (let i = 0; i < (t?.points?.count || 0); i++) {
        const p = v.objectPosition(key, i);
        if (p) pts.push(p);
      }
    }
    return pts;
  }

  async camera(step, r) {
    const v = this.viewer;
    let pose = clonePose(v.pose);
    const cam = step.camera || null;
    if (cam?.position && cam?.target) {
      pose = poseFromEyeTarget([cam.position.x, cam.position.y, cam.position.z], [cam.target.x, cam.target.y, cam.target.z], Number(cam.fov) || pose.fov);
    } else if (step.target?.kind === "point") {
      pose.target = [step.target.x, step.target.y, step.target.z];
      pose.distance *= step.distance_scale || 1;
    } else {
      const pts = this.targetPoints(step.target || {});
      if (pts.length) {
        const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
        for (const p of pts) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); }
        pose.target = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
        const radius = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2;
        pose.distance = Math.max(40, (radius * (step.fit_padding || 1.15)) / Math.tan((pose.fov * Math.PI) / 360)) * (step.distance_scale || 1);
      }
    }
    if (Number.isFinite(cam?.fov)) pose.fov = clamp(cam.fov, 1, 120);
    const res = await v.animateTo(pose, { duration: step.duration_ms ?? 1200, ease: easing(step.easing) });
    if (r.cancelled || !res?.done) return;
    const anchor = step.anchor_target;
    if (anchor?.kind === "trace") v.setAnchor({ kind: "layer", trace: anchor.key || v.traces.find((t) => t.name === anchor.name)?.key }, { fly: false });
    const orbit = step.orbit;
    if (orbit?.enabled) {
      const g = v.state.global;
      g.autoOrbitSpeed = orbit.speed_multiplier || 1;
      g.autoOrbitDirection = orbit.direction < 0 ? -1 : 1;
      this.ui.setAutoOrbit(true);
      r.orbitOwned = true;
      r.orbitPersist = orbit.persist !== false;
    }
  }

  time(step, r) {
    const v = this.viewer;
    const tl = v.timeline;
    const dir = step.direction === "backward" ? -1 : 1;
    const base = v.manifest.time?.playbackIntervalMs || 240;
    tl.speed = clamp(base / Math.max(step.interval_ms || base, 1), 0.05, 20);
    tl.play(dir);
    return new Promise((resolve) => {
      const start = performance.now();
      const f0 = tl.frame;
      const stop = () => { off(); tl.pause(); resolve(); };
      const off = tl.on(() => {
        if (r.cancelled) { off(); resolve(); return; }
        const t = tl.time;
        if (Number.isFinite(step.stop_at_time_myr) && (dir > 0 ? t >= step.stop_at_time_myr : t <= step.stop_at_time_myr)) { tl.setTime(step.stop_at_time_myr); stop(); }
        else if (Number.isFinite(step.stop_after_frames) && Math.abs(tl.frame - f0) >= step.stop_after_frames) stop();
        else if (Number.isFinite(step.stop_after_ms) && performance.now() - start >= step.stop_after_ms) stop();
        else if (!tl.playing) { off(); resolve(); }
      });
      r.onCancel.push(() => { off(); tl.pause(); resolve(); });
    });
  }
}
