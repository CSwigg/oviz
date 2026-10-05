// States: capture, exact apply, smooth transitions, legacy import.
//
// A State captures the whole viewer (time, camera + view mode, layer
// visibility/styles, volumes, sky stack, global display settings). Applying
// one animates every continuous property together and then assigns the
// target exactly, so the arrival is bit-for-bit the captured view.

import { clonePose, lerpPose, makePose, sanitizePose } from "../engine/camera.js";
import { normalizeAnchor } from "../engine/anchor.js";
import { cloneJson, legacySnapshotToState, STATE_VERSION, autoOrbitRate } from "./state.js";
import { clamp, lerp, easeInOutCubic } from "../core/math.js";

const EASINGS = {
  linear: (t) => t,
  easeInOutCubic,
  easeInOutQuad: (t) => (t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2),
  easeOutCubic: (t) => 1 - Math.pow(1 - t, 3),
};

export function uid() {
  return `st-${Math.random().toString(36).slice(2, 8)}${Date.now().toString(36).slice(-4)}`;
}

/** Snapshot the live viewer. */
export function captureState(viewer, sky) {
  const s = cloneJson(viewer.state);
  s.version = STATE_VERSION;
  s.view = { mode: viewer.state.view.mode, pose: clonePose(viewer.pose), anchor: cloneJson(viewer.state.view.anchor || { kind: "free" }) };
  // A Sky view remembers the 3D view it came from, so V returns there.
  if (s.view.mode === "sky" && viewer.returnPose) s.view.returnPose = clonePose(viewer.returnPose);
  const tl = viewer.timeline;
  s.time = { frame: tl.frame, speed: tl.speed };
  // Playing time is part of the view (classic playback_state): it resumes on arrival.
  if (tl.playing) { s.time.playing = true; s.time.direction = tl.direction < 0 ? -1 : 1; }
  if (sky) s.sky = sky.captureState();
  // Plugins (e.g. the distribution filter) contribute their own state.
  s.ext = {};
  for (const [name, ext] of viewer.stateExtensions || []) s.ext[name] = ext.capture();
  return s;
}

/** A transition spec, accepting `duration` as an alias (as Python does). */
export function readTransition(t, fallback = null) {
  if (!t || typeof t !== "object") return fallback;
  const ms = Number(t.duration_ms ?? t.duration);
  return {
    duration_ms: clamp(Number.isFinite(ms) ? ms : fallback?.duration_ms ?? 1200, 0, 60000),
    easing: EASINGS[t.easing] ? String(t.easing) : fallback?.easing ?? "easeInOutCubic",
  };
}

/**
 * Read embedded States (v2 or legacy) into v2 items. Accepts the same
 * aliases as the Python normaliser (oviz/threejs_states.py): `ordered_states`,
 * `preserve_camera_on_navigation`, `camera_behavior`/`cameraBehavior`,
 * `state` holding a legacy snapshot, and `transition.duration`.
 */
export function importStates(manifest, baseState) {
  const block = manifest.states || {};
  const raw = block.items ?? block.ordered_states;
  const items = Array.isArray(raw) ? raw : [];
  const assets = block.assets && typeof block.assets === "object" ? block.assets : {};
  const defaultCamera = block.preserve_camera_on_navigation ? "keep" : "follow";
  const out = [];
  for (const item of items) {
    if (!item || typeof item !== "object") continue;
    let state = null;
    const snap = item.snapshot ?? item.state;
    if (snap && snap.version === STATE_VERSION) state = cloneJson(snap);
    else if (snap && typeof snap === "object") state = legacySnapshotToState(snap, manifest, baseState);
    if (!state) continue;
    const cam = String(item.camera_behavior ?? item.cameraBehavior ?? item.camera ?? "").trim().toLowerCase();
    const thumb = typeof item.thumb === "string" ? item.thumb : assets[item.thumb?.__oviz_asset_ref__];
    out.push({
      id: String(item.id || uid()),
      name: String(item.name || `View ${out.length + 1}`),
      caption: String(item.caption || ""),
      thumb: typeof thumb === "string" ? thumb : "",
      camera: cam === "keep" || cam === "follow" ? cam : defaultCamera,
      transition: readTransition(item.transition),
      state,
    });
  }
  return {
    items: out,
    presentOnly: !!block.present_only,
    defaultMode: block.default_mode === "present" ? "present" : "edit",
    defaultTransition: readTransition(block.default_transition, { duration_ms: 1200, easing: "easeInOutCubic" }),
    projectId: String(block.project_id || ""),
    autosave: block.autosave_drafts !== false,
  };
}

export function exportStatesBlock(project, { presentOnly = false, defaultMode = "edit" } = {}) {
  return {
    schema_version: 2,
    project_id: project.projectId,
    present_only: !!presentOnly,
    default_mode: defaultMode,
    autosave_drafts: project.autosave !== false,
    default_transition: project.defaultTransition,
    items: project.items.map((it) => ({
      id: it.id,
      name: it.name,
      caption: it.caption || "",
      thumb: it.thumb || "",
      camera: it.camera,
      transition: it.transition,
      state: it.state,
    })),
  };
}

// ---------------------------------------------------------------- transitions

function resolvedTraceOpacity(viewer, st, key) {
  const t = viewer.traceByKey.get(key);
  const ts = st.traces?.[key] || {};
  return ts.opacity ?? t?.opacity ?? 1;
}

/**
 * Build an interpolator from the live state to `target`. Returns
 * {step(t, raw) → applies the blended state (t eased, raw the linear
 * progress), finish() → applies target exactly}.
 */
export function makeTransition(viewer, sky, target, { keepCamera = false } = {}) {
  // A camera flight or fling still in progress would override the State's
  // camera once the transition lands, so the transition takes it over.
  if (!keepCamera) {
    viewer._cancelTween?.();
    viewer.controls?.stop?.();
  }
  const from = captureState(viewer, sky);
  const to = cloneJson(target);
  if (target.view?.pose) to.view.pose = sanitizePose(clonePose(target.view.pose));
  if (keepCamera) {
    to.view = { mode: from.view.mode, pose: clonePose(viewer.pose), anchor: from.view.anchor };
  }
  // The camera anchor is part of the camera. States saved before anchors
  // existed get the one their pose implies (the LSR if they orbit it).
  const anchorB = normalizeAnchor(to.view?.anchor)
    || (to.view?.pose && viewer.inferAnchor ? viewer.inferAnchor(to.view.pose) : null)
    || from.view.anchor;
  const traceKeys = new Set([...Object.keys(from.traces || {}), ...Object.keys(to.traces || {})]);
  const volKeys = new Set([...Object.keys(from.volumes || {}), ...Object.keys(to.volumes || {})]);
  const poseA = clonePose(from.view.pose || viewer.pose);
  const poseB = to.view?.pose ? clonePose(to.view.pose) : clonePose(poseA);
  const modeChange = from.view.mode !== to.view.mode;
  // Arc the camera out when the view jumps far, like a real camera move.
  const jump = Math.hypot(poseB.target[0] - poseA.target[0], poseB.target[1] - poseA.target[1], poseB.target[2] - poseA.target[2]);
  const lift = modeChange ? 0 : Math.min(jump * 0.3, Math.max(poseA.distance, poseB.distance) * 1.2);
  const tmpPose = makePose();
  const frameA = from.time.frame, frameB = to.time?.frame ?? frameA;
  const skyA = from.sky?.layers || [];
  const skyB = to.sky?.layers || [];
  const skyByKeyA = new Map(skyA.map((l) => [l.key, l]));
  const skyByKeyB = new Map(skyB.map((l) => [l.key, l]));
  const skyKeys = [...new Set([...skyB.map((l) => l.key), ...skyA.map((l) => l.key)])];
  if (modeChange && to.view.mode === "sky") viewer.controls.mode = "sky";
  // Leaving Sky: member stars fade out as the camera leaves the Sun, not
  // after it has flown all the way into 3D.
  if (modeChange && from.view.mode === "sky") sky?.revealMembers?.(false);
  // Plugins may animate their part of the change too (the lasso fades from
  // one selection to the next); their apply still lands it exactly.
  const extSteps = [];
  for (const [name, ext] of viewer.stateExtensions || []) {
    const tr = ext.transition?.(from.ext?.[name] ?? null, to.ext ? to.ext[name] ?? null : null);
    if (tr?.step) extSteps.push(tr);
  }

  function step(t, raw = t) {
    const st = viewer.state;
    // Camera ("keep" States leave the live camera alone entirely).
    if (!keepCamera) {
      lerpPose(poseA, poseB, t, tmpPose);
      if (lift > 0) tmpPose.distance += Math.sin(Math.PI * t) * lift;
      Object.assign(viewer.renderer.camera.pose, clonePose(tmpPose));
    }
    // Time
    viewer.timeline.setFrame(lerp(frameA, frameB, t), { silent: t < 1 });
    // Traces: numeric lerps; visibility changes fade through opacity.
    for (const key of traceKeys) {
      const a = from.traces?.[key] || {}, b = to.traces?.[key] || {};
      const visA = !!a.visible, visB = !!b.visible;
      const opA = resolvedTraceOpacity(viewer, from, key), opB = resolvedTraceOpacity(viewer, to, key);
      const cur = (st.traces[key] ||= {});
      const discrete = t < 0.5 ? a : b;
      cur.colorMode = discrete.colorMode;
      cur.colormap = discrete.colormap;
      cur.color = discrete.color;
      cur.inGroup = visA || visB ? true : discrete.inGroup;
      cur.sizeScale = lerp(a.sizeScale ?? 1, b.sizeScale ?? 1, t);
      if (visA && visB) { cur.visible = true; cur.opacity = lerp(opA, opB, t); }
      else if (visA) { cur.visible = t < 1; cur.opacity = opA * (1 - t); }
      else if (visB) { cur.visible = t > 0; cur.opacity = opB * t; }
      else cur.visible = false;
    }
    // Volumes
    for (const key of volKeys) {
      const a = from.volumes?.[key], b = to.volumes?.[key];
      if (!a || !b) continue;
      const cur = st.volumes[key];
      if (!cur) continue;
      const d = t < 0.5 ? a : b;
      Object.assign(cur, { stretch: d.stretch, colormap: d.colormap, showAllTimes: d.showAllTimes, lightingMode: d.lightingMode });
      if (d.kdeTimeWindowMyr != null) cur.kdeTimeWindowMyr = d.kdeTimeWindowMyr;
      for (const f of ["vmin", "vmax", "alphaCoef", "steps"]) cur[f] = lerp(Number(a[f]), Number(b[f]), t);
      if (a.visible && b.visible) { cur.visible = true; cur.opacity = lerp(a.opacity, b.opacity, t); }
      else if (a.visible) { cur.visible = t < 1; cur.opacity = a.opacity * (1 - t); }
      else if (b.visible) { cur.visible = t > 0; cur.opacity = b.opacity * t; }
      else cur.visible = false;
    }
    viewer.volumes.invalidate();
    // Globals
    const ga = from.global || {}, gb = to.global || {};
    for (const f of ["pointSize", "pointOpacity", "glow", "fadeTime", "gridOpacity", "trails"]) {
      if (Number.isFinite(ga[f]) && Number.isFinite(gb[f])) st.global[f] = lerp(ga[f], gb[f], t);
    }
    for (const f of ["fadeInOut", "fadeByOpacity", "sizeByStars", "labels"]) st.global[f] = (t < 0.5 ? ga : gb)[f];
    if (ga.grid !== gb.grid) {
      st.global.grid = true;
      st.global.gridOpacity = ga.grid ? 1 - t : t;
    } else st.global.grid = gb.grid;
    // Images
    for (const [k, b] of Object.entries(to.images || {})) {
      const a = from.images?.[k] || b;
      const cur = st.images[k] || (st.images[k] = {});
      if (a.visible !== b.visible) { cur.visible = true; cur.opacity = a.visible ? 1 - t : t; }
      else { cur.visible = b.visible; cur.opacity = lerp(a.opacity ?? 1, b.opacity ?? 1, t); }
    }
    // Sky layers: the target's layers (on top of the union stack) fade in
    // over the outgoing ones, which hold until the last stretch, so a change
    // of background is (1 − t)·A + t·B instead of dimming through black.
    if (sky && (skyA.length || skyB.length)) {
      const onA = from.sky?.backgroundVisible !== false, onB = to.sky?.backgroundVisible !== false;
      const hold = 1 - smoothstepT(0.7, 1, t);
      const layers = skyKeys.map((k) => {
        const a = skyByKeyA.get(k), b = skyByKeyB.get(k);
        const base = cloneJson(b || a);
        const oa = onA && a && a.visible !== false ? a.opacity ?? 1 : 0;
        const ob = onB && b && b.visible !== false ? b.opacity ?? 1 : 0;
        base.opacity = b ? lerp(oa, ob, t) : oa * hold;
        base.visible = base.opacity > 0.001;
        return base;
      });
      viewer.state.sky.layers = layers;
      viewer.state.sky.backgroundVisible = onA || onB;
      sky.throttledApply();
    }
    for (const tr of extSteps) tr.step(t, raw);
    viewer.renderer.markInteraction();
    viewer.invalidatePick();
  }

  function finish() {
    // Exact arrival: assign the captured target verbatim.
    const st = viewer.state;
    st.traces = cloneJson(to.traces || {});
    st.volumes = cloneJson(to.volumes || st.volumes);
    st.images = cloneJson(to.images || st.images);
    st.global = cloneJson(to.global || st.global);
    st.group = to.group ?? st.group;
    st.view.mode = to.view.mode;
    viewer.controls.mode = to.view.mode === "sky" ? "sky" : "galactic";
    if (!keepCamera) Object.assign(viewer.renderer.camera.pose, clonePose(poseB));
    st.view.pose = viewer.renderer.camera.pose;
    const anchorChanged = JSON.stringify(st.view.anchor || null) !== JSON.stringify(anchorB || null);
    st.view.anchor = cloneJson(anchorB || { kind: "free" });
    // Arrive exactly: the anchor starts riding from here, it never snaps.
    viewer._anchorLast = null;
    viewer.timeline.setFrame(frameB);
    viewer._anchorLast = viewer.anchorPoint?.() ?? null;
    if (anchorChanged) viewer.emit?.("anchor", { anchor: st.view.anchor, released: false });
    if (Number.isFinite(to.time?.speed)) viewer.timeline.speed = to.time.speed;
    if (to.view.mode === "sky" && to.view.returnPose) viewer.returnPose = clonePose(to.view.returnPose);
    if (sky && to.sky) sky.applyState(to.sky);
    for (const [name, ext] of viewer.stateExtensions || []) ext.apply(to.ext ? to.ext[name] ?? null : null);
    // Auto-orbit is camera behaviour: "keep" views leave it as it was.
    if (keepCamera) st.global.autoOrbit = !!from.global?.autoOrbit;
    viewer.controls.autoOrbit = autoOrbitRate(st.global);
    viewer.controls.speed = st.global.navSpeed ?? 1;
    // A view saved while time played plays on from its frame.
    if (to.time?.playing) viewer.timeline.play(to.time.direction < 0 ? -1 : 1);
    viewer.volumes.invalidate();
    viewer.invalidatePick();
    viewer.renderer.invalidate();
  }

  return { step, finish, modeChange, keepCamera, from, to };
}

function smoothstepT(e0, e1, x) {
  const t = clamp((x - e0) / (e1 - e0), 0, 1);
  return t * t * (3 - 2 * t);
}

export function easing(name) {
  return EASINGS[name] || easeInOutCubic;
}

/** Numeric fidelity check used by tests and the debug overlay. */
export function stateDifferences(a, b, eps = 1e-6) {
  const diffs = [];
  const walk = (x, y, path) => {
    if (typeof x === "number" && typeof y === "number") {
      if (Math.abs(x - y) > eps * Math.max(1, Math.abs(x), Math.abs(y))) diffs.push(path);
      return;
    }
    if (x && y && typeof x === "object" && typeof y === "object") {
      for (const k of new Set([...Object.keys(x), ...Object.keys(y)])) walk(x[k], y[k], `${path}.${k}`);
      return;
    }
    if (x !== y && !(x == null && y == null)) diffs.push(path);
  };
  walk(a, b, "state");
  return diffs;
}

