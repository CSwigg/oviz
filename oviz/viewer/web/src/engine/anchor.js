// Camera anchors: the point the camera orbits and zooms about, and the body
// it rides with through time.
//
//   {kind: "lsr"}                   the Local Standard of Rest. Oviz scenes
//                                   are LSR-centred (positions are relative
//                                   to a reference orbit that starts at the
//                                   Sun with the LSR's velocity), so the LSR
//                                   is the origin at every time: an anchored
//                                   camera stays with it while the Sun and
//                                   the clusters move around it.
//   {kind: "sun"}                   the Sun's track (the figure's "Sun" trace).
//   {kind: "object", trace, index}  any object's track ("Follow").
//   {kind: "free"}                  no anchor: the wheel zooms toward the
//                                   cursor and the camera stays where it is.
//
// While anchored, zooming keeps the anchor at the centre of the view, and
// the camera moves by the anchor's own displacement as time changes, so it
// never jumps. Panning or flying elsewhere frees the camera; Home re-anchors.

/** The LSR in scene coordinates (the origin of Oviz's LSR-centred frame). */
export const LSR_POINT = Object.freeze([0, 0, 0]);

const SIMPLE_KINDS = ["lsr", "sun", "free"];

/** A canonical anchor object, or null if `a` is not one (strings allowed). */
export function normalizeAnchor(a) {
  if (typeof a === "string") return decodeAnchor(a);
  if (!a || typeof a !== "object") return null;
  const kind = String(a.kind || "").toLowerCase();
  if (kind === "object") {
    const index = Number(a.index);
    if (typeof a.trace !== "string" || !a.trace || !Number.isInteger(index) || index < 0) return null;
    return { kind, trace: a.trace, index };
  }
  return SIMPLE_KINDS.includes(kind) ? { kind } : null;
}

/** Compact text form (view links): lsr, sun, free, or o:<trace>:<index>. */
export function encodeAnchor(a) {
  const n = normalizeAnchor(a);
  if (!n) return "";
  return n.kind === "object" ? `o:${n.trace}:${n.index}` : n.kind;
}

export function decodeAnchor(s) {
  const text = String(s ?? "").trim();
  if (text.startsWith("o:")) {
    const i = text.lastIndexOf(":");
    const index = Number(text.slice(i + 1));
    return i > 2 && /^\d+$/.test(text.slice(i + 1)) && Number.isSafeInteger(index)
      ? { kind: "object", trace: text.slice(2, i), index }
      : null;
  }
  const kind = text.toLowerCase();
  return SIMPLE_KINDS.includes(kind) ? { kind } : null;
}

export function sameAnchor(a, b) {
  return encodeAnchor(a) === encodeAnchor(b);
}

/**
 * The anchor a pose implies when none was recorded (States and links from
 * before anchors existed): the LSR if the camera orbits it, otherwise free.
 * A Sky pose (eye at the Sun) implies nothing, so it gets `fallback`.
 */
export function inferAnchor(pose, lsr, tolerance, fallback = { kind: "free" }) {
  if (!pose || !Array.isArray(pose.target)) return { ...fallback };
  if (!(pose.distance > 1e-6)) return { ...fallback };
  if (!lsr) return { kind: "free" };
  const t = pose.target;
  const d = Math.hypot(t[0] - lsr[0], t[1] - lsr[1], t[2] - lsr[2]);
  return d <= tolerance ? { kind: "lsr" } : { kind: "free" };
}

/** Short name for the readout and menus. */
export function anchorLabel(a, objectName = "") {
  switch (a?.kind) {
    case "lsr": return "the LSR";
    case "sun": return "the Sun";
    case "object": return objectName || "the selected object";
    default: return "";
  }
}
