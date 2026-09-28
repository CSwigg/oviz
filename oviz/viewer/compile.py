"""Compile a legacy Oviz scene spec into a columnar viewer bundle.

The legacy spec (produced by :meth:`oviz.Animate3D.make_plot` or decoded from
an existing HTML figure) stores every frame as a list of traces, and every
trace as a list of point dictionaries. That layout repeats all static
metadata once per frame and encodes numbers as decimal text.

This compiler converts it into:

* one object table per trace (stable identities matched across frames by
  ``motion.key`` → selection name → index, exactly like the legacy runtime),
* ``(frames, objects, 3)`` float32 position arrays, collapsed to a single
  frame when static or to ``base + per-frame offset`` when every vertex moves
  rigidly (the reference circles and labels re-centred each frame),
* per-object attributes (colour, age, n_stars, sky coordinates) stored once,
* raw voxel bytes and LUTs for volumes (no nested base64),
* member stars as Cartesian position/velocity columns.

Nothing scientific is recomputed: positions, sizes, opacities, colours, ages
and volume windows are carried over value-for-value from the spec.
"""

from __future__ import annotations

import base64
import io
import math
import re
from collections import OrderedDict
from typing import Any, Iterable

import numpy as np

from .bundle import CRITICAL, DEFERRED, NORMAL, BUNDLE_FORMAT, Bundle, BundleBuilder
from .colors import css_hex, parse_color

KM_S_TO_PC_MYR = 1.0227121650537077
MAS_YR_KM_S_PER_PC = 4.740470463533348e-3  # v[km/s] = this * mu[mas/yr] * d[pc]

_ICRS_TO_GAL = np.array(
    [
        [-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
        [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
        [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669],
    ]
)

SYMBOLS = ["circle", "square", "diamond", "cross", "x", "triangle-up", "star", "circle-open"]


class CompileError(ValueError):
    pass


def _num(value: Any, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def _xyz(obj: Any) -> list[float] | None:
    if isinstance(obj, dict):
        vals = [_num(obj.get("x")), _num(obj.get("y")), _num(obj.get("z"))]
    elif isinstance(obj, (list, tuple)) and len(obj) >= 3:
        vals = [_num(v) for v in obj[:3]]
    else:
        return None
    return vals if all(math.isfinite(v) for v in vals) else None


def _strip_heavy(value: Any, limit: int = 20000) -> Any:
    """Drop long strings (data URLs, base64) from free-form dictionaries."""

    if isinstance(value, dict):
        return {k: _strip_heavy(v, limit) for k, v in value.items() if not (isinstance(v, str) and len(v) > limit)}
    if isinstance(value, list):
        return [_strip_heavy(v, limit) for v in value if not (isinstance(v, str) and len(v) > limit)]
    return value


def _clean_json(value: Any) -> Any:
    """Make a value strict-JSON safe (NaN/Inf → None, tuples → lists)."""

    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): _clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_json(v) for v in value]
    if isinstance(value, (np.floating,)):
        return _clean_json(float(value))
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.ndarray):
        return _clean_json(value.tolist())
    return value


# ---------------------------------------------------------------------------
# Position channels


def _frames_equal(arr: np.ndarray, tol: float = 1e-6) -> bool:
    """True when every frame of a (F, ...) array equals frame 0 (NaN-aware)."""

    if arr.shape[0] <= 1:
        return True
    base = arr[0]
    nan0 = np.isnan(base)
    for f in range(1, arr.shape[0]):
        cur = arr[f]
        if not np.array_equal(np.isnan(cur), nan0):
            return False
        diff = np.abs(np.where(nan0, 0.0, cur - base))
        scale = np.maximum(1.0, np.abs(np.where(nan0, 0.0, base)))
        if np.any(diff > tol * scale):
            return False
    return True


def _rigid_offsets(pos: np.ndarray, tol: float = 2e-3) -> np.ndarray | None:
    """If every frame is frame-0 translated by one vector, return (F, 3)."""

    F = pos.shape[0]
    if F <= 1:
        return None
    base = pos[0]
    valid = ~np.isnan(base[:, 0])
    if valid.sum() < 2:
        return None
    offsets = np.zeros((F, 3), dtype=np.float64)
    for f in range(F):
        cur = pos[f]
        if not np.array_equal(~np.isnan(cur[:, 0]), valid):
            return None
        delta = cur[valid] - base[valid]
        mean = delta.mean(axis=0)
        # Specs are often rounded to 0.01 pc, so allow a couple of quanta.
        if np.max(np.abs(delta - mean)) > tol * 10.0 + 1e-6 * float(np.max(np.abs(base[valid]))):
            return None
        offsets[f] = mean
    if np.max(np.abs(offsets)) < 1e-9:
        return None
    return offsets


def _fill_absent_frames(arr: np.ndarray, present: list[int] | None) -> np.ndarray:
    """Copy the nearest present frame into frames where the trace is hidden.

    Those frames are never drawn, so their values only affect compression:
    a trace that exists at a single frame collapses to one static frame.
    """

    if not present or len(present) == arr.shape[0]:
        return arr
    out = arr.copy()
    pres = np.array(sorted(present))
    for f in range(arr.shape[0]):
        if f in present:
            continue
        nearest = int(pres[np.argmin(np.abs(pres - f))])
        out[f] = arr[nearest]
    return out


def _position_channel(builder: BundleBuilder, pos: np.ndarray, hint: str, priority: int,
                      allow_rigid: bool = False) -> dict[str, Any]:
    """Encode (F, N, 3) positions as static / rigid / per-frame.

    Static collapse only happens when frames agree to float32 precision.
    The rigid (base + per-frame offset) form reproduces positions to within
    the spec's 0.01 pc rounding, so it is reserved for decorative geometry
    (reference lines and labels) and never used for data points.
    """

    if _frames_equal(pos):
        return {"blob": builder.array(pos[:1], "f32", hint=hint, priority=priority), "frames": 1}
    offsets = _rigid_offsets(pos) if allow_rigid else None
    if offsets is not None:
        return {
            "blob": builder.array(pos[:1], "f32", hint=hint, priority=priority),
            "frames": 1,
            "offset": builder.array(offsets, "f32", hint=hint + "-offset", priority=priority),
        }
    return {"blob": builder.array(pos, "f32", hint=hint, priority=priority), "frames": int(pos.shape[0])}


def _scalar_channel(
    builder: BundleBuilder,
    values: np.ndarray,
    hint: str,
    priority: int,
    default: float | None = None,
) -> dict[str, Any] | None:
    """Encode (F, N) scalars; NaN means "use default"."""

    if np.all(np.isnan(values)):
        return None
    if _frames_equal(values):
        row = values[0]
        finite = row[~np.isnan(row)]
        if finite.size and np.all(finite == finite[0]) and not np.any(np.isnan(row)):
            if default is not None and abs(float(finite[0]) - default) < 1e-9:
                return {"value": float(default)}
            return {"value": float(finite[0])}
        return {"blob": builder.array(row, "f32", hint=hint, priority=priority), "frames": 1}
    return {"blob": builder.array(values, "f32", hint=hint, priority=priority), "frames": int(values.shape[0])}


# ---------------------------------------------------------------------------
# Traces


def _point_identity(point: dict[str, Any]) -> str:
    motion = point.get("motion")
    if isinstance(motion, dict) and str(motion.get("key") or "").strip():
        return "m:" + str(motion.get("key")).strip()
    selection = point.get("selection")
    if isinstance(selection, dict):
        name = str(selection.get("cluster_name") or "").strip()
        if name:
            return "s:" + name
    return ""


def _stars_factor(n_stars: np.ndarray) -> np.ndarray:
    valid = np.isfinite(n_stars) & (n_stars > 0)
    out = np.ones_like(n_stars, dtype=np.float64)
    if not np.any(valid):
        return out
    median = float(np.median(n_stars[valid]))
    median = max(median if math.isfinite(median) else 1.0, 1e-3)
    out[valid] = np.clip(np.sqrt(n_stars[valid] / median), 0.35, 2.75)
    return out


class _TraceCompiler:
    def __init__(self, builder: BundleBuilder, times: np.ndarray, t0_index: int, colormaps: "_ColormapRegistry"):
        self.builder = builder
        self.times = times
        self.F = len(times)
        self.t0 = t0_index
        self.colormaps = colormaps

    def compile(self, key: str, entries: list[dict | None], legend_item: dict | None) -> dict[str, Any] | None:
        F = self.F
        present = [e for e in entries if isinstance(e, dict)]
        if not present:
            return None
        ref = entries[self.t0] if isinstance(entries[self.t0], dict) else present[-1]
        legend_item = legend_item or {}

        def pick(field: str, default: Any = None) -> Any:
            for src in (ref, legend_item, *present):
                if isinstance(src, dict) and src.get(field) is not None:
                    return src.get(field)
            return default

        name = str(pick("name", key))
        default_color = pick("default_color") or pick("legend_color") or pick("color") or "#ffffff"
        spec: dict[str, Any] = {
            "key": key,
            "name": name,
            "showInLegend": bool(legend_item) and legend_item.get("kind", "trace") == "trace",
            "legendColor": css_hex(pick("legend_color") or default_color),
            "color": css_hex(default_color),
            "opacity": _num(pick("default_opacity", pick("opacity", 1.0)), 1.0),
            "pointSize": _num(pick("default_point_size", 0.0), 0.0),
            "symbol": str(pick("default_symbol", "circle") or "circle"),
            "sizeByStarsDefault": bool(pick("size_by_n_stars_default", False)),
            "hasStars": bool(pick("has_n_stars", False)),
        }
        color_by = legend_item.get("color_by") or pick("color_by")
        if isinstance(color_by, dict) and color_by.get("mode"):
            options = color_by.get("colormap_options") or []
            names = [self.colormaps.add(o) for o in options if isinstance(o, dict)]
            spec["colorBy"] = {
                "mode": str(color_by.get("mode")),
                "label": str(color_by.get("label") or color_by.get("mode")),
                "cmin": _num(color_by.get("cmin"), 0.0),
                "cmax": _num(color_by.get("cmax"), 1.0),
                "colormap": str(color_by.get("colormap") or (names[0] if names else "turbo")),
                "colormaps": [n for n in names if n],
                # Legacy default: colour by value unless the spec says "fixed".
                "defaultMode": "fixed" if str(color_by.get("default_color_mode") or "by_value") == "fixed" else "by_value",
            }

        presence = [
            i for i, e in enumerate(entries)
            if isinstance(e, dict) and (e.get("points") or e.get("segments") or e.get("labels"))
        ]
        if len(presence) not in (0, F):
            spec["presence"] = presence

        has_points = any(e.get("points") for e in present)
        has_segments = any(e.get("segments") for e in present)
        has_labels = any(e.get("labels") for e in present)
        if has_points:
            spec["points"] = self._points(key, entries, spec)
        if has_segments:
            spec["lines"] = self._lines(key, entries, ref, spec)
        if has_labels:
            spec["labels"] = self._labels(key, entries)
        if not (has_points or has_segments or has_labels):
            return None
        return spec

    # -- points --------------------------------------------------------------

    def _points(self, key: str, entries: list[dict | None], spec: dict[str, Any]) -> dict[str, Any]:
        F = self.F
        ids: "OrderedDict[str, int]" = OrderedDict()
        per_frame_index: list[list[int]] = []
        for e in entries:
            idx_list: list[int] = []
            if isinstance(e, dict):
                occurrences: dict[str, int] = {}
                for i, p in enumerate(e.get("points") or []):
                    ident = _point_identity(p) if isinstance(p, dict) else ""
                    if ident:
                        n = occurrences.get(ident, 0)
                        occurrences[ident] = n + 1
                        ident = f"{ident}#{n}"
                    else:
                        ident = f"i:{i}"
                    if ident not in ids:
                        ids[ident] = len(ids)
                    idx_list.append(ids[ident])
            per_frame_index.append(idx_list)
        N = len(ids)
        pos = np.full((F, N, 3), np.nan, dtype=np.float64)
        size = np.full((F, N), np.nan)
        opacity = np.full((F, N), np.nan)
        scalar = np.full((F, N), np.nan)
        rgba = np.full((N, 4), np.nan)
        symbol = np.zeros(N, dtype=np.uint8)
        age_now = np.full(N, np.nan)
        n_stars = np.full(N, np.nan)
        info: list[dict[str, Any] | None] = [None] * N
        info_frame = np.full(N, -1)
        default_size = spec["pointSize"]
        default_opacity = spec["opacity"]
        sym_index = {s: i for i, s in enumerate(SYMBOLS)}
        color_cache: dict[Any, tuple] = {}

        # Prefer metadata from the frame closest to t = 0.
        frame_order = sorted(range(F), key=lambda f: abs(f - self.t0))
        rank = {f: r for r, f in enumerate(frame_order)}

        for f, e in enumerate(entries):
            if not isinstance(e, dict):
                continue
            pts = e.get("points") or []
            idx = per_frame_index[f]
            for i, p in enumerate(pts):
                if not isinstance(p, dict):
                    continue
                j = idx[i]
                pos[f, j, 0] = _num(p.get("x"))
                pos[f, j, 1] = _num(p.get("y"))
                pos[f, j, 2] = _num(p.get("z"))
                size[f, j] = _num(p.get("size"), default_size)
                opacity[f, j] = _num(p.get("opacity"), default_opacity)
                cs = p.get("color_scalar")
                if cs is not None:
                    scalar[f, j] = _num(cs)
                col = p.get("color")
                if col is not None and np.isnan(rgba[j, 0]):
                    c = color_cache.get(col)
                    if c is None:
                        c = parse_color(col)
                        color_cache[col] = c
                    rgba[j] = c
                sym = p.get("symbol")
                if sym:
                    symbol[j] = sym_index.get(str(sym), 0)
                motion = p.get("motion") if isinstance(p.get("motion"), dict) else None
                selection = p.get("selection") if isinstance(p.get("selection"), dict) else None
                if motion is not None and np.isnan(age_now[j]):
                    age_now[j] = _num(motion.get("age_now_myr"))
                if np.isnan(n_stars[j]):
                    n_stars[j] = _num(p.get("n_stars"), _num((selection or {}).get("n_stars")))
                has_meta = selection is not None or bool(p.get("hovertext")) or motion is not None
                if has_meta and (info_frame[j] < 0 or rank[f] < rank.get(int(info_frame[j]), 1 << 30)):
                    if info_frame[j] < 0 or selection is not None or info[j] is None or not info[j].get("sel"):
                        info[j] = {"sel": selection, "motion": motion, "hover": p.get("hovertext")}
                        info_frame[j] = f

        out: dict[str, Any] = {"count": N}
        b = self.builder
        present = spec.get("presence")
        pos = _fill_absent_frames(pos, present)
        size = _fill_absent_frames(size, present)
        opacity = _fill_absent_frames(opacity, present)
        scalar = _fill_absent_frames(scalar, present)
        out["position"] = _position_channel(b, pos, f"{key}-pos", CRITICAL)

        size_ch = _scalar_channel(b, size, f"{key}-size", CRITICAL, default_size)
        out["size"] = size_ch or {"value": default_size}
        if "blob" in out["size"] and out["size"].get("frames", 1) > 1:
            size_max = np.nanmax(np.where(np.isnan(size), -np.inf, size), axis=0)
            size_max[~np.isfinite(size_max)] = 0.0
            out["sizeMax"] = {"blob": b.array(size_max, "f32", hint=f"{key}-sizemax", priority=CRITICAL)}

        op_ch = _scalar_channel(b, opacity, f"{key}-opacity", CRITICAL, default_opacity)
        out["opacity"] = op_ch or {"value": default_opacity}

        sc_ch = _scalar_channel(b, scalar, f"{key}-scalar", CRITICAL)
        if sc_ch is not None:
            out["colorScalar"] = sc_ch

        if not np.all(np.isnan(rgba[:, 0])):
            base = parse_color(spec["color"])
            fill = np.where(np.isnan(rgba), np.array(base)[None, :], rgba)
            rgb8 = np.clip(np.round(fill[:, :3] * 255), 0, 255).astype(np.uint8)
            if not np.all(rgb8 == rgb8[0]):
                out["rgb"] = {"blob": b.array(rgb8, "u8", hint=f"{key}-rgb", priority=CRITICAL)}
            else:
                spec["color"] = "#{:02x}{:02x}{:02x}".format(*[int(v) for v in rgb8[0]])
            # Colour alpha is deliberately ignored: the scene builder already
            # folds it into each point's opacity, and the classic runtime
            # never re-applied it.
        if np.any(symbol):
            out["symbol"] = {"blob": b.array(symbol, "u8", hint=f"{key}-symbol", priority=CRITICAL)}
            out["symbols"] = SYMBOLS
        if not np.all(np.isnan(age_now)):
            out["ageNow"] = {"blob": b.array(np.where(np.isnan(age_now), np.nan, age_now), "f32", hint=f"{key}-age", priority=CRITICAL)}
        if spec["hasStars"] and not np.all(np.isnan(n_stars)):
            out["nStars"] = {"blob": b.array(np.nan_to_num(n_stars, nan=-1.0), "f32", hint=f"{key}-nstars", priority=CRITICAL)}
            out["starsFactor"] = {"blob": b.array(_stars_factor(n_stars), "f32", hint=f"{key}-stars", priority=CRITICAL)}

        out["meta"] = self._meta_columns(key, ids, info, age_now, n_stars)
        return out

    def _meta_columns(self, key, ids, info, age_now, n_stars) -> dict[str, Any]:
        N = len(ids)
        names: list[str] = []
        aliases: list[str] = []
        lon = np.full(N, np.nan)
        lat = np.full(N, np.nan)
        dist = np.full(N, np.nan)
        ra = np.full(N, np.nan)
        dec = np.full(N, np.nan)
        hover: list[str] = []
        need_hover = False
        for ident, j in ids.items():
            rec = info[j] or {}
            sel = rec.get("sel") or {}
            motion = rec.get("motion") or {}
            raw = str(sel.get("cluster_name") or motion.get("key") or "")
            if not raw and ident.startswith(("m:", "s:")):
                raw = ident[2:].rsplit("#", 1)[0]
            names.append(raw)
            aliases.append(str(sel.get("name_all") or ""))
            lon[j] = _num(sel.get("l_deg"))
            lat[j] = _num(sel.get("b_deg"))
            dist[j] = _num(sel.get("dist_pc"))
            ra[j] = _num(sel.get("ra_deg"))
            dec[j] = _num(sel.get("dec_deg"))
            h = rec.get("hover") or ""
            # Structured fields rebuild the standard Oviz tooltip; keep raw
            # hover HTML only when there is nothing structured to show.
            if h and not sel:
                need_hover = True
            hover.append(str(h) if (h and not sel) else "")
        cols: dict[str, Any] = {"name": names}
        if any(aliases):
            cols["aliases"] = aliases
        for label, arr in (("l", lon), ("b", lat), ("dist", dist), ("ra", ra), ("dec", dec)):
            if not np.all(np.isnan(arr)):
                cols[label] = [None if not math.isfinite(v) else round(float(v), 6) for v in arr]
        if need_hover:
            cols["hover"] = hover
        return {"blob": self.builder.json(cols, hint=f"{key}-meta", priority=NORMAL)}

    # -- lines ---------------------------------------------------------------

    def _lines(self, key: str, entries: list[dict | None], ref: dict, spec: dict[str, Any]) -> dict[str, Any]:
        F = self.F
        counts = [len(e.get("segments") or []) if isinstance(e, dict) else 0 for e in entries]
        S = max(counts)
        pos = np.full((F, 2 * S, 3), np.nan)
        for f, e in enumerate(entries):
            if not isinstance(e, dict):
                continue
            segs = e.get("segments") or []
            for s, seg in enumerate(segs):
                if isinstance(seg, (list, tuple)) and len(seg) >= 6:
                    pos[f, 2 * s] = [_num(seg[0]), _num(seg[1]), _num(seg[2])]
                    pos[f, 2 * s + 1] = [_num(seg[3]), _num(seg[4]), _num(seg[5])]
        # Arc length at the reference frame (for dash patterns).
        fref = self.t0 if counts[self.t0] else int(np.argmax(counts))
        p = pos[fref]
        arc = np.zeros(2 * S)
        acc = 0.0
        for s in range(S):
            a, b_ = p[2 * s], p[2 * s + 1]
            if s > 0 and not np.allclose(a, p[2 * s - 1], equal_nan=False):
                acc += 0.0
            arc[2 * s] = acc
            seg_len = float(np.linalg.norm(b_ - a)) if np.all(np.isfinite(a)) and np.all(np.isfinite(b_)) else 0.0
            acc += seg_len
            arc[2 * s + 1] = acc
        line = ref.get("line") if isinstance(ref.get("line"), dict) else {}
        color = line.get("color") or spec["color"]
        r, g, b_, a = parse_color(color)
        opacity = _num(ref.get("opacity"), _num(ref.get("default_opacity"), 1.0)) * a
        pos = _fill_absent_frames(pos, spec.get("presence"))
        return {
            "count": S,
            "position": _position_channel(self.builder, pos, f"{key}-lines", CRITICAL, allow_rigid=True),
            "arc": {"blob": self.builder.array(arc, "f32", hint=f"{key}-arc", priority=CRITICAL)},
            "color": css_hex(color),
            "width": _num(line.get("width"), 1.0),
            "dash": str(line.get("dash") or "solid"),
            "opacity": opacity,
        }

    # -- labels ---------------------------------------------------------------

    def _labels(self, key: str, entries: list[dict | None]) -> dict[str, Any]:
        F = self.F
        ids: "OrderedDict[str, int]" = OrderedDict()
        styles: dict[int, dict[str, Any]] = {}
        per_frame: list[list[tuple[int, dict]]] = []
        for e in entries:
            rows = []
            if isinstance(e, dict):
                seen: dict[str, int] = {}
                for lab in e.get("labels") or []:
                    if not isinstance(lab, dict):
                        continue
                    text = str(lab.get("text") or "")
                    n = seen.get(text, 0)
                    seen[text] = n + 1
                    ident = f"{text}#{n}"
                    if ident not in ids:
                        ids[ident] = len(ids)
                    j = ids[ident]
                    rows.append((j, lab))
                    styles.setdefault(j, lab)
            per_frame.append(rows)
        L = len(ids)
        pos = np.full((F, L, 3), np.nan)
        for f, rows in enumerate(per_frame):
            for j, lab in rows:
                xyz = _xyz(lab)
                if xyz:
                    pos[f, j] = xyz
        present = [f for f, rows in enumerate(per_frame) if rows]
        pos = _fill_absent_frames(pos, present if len(present) != F else None)
        texts = [ident.rsplit("#", 1)[0] for ident in ids]
        style_rows = [styles.get(j, {}) for j in range(L)]
        return {
            "count": L,
            "text": texts,
            "position": _position_channel(self.builder, pos, f"{key}-labels", CRITICAL, allow_rigid=True),
            "color": [css_hex(s.get("color") or "#e2e8f0") for s in style_rows],
            "size": [_num(s.get("screen_px"), _num(s.get("size"), 14.0)) for s in style_rows],
            "screenStable": [bool(s.get("screen_stable", True)) for s in style_rows],
            "family": str((style_rows[0] if style_rows else {}).get("family") or ""),
        }


# ---------------------------------------------------------------------------
# Colormaps


class _ColormapRegistry:
    def __init__(self, builder: BundleBuilder):
        self.builder = builder
        self.entries: dict[str, dict[str, Any]] = {}
        self._by_blob: dict[str, str] = {}

    def add(self, option: dict[str, Any]) -> str:
        name = str(option.get("name") or "").strip()
        lut = option.get("lut_b64")
        if not name:
            return ""
        if not lut:
            return name if name in self.entries else ""
        try:
            data = np.frombuffer(base64.b64decode(lut), dtype=np.uint8)
        except Exception:
            return ""
        width = data.size // 4
        if width < 2:
            return ""
        blob = self.builder.array(data.reshape(width, 4), "u8", hint=f"lut-{name}", priority=CRITICAL, shuffle=False)
        existing = self._by_blob.get(blob)
        if name in self.entries and self.entries[name]["blob"] != blob:
            # Same name, different LUT (e.g. an opacity function baked into a
            # volume colormap): register under a disambiguated name.
            name = existing or f"{name}~{len(self.entries)}"
        if name not in self.entries:
            self.entries[name] = {
                "blob": blob,
                "width": int(width),
                "label": str(option.get("label") or name),
                "legendColor": css_hex(option.get("legend_color") or "#888888"),
            }
            self._by_blob.setdefault(blob, name)
        return name


# ---------------------------------------------------------------------------
# Volumes


def _decode_png_atlas(b64: str, nx: int, ny: int, nz: int, tiles_x: int, z_start: int, out: np.ndarray) -> int:
    from PIL import Image  # pillow ships with the scientific stack

    img = Image.open(io.BytesIO(base64.b64decode(b64)))
    arr = np.asarray(img.convert("L"), dtype=np.uint8)
    rows = arr.shape[0] // ny
    cols = arr.shape[1] // nx
    count = 0
    for k in range(rows * cols):
        z = z_start + k
        if z >= nz:
            break
        r, c = divmod(k, cols)
        out[z] = arr[r * ny:(r + 1) * ny, c * nx:(c + 1) * nx]
        count += 1
    return count


def _volume_bytes(layer: dict[str, Any]) -> tuple[np.ndarray, np.ndarray | None]:
    """Return uint8 voxels (z, y, x) and optional validity (z, y, x)."""

    shape = layer.get("shape") or {}
    nx, ny, nz = int(shape.get("x", 0)), int(shape.get("y", 0)), int(shape.get("z", 0))
    if min(nx, ny, nz) <= 0:
        raise CompileError(f"Volume {layer.get('key')!r} has no shape")

    def decode_payload(payload: dict[str, Any] | None, root_style: dict[str, Any]) -> np.ndarray | None:
        src = payload if isinstance(payload, dict) else root_style
        encoding = str(src.get("encoding") or src.get("data_encoding") or "uint8")
        b64 = src.get("data_b64")
        slabs = src.get("data_b64_slabs") or src.get("slabs")
        tiles = src.get("atlas_tiles") or src.get("data_atlas_tiles") or {}
        if encoding.startswith("png_atlas"):
            out = np.zeros((nz, ny, nx), dtype=np.uint8)
            if b64:
                _decode_png_atlas(b64, nx, ny, nz, int(tiles.get("x", 1) or 1), 0, out)
            elif isinstance(slabs, list):
                z = 0
                for slab in slabs:
                    data = slab.get("data_b64") if isinstance(slab, dict) else slab
                    z_start = int(slab.get("z_start", z)) if isinstance(slab, dict) else z
                    z += _decode_png_atlas(data, nx, ny, nz, int(tiles.get("x", 1) or 1), z_start, out)
            else:
                return None
            return out
        if b64:
            data = np.frombuffer(base64.b64decode(b64), dtype=np.uint8)
        elif isinstance(slabs, list) and slabs:
            parts = [base64.b64decode(s.get("data_b64") if isinstance(s, dict) else s) for s in slabs]
            data = np.frombuffer(b"".join(parts), dtype=np.uint8)
        else:
            return None
        if data.size != nx * ny * nz:
            raise CompileError(
                f"Volume {layer.get('key')!r}: {data.size} bytes for shape {nx}x{ny}x{nz}"
            )
        return data.reshape(nz, ny, nx)

    scalar = decode_payload(layer.get("scalar_payload"), layer)
    if scalar is None:
        raise CompileError(f"Volume {layer.get('key')!r} has no decodable voxel payload")
    validity = None
    vp = layer.get("validity_payload")
    if isinstance(vp, dict):
        validity = decode_payload(vp, {})
    elif layer.get("validity_b64"):
        validity = decode_payload(None, {"data_b64": layer.get("validity_b64"), "encoding": "uint8"})
    return scalar, validity


def occupancy_grid(vol: np.ndarray, max_blocks: int = 64) -> np.ndarray | None:
    """Conservative block maxima of a (z, y, x) uint8 volume.

    Matches the viewer's empty-space-skipping contract: at most
    ``max_blocks`` blocks per axis, voxel ``i`` belongs to block
    ``floor((i + 0.5) * g / n)``, and every voxel's value is spread to its
    one-voxel neighbourhood first so trilinear tails at block borders are
    always covered. Dilation and block reduction are separable max
    operations, so each axis is dilated and reduced in turn (memory shrinks
    as it goes).
    """

    if vol.size < (1 << 18) or not np.any(vol):
        return None
    out = vol
    for axis in range(3):
        n = out.shape[axis]
        g = min(max_blocks, n)
        dil = out.copy()
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis], hi[axis] = slice(0, -1), slice(1, None)
        np.maximum(dil[tuple(hi)], out[tuple(lo)], out=dil[tuple(hi)])
        np.maximum(dil[tuple(lo)], out[tuple(hi)], out=dil[tuple(lo)])
        block = np.minimum(g - 1, np.floor((np.arange(n) + 0.5) * g / n).astype(np.int64))
        starts = np.searchsorted(block, np.arange(g))
        out = np.maximum.reduceat(dil, starts, axis=axis)
    return np.ascontiguousarray(out)


def _compile_volume(builder: BundleBuilder, layer: dict[str, Any], colormaps: _ColormapRegistry,
                    presence: list[int] | None, co_rotation: float) -> dict[str, Any]:
    scalar, validity = _volume_bytes(layer)
    nz, ny, nx = scalar.shape
    bounds = layer.get("bounds") or {}
    options = layer.get("colormap_options") or []
    names = [colormaps.add(o) for o in options if isinstance(o, dict)]
    names = [n for n in names if n]
    defaults = dict(layer.get("default_controls") or {})
    data_range = layer.get("data_range") or [0.0, 1.0]
    out: dict[str, Any] = {
        "key": str(layer.get("key")),
        "stateKey": str(layer.get("state_key") or layer.get("key")),
        "name": str(layer.get("state_name") or layer.get("name") or layer.get("key")),
        "dims": [int(nx), int(ny), int(nz)],
        "bounds": [
            [float(bounds.get("x", [-0.5, 0.5])[0]), float(bounds.get("y", [-0.5, 0.5])[0]), float(bounds.get("z", [-0.5, 0.5])[0])],
            [float(bounds.get("x", [-0.5, 0.5])[1]), float(bounds.get("y", [-0.5, 0.5])[1]), float(bounds.get("z", [-0.5, 0.5])[1])],
        ],
        "data": {"blob": builder.array(scalar, "u8", hint=f"vol-{layer.get('key')}", priority=DEFERRED, shuffle=False)},
        "dataRange": [_num(data_range[0], 0.0), _num(data_range[1], 1.0)],
        "unit": str(layer.get("value_unit") or ""),
        "legendColor": css_hex(layer.get("legend_color") or "#888888"),
        "colormaps": names,
        "defaults": _clean_json(defaults),
        "visible": bool(layer.get("visible", True)),
        "timeMyr": None if layer.get("time_myr") is None else _num(layer.get("time_myr")),
        "onlyAtT0": bool(layer.get("only_at_t0", True)),
        "supportsShowAllTimes": bool(layer.get("supports_show_all_times", False)),
        "coRotate": bool(layer.get("co_rotate_with_frame", False)),
        "coRotationRate": co_rotation,
        "referenceTimeMyr": _num(layer.get("reference_time_myr"), 0.0),
        "interpolation": layer.get("interpolation", True) is not False,
        "variantGroup": layer.get("variant_group"),
        "exclusiveGroup": layer.get("exclusive_group"),
        "renderOrder": _num(layer.get("render_order"), 0.0),
        "opticalModel": "legacy-texture-space",
    }
    occ = occupancy_grid(scalar)
    if occ is not None:
        gz, gy, gx = occ.shape
        out["occupancy"] = {
            "blob": builder.array(occ, "u8", hint=f"vol-{layer.get('key')}-occ", priority=DEFERRED, shuffle=False),
            "dims": [int(gx), int(gy), int(gz)],
        }
    if validity is not None:
        out["validity"] = {
            "blob": builder.array(validity, "u8", hint=f"vol-{layer.get('key')}-valid", priority=DEFERRED, shuffle=False)
        }
    if presence is not None:
        out["presence"] = presence
    return out


# ---------------------------------------------------------------------------
# Images


def _image_media(data_url: str) -> tuple[bytes, str] | None:
    m = re.match(r"^data:([\w/+.-]+);base64,(.*)$", data_url or "", re.S)
    if not m:
        return None
    return base64.b64decode(m.group(2)), m.group(1)


# ---------------------------------------------------------------------------
# Sky members


def _radec_to_gal(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    e = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    g = _ICRS_TO_GAL @ e
    lon = np.degrees(np.arctan2(g[1], g[0])) % 360.0
    lat = np.degrees(np.arcsin(np.clip(g[2], -1, 1)))
    return lon, lat


def _compile_members(builder: BundleBuilder, sky_panel: dict[str, Any]) -> tuple[dict[str, Any] | None, dict[str, int]]:
    members = sky_panel.get("members_by_cluster") if isinstance(sky_panel, dict) else None
    if not isinstance(members, dict) or not members:
        return None, {}
    clusters: list[str] = []
    offsets = [0]
    signature_to_cluster: dict[tuple, int] = {}
    key_to_cluster: dict[str, int] = {}
    cols: dict[str, list[float]] = {k: [] for k in ("ra", "dec", "l", "b", "dist", "pmra", "pmdec")}
    for name, rows in members.items():
        if not isinstance(rows, list):
            continue
        stars = [r for r in rows if isinstance(r, dict) and r.get("is_cluster_member") is True]
        if not stars:
            continue
        sig = (len(stars), str(stars[0].get("source_id") or ""), round(_num(stars[0].get("ra"), 0.0), 6),
               str(stars[-1].get("source_id") or ""))
        if sig in signature_to_cluster:
            key_to_cluster[str(name)] = signature_to_cluster[sig]
            continue
        idx = len(clusters)
        clusters.append(str(name))
        signature_to_cluster[sig] = idx
        key_to_cluster[str(name)] = idx
        for s in stars:
            cols["ra"].append(_num(s.get("ra")))
            cols["dec"].append(_num(s.get("dec")))
            cols["l"].append(_num(s.get("l")))
            cols["b"].append(_num(s.get("b")))
            cols["dist"].append(_num(s.get("distance_pc")))
            cols["pmra"].append(_num(s.get("pmra_masyr")))
            cols["pmdec"].append(_num(s.get("pmdec_masyr")))
        offsets.append(len(cols["ra"]))
    if not clusters:
        return None, {}
    ra = np.array(cols["ra"])
    dec = np.array(cols["dec"])
    lon = np.array(cols["l"])
    lat = np.array(cols["b"])
    need = ~np.isfinite(lon) | ~np.isfinite(lat)
    if np.any(need & np.isfinite(ra) & np.isfinite(dec)):
        gl, gb = _radec_to_gal(ra[need], dec[need])
        lon[need], lat[need] = gl, gb
    dist = np.array(cols["dist"])
    pmra = np.nan_to_num(np.array(cols["pmra"]), nan=0.0)
    pmdec = np.nan_to_num(np.array(cols["pmdec"]), nan=0.0)
    # Fill missing distances with the cluster median so stars stay on the sky.
    off = np.array(offsets)
    for c in range(len(clusters)):
        seg = dist[off[c]:off[c + 1]]
        bad = ~np.isfinite(seg) | (seg <= 0)
        if np.any(bad):
            good = seg[~bad]
            seg[bad] = float(np.median(good)) if good.size else np.nan
    lr = np.radians(lon)
    br = np.radians(lat)
    unit = np.stack([np.cos(br) * np.cos(lr), np.cos(br) * np.sin(lr), np.sin(br)], axis=1)
    xyz = unit * dist[:, None]
    # Tangential velocity from proper motions, rotated ICRS → Galactic.
    rar = np.radians(np.nan_to_num(ra, nan=0.0))
    der = np.radians(np.nan_to_num(dec, nan=0.0))
    e_ra = np.stack([-np.sin(rar), np.cos(rar), np.zeros_like(rar)], axis=1)
    e_de = np.stack([-np.sin(der) * np.cos(rar), -np.sin(der) * np.sin(rar), np.cos(der)], axis=1)
    mu_icrs = pmra[:, None] * e_ra + pmdec[:, None] * e_de
    mu_gal = mu_icrs @ _ICRS_TO_GAL.T
    vel = mu_gal * (MAS_YR_KM_S_PER_PC * np.nan_to_num(dist, nan=0.0))[:, None] * KM_S_TO_PC_MYR
    good = np.all(np.isfinite(xyz), axis=1)
    xyz[~good] = 0.0
    vel[~np.isfinite(vel)] = 0.0
    spec = {
        "count": int(len(ra)),
        "clusters": clusters,
        "offsets": {"blob": builder.array(np.array(offsets, dtype=np.uint32), "u32", hint="members-offsets", priority=NORMAL)},
        "position": {"blob": builder.array(xyz.astype(np.float32), "f32", hint="members-xyz", priority=NORMAL)},
        "velocity": {"blob": builder.array(vel.astype(np.float32), "f32", hint="members-vel", priority=NORMAL)},
        "valid": {"blob": builder.array(good.astype(np.uint8), "u8", hint="members-valid", priority=NORMAL)},
    }
    return spec, key_to_cluster


# ---------------------------------------------------------------------------
# Top level


def _legend_items(spec: dict[str, Any]) -> list[dict[str, Any]]:
    items = (spec.get("legend") or {}).get("items") or []
    return [i for i in items if isinstance(i, dict)]


def _camera_from_spec(spec: dict[str, Any], center: list[float], max_span: float) -> dict[str, Any]:
    init = spec.get("initial_state") or {}
    cam = init.get("camera") if isinstance(init.get("camera"), dict) else None
    fov = _num(((init.get("global_controls") or {}).get("camera_fov")), 60.0)
    if cam and _xyz(cam.get("position")) and _xyz(cam.get("target")):
        return {"position": _xyz(cam.get("position")), "target": _xyz(cam.get("target")), "fov": fov}
    scene = (spec.get("layout") or {}).get("scene") or {}
    eye = (scene.get("camera") or {}).get("eye") or {"x": 0.0, "y": -0.8, "z": 1.5}
    scale = max(max_span, 1.0)
    pos = [
        center[0] + _num(eye.get("x"), 0.0) * scale,
        center[1] + _num(eye.get("y"), -0.8) * scale,
        center[2] + _num(eye.get("z"), 1.5) * scale,
    ]
    return {"position": pos, "target": list(center), "fov": fov}


def compile_scene_spec(spec: dict[str, Any], *, compress_level: int = 9) -> Bundle:
    """Convert a legacy Three.js scene spec into a viewer :class:`Bundle`."""

    frames = [f for f in (spec.get("frames") or []) if isinstance(f, dict)]
    if not frames:
        raise CompileError("Scene spec has no frames")
    order = sorted(range(len(frames)), key=lambda i: _num(frames[i].get("time"), float(i)))
    frames = [frames[i] for i in order]
    times = np.array([_num(f.get("time"), float(i)) for i, f in enumerate(frames)], dtype=np.float64)
    init_idx_raw = int(_num(spec.get("initial_frame_index"), -1))
    if 0 <= init_idx_raw < len(order):
        init_idx = order.index(init_idx_raw)
    else:
        zero = np.where(np.abs(times) < 1e-9)[0]
        init_idx = int(zero[0]) if zero.size else len(frames) - 1
    t0_candidates = np.where(np.abs(times) < 1e-9)[0]
    t0_index = int(t0_candidates[0]) if t0_candidates.size else init_idx

    builder = BundleBuilder(compress_level=compress_level)
    colormaps = _ColormapRegistry(builder)

    # ---- traces
    legend = _legend_items(spec)
    legend_by_key = {str(i.get("key")): i for i in legend}
    trace_keys: list[str] = []
    for item in legend:
        if item.get("kind", "trace") == "trace":
            trace_keys.append(str(item.get("key")))
    for f in frames:
        for t in f.get("traces") or []:
            k = str(t.get("key") or t.get("name"))
            if k not in trace_keys:
                trace_keys.append(k)
    by_frame: dict[str, list[dict | None]] = {k: [None] * len(frames) for k in trace_keys}
    for fi, f in enumerate(frames):
        for t in f.get("traces") or []:
            k = str(t.get("key") or t.get("name"))
            by_frame[k][fi] = t
    tc = _TraceCompiler(builder, times, t0_index, colormaps)
    traces = []
    for k in trace_keys:
        compiled = tc.compile(k, by_frame[k], legend_by_key.get(k))
        if compiled is not None:
            traces.append(compiled)

    # ---- decorations (images + volume presence)
    image_specs = {str(p.get("key")): p for p in (spec.get("image_planes") or []) if isinstance(p, dict)}
    image_centers: dict[str, np.ndarray] = {}
    image_opacity: dict[str, np.ndarray] = {}
    volume_presence: dict[str, list[int]] = {}
    extra_decorations: list[dict[str, Any]] = []
    for fi, f in enumerate(frames):
        for d in f.get("decorations") or []:
            if not isinstance(d, dict):
                continue
            kind = d.get("kind")
            if kind == "image_plane":
                key = str(d.get("key"))
                image_centers.setdefault(key, np.full((len(frames), 3), np.nan))
                image_opacity.setdefault(key, np.zeros(len(frames)))
                c = _xyz(d.get("center"))
                if c:
                    image_centers[key][fi] = c
                image_opacity[key][fi] = _num(d.get("opacity"), 1.0) * _num(d.get("opacity_scale"), 1.0)
            elif kind == "volume_layer":
                volume_presence.setdefault(str(d.get("key")), []).append(fi)
            elif fi == t0_index:
                extra_decorations.append(_clean_json(_strip_heavy(d)))
    images = []
    for key, plane in image_specs.items():
        media = _image_media(str(plane.get("image_data_url") or ""))
        if media is None:
            continue
        data, mime = media
        centers = image_centers.get(key)
        entry = {
            "key": key,
            "image": {"blob": builder.media(data, mime, hint=f"img-{key}", priority=NORMAL), "mime": mime},
            "widthPc": _num(plane.get("width_pc"), 40000.0),
            "heightPc": _num(plane.get("height_pc"), 40000.0),
            "hideBelowScalePc": _num(plane.get("hide_below_scale_bar_pc"), 0.0),
            "fadeStartScalePc": _num(plane.get("fade_start_scale_bar_pc"), 0.0),
            "renderOrder": _num(plane.get("render_order"), -20.0),
        }
        if centers is not None and not np.all(np.isnan(centers)):
            filled = centers.copy()
            last = None
            for i in range(len(filled)):
                if np.all(np.isfinite(filled[i])):
                    last = filled[i]
                elif last is not None:
                    filled[i] = last
            filled[np.isnan(filled)] = 0.0
            entry["center"] = {"blob": builder.array(filled, "f32", hint=f"img-{key}-center", priority=CRITICAL)}
            entry["opacity"] = {"blob": builder.array(image_opacity[key], "f32", hint=f"img-{key}-opacity", priority=CRITICAL)}
        else:
            entry["opacity"] = {"value": 1.0}
        images.append(entry)

    # ---- volumes
    vol_block = spec.get("volumes") or {}
    co_rotation = _num(vol_block.get("co_rotation_rate_rad_per_myr"), 0.0)
    volumes = []
    for layer in vol_block.get("layers") or []:
        if not isinstance(layer, dict):
            continue
        pres = volume_presence.get(str(layer.get("key")))
        presence = None if pres is None or len(pres) == len(frames) else sorted(pres)
        volumes.append(_compile_volume(builder, layer, colormaps, presence, co_rotation))

    # ---- sky
    sky_panel = spec.get("sky_panel") or {}
    sky_dome = spec.get("sky_dome") or {}
    members, member_keys = _compile_members(builder, sky_panel)
    if members:
        # Link member clusters to trace objects by name and aliases.
        for trace in traces:
            pts = trace.get("points")
            if not pts:
                continue
            meta = builder_blob_json(builder, pts["meta"]["blob"])
            names = meta.get("name") or []
            aliases = meta.get("aliases") or [""] * len(names)
            link = np.full(len(names), -1, dtype=np.int32)
            for j, nm in enumerate(names):
                cands = [nm, nm.replace("_", " "), nm.replace(" ", "_")]
                for alias in str(aliases[j] or "").split(","):
                    alias = alias.strip()
                    if alias:
                        cands.extend([alias, alias.replace("_", " "), alias.replace(" ", "_")])
                for c in cands:
                    if c in member_keys:
                        link[j] = member_keys[c]
                        break
            if np.any(link >= 0):
                pts["members"] = {"blob": builder.array(link, "i32", hint=f"{trace['key']}-members", priority=NORMAL)}
    init = spec.get("initial_state") or {}
    sky = None
    if sky_panel.get("enabled") or sky_dome.get("enabled"):
        sky = {
            "enabled": True,
            "survey": str(sky_panel.get("survey") or sky_dome.get("hips_survey") or "P/DSS2/color"),
            "frame": str(sky_panel.get("frame") or "galactic"),
            "layers": _clean_json(init.get("sky_layers") or sky_dome.get("sky_layers") or sky_dome.get("layers") or []),
            "layerGroups": _clean_json(sky_dome.get("layer_groups") or []),
            "activeLayerKey": init.get("active_sky_layer_key"),
            "activeGroupKey": init.get("active_sky_layer_group_key"),
            "memberPointSizeDenominator": _num(sky_panel.get("member_point_size_denominator"), 30.0),
            "memberMinScreenPx": _num(sky_panel.get("member_min_screen_size_px"), 0.0),
            "showMembers": bool(sky_panel.get("show_cluster_members_in_sky", bool(members))),
            "fullOpacityScalePc": _num(sky_dome.get("full_opacity_scale_bar_pc"), 120.0),
            "fadeOutScalePc": _num(sky_dome.get("fade_out_scale_bar_pc"), 360.0),
            "members": members,
        }

    # ---- world
    center = _xyz(spec.get("center")) or [0.0, 0.0, 0.0]
    max_span = _num(spec.get("max_span"), 1000.0)
    ref_span = _num(spec.get("point_size_reference_span_pc"), max_span)
    baseline = _num(spec.get("point_size_baseline_scale"), 4.0 / 3.0)
    ranges = spec.get("ranges") or {}
    sun_key = next((t["key"] for t in traces if t["name"].strip().lower() == "sun"), None)
    gc = None
    for layer in vol_block.get("layers") or []:
        c = ((layer or {}).get("default_controls") or {}).get("galactic_center")
        if isinstance(c, list) and len(c) == 3:
            gc = [float(v) for v in c]
            break
    world = {
        "center": center,
        "maxSpan": max_span,
        "ranges": _clean_json(ranges),
        "up": [0.0, 0.0, 1.0],
        "pointScale": (max(ref_span, 1.0) / 2600.0) * max(baseline, 0.01),
        "sunTrace": sun_key,
        "galacticCenter": gc or [8122.0, 0.0, 0.0],
        "showAxes": bool(spec.get("show_axes", False)),
    }

    animation = spec.get("animation") or {}
    states = spec.get("states") or {}
    manifest: dict[str, Any] = {
        "format": BUNDLE_FORMAT,
        "title": str(spec.get("title") or ""),
        "note": spec.get("note") or "",
        "theme": _clean_json(spec.get("theme") or {}),
        "time": {
            "values": [float(t) for t in times],
            "initialIndex": int(init_idx),
            "zeroIndex": int(t0_index),
            "unit": "Myr",
            "playbackIntervalMs": int(_num(spec.get("playback_interval_ms"), 240.0)),
            "enabled": bool((spec.get("timeline") or {}).get("enabled", len(frames) > 1)),
        },
        "world": world,
        "camera": _camera_from_spec(spec, center, max_span),
        "groups": {
            "order": [str(g) for g in spec.get("group_order") or []],
            "default": str(spec.get("default_group") or ""),
            "visibility": _clean_json(spec.get("group_visibility") or {}),
        },
        "animation": {
            "fadeTimeMyr": _num(animation.get("fade_in_time_default"), 0.0),
            "fadeInOut": bool(animation.get("fade_in_and_out_default", False)),
            "fadeByOpacity": bool(animation.get("fade_opacity_by_birth_time_default", False)),
        },
        "traces": traces,
        "images": images,
        "volumes": volumes,
        "colormaps": colormaps.entries,
        "sky": sky,
        "initialState": _remap_frame_fields(_clean_json(_strip_heavy(init)), order),
        "states": _remap_states_frames(_clean_json(states), order),
        "legacy": {
            "decorations": extra_decorations,
            "exportProfile": spec.get("export_profile"),
            "mobile": bool((spec.get("mobile") or {}).get("enabled")),
        },
    }
    if isinstance(spec.get("provenance"), dict):
        manifest["provenance"] = _clean_json(_strip_heavy(spec["provenance"]))
    return Bundle(manifest=manifest, blobs=builder.blobs)


def _remap_frame_value(value: Any, order: list[int]) -> Any:
    """Map a legacy (spec-order) frame index/value to the time-sorted order.

    The compiler sorts frames by time; legacy snapshots index frames in the
    order the spec stored them, which can be descending. Fractional values
    are mapped linearly between their neighbouring frames.
    """

    v = _num(value)
    n = len(order)
    if not math.isfinite(v) or n == 0:
        return value
    inverse = {old: new for new, old in enumerate(order)}
    v = min(max(v, 0.0), float(n - 1))
    lo = int(math.floor(v))
    hi = min(lo + 1, n - 1)
    t = v - lo
    a, b = inverse[lo], inverse[hi]
    out = a + (b - a) * t
    return int(round(out)) if isinstance(value, int) else float(out)


def _remap_frame_fields(snapshot: Any, order: list[int]) -> Any:
    if not isinstance(snapshot, dict) or order == sorted(order):
        return snapshot
    for field in ("current_frame_index", "current_frame_value"):
        if field in snapshot and snapshot[field] is not None:
            snapshot[field] = _remap_frame_value(snapshot[field], order)
    return snapshot


def _remap_states_frames(states: Any, order: list[int]) -> Any:
    if not isinstance(states, dict) or order == sorted(order):
        return states
    for item in states.get("items") or []:
        if isinstance(item, dict) and isinstance(item.get("snapshot"), dict):
            _remap_frame_fields(item["snapshot"], order)
    return states


def builder_blob_json(builder: BundleBuilder, blob_id: str) -> Any:
    return Bundle(manifest={}, blobs=builder.blobs).decode(blob_id)
