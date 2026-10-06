"""Annotations for Oviz figures: labels, curves and arrows, bubbles and shells.

Pass a list of dicts to ``make_plot(annotations=...)``. Positions are in pc,
in the figure's frame (heliocentric Galactic x, y, z at the present day)::

    annotations = [
        {"kind": "shell", "center": [0, 0, 0], "radius": 150, "label": "Local Bubble"},
        {"kind": "arrow", "points": [[100, -300, 0], [50, -150, 20]], "color": "#ffd27a"},
        {"kind": "curve", "points": [[0, 0, 0], [200, 100, 50], [400, 0, 0]], "dash": "dash"},
        {"kind": "text", "at": [300, 200, 0], "text": "Sco-Cen", "size": 16},
    ]

Kinds are ``text``, ``curve``, ``arrow`` (a curve with an arrowhead),
``bubble``, ``shell`` and ``wire`` (ellipsoids: ``radius`` or ``radii`` in pc,
``rot`` in degrees about x, y, z). A shape's ``label`` adds a label at its
centre. Items with the same ``group`` share one entry in the figure's key;
otherwise the viewer groups them as it does what readers draw. Every item
takes ``color``, ``opacity`` and ``present`` (show only at the present day).
Readers can edit all of them in the figure (press K).
"""

from __future__ import annotations

from typing import Any, Iterable

_SPHERES = {"bubble": "bubble", "shell": "shell", "wire": "wire", "sphere": "bubble"}
_STYLE_KEYS = ("color", "opacity", "present")


def _point(value: Any, what: str) -> dict[str, list[float]]:
    try:
        xyz = [float(v) for v in value]
    except TypeError as exc:
        raise ValueError(f"{what} must be [x, y, z] in pc") from exc
    if len(xyz) != 3:
        raise ValueError(f"{what} must be [x, y, z] in pc")
    return {"pos": xyz}


def normalize_annotations(items: Iterable[dict] | None) -> dict[str, list[dict]] | None:
    """Turn ``make_plot(annotations=...)`` items into the viewer's document.

    Returns ``{"groups": [...], "items": [...]}``, or None when there are none.
    Raises ``ValueError`` for an unknown kind or a missing position.
    """

    if not items:
        return None
    groups: dict[str, dict] = {}
    out: list[dict] = []

    def group_for(name: Any, fallback: str) -> str:
        key = str(name) if name is not None else fallback
        if key not in groups:
            gid = f"py-g{len(groups)}"
            groups[key] = {"id": gid, "name": str(name) if name is not None else "", "auto": name is None, "visible": True}
        return groups[key]["id"]

    for n, raw in enumerate(items):
        if not isinstance(raw, dict):
            raise ValueError(f"annotation {n} must be a dict")
        kind = str(raw.get("kind", "")).lower()
        style = {k: raw[k] for k in _STYLE_KEYS if k in raw}
        item_id = str(raw.get("id") or f"py-i{n}")
        if kind == "text":
            if "at" not in raw:
                raise ValueError(f"annotation {n} (text) needs 'at'")
            gid = group_for(raw.get("group"), f"item-{n}")
            item = {"kind": "text", "id": item_id, "group": gid, "at": _point(raw["at"], f"annotation {n} 'at'"), "text": str(raw.get("text", "Label")), **style}
            for k in ("size", "weight", "bg", "dx", "dy"):
                if k in raw:
                    item[k] = raw[k]
            out.append(item)
        elif kind in ("curve", "arrow"):
            pts = raw.get("points")
            if not pts or len(pts) < 2:
                raise ValueError(f"annotation {n} ({kind}) needs at least two 'points'")
            gid = group_for(raw.get("group"), f"item-{n}")
            item = {"kind": "curve", "id": item_id, "group": gid, "points": [_point(p, f"annotation {n} point") for p in pts],
                    "arrow": raw.get("arrow", "end" if kind == "arrow" else "none"), **style}
            for k in ("smooth", "closed", "width", "dash"):
                if k in raw:
                    item[k] = raw[k]
            out.append(item)
        elif kind in _SPHERES:
            if "center" not in raw:
                raise ValueError(f"annotation {n} ({kind}) needs 'center'")
            radii = raw.get("radii")
            if radii is None:
                if "radius" not in raw:
                    raise ValueError(f"annotation {n} ({kind}) needs 'radius' or 'radii'")
                radii = [raw["radius"]] * 3
            radii = [float(r) for r in radii]
            if len(radii) != 3 or min(radii) <= 0:
                raise ValueError(f"annotation {n} ({kind}) needs three positive radii")
            label = raw.get("label")
            gid = group_for(raw.get("group", label), f"item-{n}")
            center = _point(raw["center"], f"annotation {n} 'center'")
            out.append({"kind": "sphere", "id": item_id, "group": gid, "center": center, "radii": radii,
                        "rot": [float(r) for r in raw.get("rot", (0, 0, 0))], "style": raw.get("style", _SPHERES[kind]), **style})
            if label:
                out.append({"kind": "text", "id": f"{item_id}-label", "group": gid, "at": center, "text": str(label),
                            **{k: raw[k] for k in ("present",) if k in raw}})
        else:
            raise ValueError(f"annotation {n}: unknown kind {raw.get('kind')!r} (text, curve, arrow, bubble, shell, wire)")
    return {"groups": list(groups.values()), "items": out}
