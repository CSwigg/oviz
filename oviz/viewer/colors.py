"""CSS colour parsing shared by the bundle compiler."""

from __future__ import annotations

import re
from functools import lru_cache

_NAMED = {
    "black": (0, 0, 0), "white": (255, 255, 255), "red": (255, 0, 0),
    "green": (0, 128, 0), "lime": (0, 255, 0), "blue": (0, 0, 255),
    "yellow": (255, 255, 0), "cyan": (0, 255, 255), "aqua": (0, 255, 255),
    "magenta": (255, 0, 255), "fuchsia": (255, 0, 255), "gray": (128, 128, 128),
    "grey": (128, 128, 128), "orange": (255, 165, 0), "purple": (128, 0, 128),
    "pink": (255, 192, 203), "brown": (165, 42, 42), "gold": (255, 215, 0),
    "silver": (192, 192, 192), "navy": (0, 0, 128), "teal": (0, 128, 128),
    "maroon": (128, 0, 0), "olive": (128, 128, 0), "violet": (238, 130, 238),
    "indigo": (75, 0, 130), "coral": (255, 127, 80), "salmon": (250, 128, 114),
    "crimson": (220, 20, 60), "tomato": (255, 99, 71), "orchid": (218, 112, 214),
    "skyblue": (135, 206, 235), "steelblue": (70, 130, 180), "khaki": (240, 230, 140),
    "turquoise": (64, 224, 208), "lightgray": (211, 211, 211), "lightgrey": (211, 211, 211),
    "darkgray": (169, 169, 169), "darkgrey": (169, 169, 169), "dodgerblue": (30, 144, 255),
    "deepskyblue": (0, 191, 255), "limegreen": (50, 205, 50), "firebrick": (178, 34, 34),
    "darkorange": (255, 140, 0), "hotpink": (255, 105, 180), "slategray": (112, 128, 144),
    "transparent": (0, 0, 0),
}

_FUNC_RE = re.compile(r"^(rgba?|hsla?)\((.*)\)$")


def _hsl_to_rgb(h: float, s: float, l: float) -> tuple[float, float, float]:
    import colorsys

    r, g, b = colorsys.hls_to_rgb((h % 360) / 360.0, l, s)
    return r * 255, g * 255, b * 255


@lru_cache(maxsize=4096)
def parse_color(value: object, default: tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)):
    """Return an (r, g, b, a) tuple in 0..1 for any CSS-ish colour value."""

    if value is None:
        return default
    if isinstance(value, (tuple, list)) and len(value) in (3, 4):
        vals = [float(v) for v in value]
        scale = 255.0 if max(vals[:3]) > 1.0 else 1.0
        a = vals[3] if len(vals) == 4 else 1.0
        return (vals[0] / scale, vals[1] / scale, vals[2] / scale, a)
    text = str(value).strip().lower()
    if not text:
        return default
    if text.startswith("#"):
        h = text[1:]
        if len(h) in (3, 4):
            h = "".join(c * 2 for c in h)
        try:
            if len(h) == 6:
                return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, 1.0)
            if len(h) == 8:
                return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, int(h[6:8], 16) / 255)
        except ValueError:
            return default
        return default
    m = _FUNC_RE.match(text.replace(" ", ""))
    if m:
        kind, body = m.groups()
        parts = [p for p in re.split(r"[,/]", body) if p]
        try:
            nums = []
            for p in parts:
                if p.endswith("%"):
                    nums.append(float(p[:-1]) / 100.0)
                else:
                    nums.append(float(p.replace("deg", "")))
        except ValueError:
            return default
        if kind.startswith("rgb") and len(nums) >= 3:
            r, g, b = nums[:3]
            # Percent components were already divided to 0..1.
            if any(p.endswith("%") for p in parts[:3]):
                r, g, b = r * 255, g * 255, b * 255
            a = nums[3] if len(nums) > 3 else 1.0
            return (r / 255, g / 255, b / 255, max(0.0, min(1.0, a)))
        if kind.startswith("hsl") and len(nums) >= 3:
            r, g, b = _hsl_to_rgb(nums[0], nums[1], nums[2])
            a = nums[3] if len(nums) > 3 else 1.0
            return (r / 255, g / 255, b / 255, max(0.0, min(1.0, a)))
        return default
    if text in _NAMED:
        r, g, b = _NAMED[text]
        return (r / 255, g / 255, b / 255, 0.0 if text == "transparent" else 1.0)
    try:  # pragma: no cover - optional dependency
        from matplotlib.colors import to_rgba

        return tuple(float(c) for c in to_rgba(text))
    except Exception:
        return default


def css_hex(value: object, default: str = "#ffffff") -> str:
    """``#rrggbb`` for any CSS-ish colour value, or ``default`` when it does not parse."""
    r, g, b, _ = parse_color(value, (-1.0, -1.0, -1.0, 1.0))
    if r < 0:
        return default
    return "#{:02x}{:02x}{:02x}".format(
        int(round(max(0, min(1, r)) * 255)),
        int(round(max(0, min(1, g)) * 255)),
        int(round(max(0, min(1, b)) * 255)),
    )
