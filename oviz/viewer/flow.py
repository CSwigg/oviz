"""Streamlines of steady 3D vector fields, drawn as animated flow layers.

A gas-velocity field sampled on a regular grid (for example a NIFTy
reconstruction on 50 pc cells) becomes a set of evenly spaced streamlines:

1. Candidate seeds are drawn inside the supported volume, weighted by
   ``seed_weight`` (usually gas density).
2. Every candidate is integrated forwards and backwards at once with a
   fourth-order Runge-Kutta step along the normalised velocity, through
   trilinear interpolation. A path stops where it leaves the grid or the
   support mask, where the speed drops below ``min_speed_kms``, or where it
   turns sharply (a stagnation point).
3. Paths are accepted greedily in seed-priority order. An occupancy grid
   truncates each path where it comes within about ``spacing_pc`` of a path
   already kept, so dense and diffuse regions both read clearly instead of
   piling up where seeds cluster.

Each vertex stores the integration time from the start of its line (Myr),
the local speed (km/s), a trust value and an end fade. The viewer moves
light pulses along each line in integration time, so on screen the pulses
travel at the gas's own relative speed. Nothing about the field is changed:
the geometry follows the velocities exactly and the extra columns are only
used for display.

:func:`flow_layer` wraps the result as a scene-spec entry. Pass it to
``make_plot(flows=[...])`` or append it to ``figure.scene_spec["flows"]``.
Only the Oviz viewer draws flows; the classic viewer ignores them.
"""

from __future__ import annotations

import base64
import math
from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

#: 1 km/s expressed in pc/Myr.
KMS_TO_PC_PER_MYR = 1.0227121650537077

#: Built-in colour ramps (position, hex) for flow speed.
FLOW_COLORMAPS: dict[str, list[tuple[float, str]]] = {
    # Deep blue → cyan → ice white: reads over warm (gist_heat / inferno) gas.
    "flow-ice": [(0.0, "#193f86"), (0.3, "#1f7cc6"), (0.6, "#2fc2e3"), (0.85, "#8ee6f6"), (1.0, "#d6f8ff")],
    # Teal → pale gold, for dark blue or grey backdrops.
    "flow-dawn": [(0.0, "#16505c"), (0.4, "#2fa39a"), (0.75, "#b9e28c"), (1.0, "#fff4c2")],
}


@dataclass
class Streamlines:
    """Evenly spaced streamlines through a steady vector field."""

    positions: np.ndarray  #: (V, 3) float32 vertices in pc.
    offsets: np.ndarray  #: (L + 1,) int64 start of each line in ``positions``.
    tau_myr: np.ndarray  #: (V,) integration time from the line's start.
    speed_kms: np.ndarray  #: (V,) speed at each vertex.
    trust: np.ndarray  #: (V,) trust in [0, 1] (1 when none was given).
    edge: np.ndarray  #: (V,) end fade in [0, 1].
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def count(self) -> int:
        return int(len(self.offsets) - 1)

    def lines(self) -> list[np.ndarray]:
        return [self.positions[a:b] for a, b in zip(self.offsets[:-1], self.offsets[1:])]


class _Grid:
    """Trilinear sampling on a regular (z, y, x) grid given by cell centres."""

    def __init__(self, axes: Sequence[np.ndarray]):
        x, y, z = (np.asarray(a, dtype=float) for a in axes)
        self.origin = np.array([x[0], y[0], z[0]])
        steps = []
        for a in (x, y, z):
            if a.size < 2:
                raise ValueError("each axis needs at least two cell centres")
            d = np.diff(a)
            if not np.allclose(d, d[0], rtol=1e-4, atol=1e-6) or d[0] <= 0:
                raise ValueError("axes must be increasing and uniformly spaced")
            steps.append(float(d[0]))
        self.step = np.array(steps)
        self.shape_xyz = np.array([x.size, y.size, z.size])

    def coords(self, p: np.ndarray) -> np.ndarray:
        return (p - self.origin) / self.step

    def inside(self, g: np.ndarray) -> np.ndarray:
        return np.all((g >= 0) & (g <= self.shape_xyz - 1), axis=-1)

    def nearest(self, g: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        i = np.clip(np.rint(g), 0, self.shape_xyz - 1).astype(np.int64)
        return i[..., 2], i[..., 1], i[..., 0]

    def trilinear(self, field: np.ndarray, g: np.ndarray) -> np.ndarray:
        """Sample ``field`` (nz, ny, nx[, c]) at grid coordinates ``g`` (N, 3)."""

        hi = self.shape_xyz - 1
        g = np.clip(g, 0, hi)
        i0 = np.minimum(np.floor(g).astype(np.int64), np.maximum(hi - 1, 0))
        t = g - i0
        ix, iy, iz = i0[:, 0], i0[:, 1], i0[:, 2]
        tx, ty, tz = t[:, 0], t[:, 1], t[:, 2]
        if field.ndim == 4:
            tx, ty, tz = tx[:, None], ty[:, None], tz[:, None]
        out = 0.0
        for dz, wz in ((0, 1 - tz), (1, tz)):
            for dy, wy in ((0, 1 - ty), (1, ty)):
                for dx, wx in ((0, 1 - tx), (1, tx)):
                    out = out + field[iz + dz, iy + dy, ix + dx] * (wz * wy * wx)
        return out


def _priority_order(weights: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Weighted random order without replacement (Efraimidis-Spirakis keys)."""

    w = np.maximum(np.asarray(weights, dtype=float), 1e-300)
    keys = np.log(rng.random(w.size)) / w
    return np.argsort(-keys, kind="stable")


def trace_streamlines(
    velocity: Sequence[np.ndarray],
    axes: Sequence[np.ndarray],
    *,
    support: np.ndarray | None = None,
    seed_weight: np.ndarray | None = None,
    trust: np.ndarray | None = None,
    spacing_pc: float | None = None,
    step_pc: float | None = None,
    max_length_pc: float = 900.0,
    min_length_pc: float | None = None,
    min_speed_kms: float = 0.3,
    max_turn_deg: float = 40.0,
    max_lines: int = 4000,
    n_candidates: int | None = None,
    edge_fade_pc: float | None = None,
    seed: int = 0,
) -> Streamlines:
    """Trace evenly spaced streamlines through a steady 3D velocity field.

    Parameters
    ----------
    velocity
        ``(vx, vy, vz)`` arrays of shape ``(nz, ny, nx)`` in km/s. Non-finite
        cells are treated as unsupported.
    axes
        ``(x, y, z)`` cell-centre coordinates in pc (uniformly spaced).
    support
        Boolean ``(nz, ny, nx)`` mask of cells where lines may run. A path
        stops on entering an unsupported cell (nearest-cell test).
    seed_weight
        Non-negative ``(nz, ny, nx)`` weights for where lines start; zero
        cells never seed. Defaults to uniform over ``support``.
    trust
        Optional ``(nz, ny, nx)`` values in [0, 1], interpolated to vertices
        so a viewer can fade poorly constrained flow.
    spacing_pc
        Target separation between neighbouring lines (default: one cell).
    step_pc
        Integration step (default: a quarter of the smaller of one cell and
        ``spacing_pc``).
    max_length_pc, min_length_pc
        Longest line kept (each direction integrates half of it) and the
        shortest line kept after spacing truncation (default 2.5 spacings).
    min_speed_kms
        Paths stop where the speed falls below this.
    max_turn_deg
        Paths stop where the direction turns more than this in one step.
    max_lines
        Stop once this many lines are kept.
    n_candidates
        Number of candidate seeds (default ``8 * max_lines``).
    edge_fade_pc
        Length over which the ``edge`` column ramps up at both line ends.
    seed
        Random seed; the result is deterministic for a given seed.
    """

    vx, vy, vz = (np.asarray(v, dtype=np.float64) for v in velocity)
    if not (vx.shape == vy.shape == vz.shape) or vx.ndim != 3:
        raise ValueError("velocity components must share one (nz, ny, nx) shape")
    grid = _Grid(axes)
    if tuple(grid.shape_xyz[::-1]) != vx.shape:
        raise ValueError(f"axes imply shape {tuple(grid.shape_xyz[::-1])}, velocity has {vx.shape}")
    finite = np.isfinite(vx) & np.isfinite(vy) & np.isfinite(vz)
    mask = finite if support is None else (np.asarray(support, dtype=bool) & finite)
    if not mask.any():
        raise ValueError("no supported cells to trace through")
    vel = np.stack([np.where(finite, vx, 0.0), np.where(finite, vy, 0.0), np.where(finite, vz, 0.0)], axis=-1)
    trust_field = None if trust is None else np.nan_to_num(np.clip(np.asarray(trust, dtype=np.float64), 0, 1))

    cell = float(grid.step.min())
    spacing = float(spacing_pc or cell)
    h = float(step_pc or 0.25 * min(cell, spacing))
    min_len = float(min_length_pc if min_length_pc is not None else 2.5 * spacing)
    fade = float(edge_fade_pc if edge_fade_pc is not None else 1.5 * spacing)
    half_steps = max(2, int(math.ceil(0.5 * max_length_pc / h)))
    cos_turn = math.cos(math.radians(max_turn_deg))
    rng = np.random.default_rng(seed)

    # ---- candidate seeds: weighted cells, jittered inside each cell
    w = np.ones(vx.shape) if seed_weight is None else np.nan_to_num(np.asarray(seed_weight, dtype=np.float64))
    w = np.where(mask, np.maximum(w, 0.0), 0.0)
    cells = np.flatnonzero(w > 0)
    if cells.size == 0:
        raise ValueError("seed_weight is zero everywhere inside the support")
    n_cand = int(n_candidates or 8 * max_lines)
    p = w.ravel()[cells]
    pick = rng.choice(cells.size, size=n_cand, replace=True, p=p / p.sum())
    iz, iy, ix = np.unravel_index(cells[pick], vx.shape)
    jitter = rng.random((n_cand, 3)) - 0.5
    seeds = grid.origin + (np.stack([ix, iy, iz], axis=1) + jitter) * grid.step
    order = _priority_order(p[pick], rng)
    seeds = seeds[order]

    def direction(pts: np.ndarray, alive: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        g = grid.coords(pts)
        ok = alive & grid.inside(g)
        zi, yi, xi = grid.nearest(g)
        ok &= mask[zi, yi, xi]
        v = grid.trilinear(vel, g)
        s = np.linalg.norm(v, axis=1)
        ok &= s >= min_speed_kms
        d = np.where(ok[:, None], v / np.maximum(s, 1e-12)[:, None], 0.0)
        return d, ok

    def integrate(sign: float) -> tuple[np.ndarray, np.ndarray]:
        n = len(seeds)
        path = np.full((half_steps + 1, n, 3), np.nan, dtype=np.float64)
        path[0] = seeds
        pts = seeds.copy()
        alive = np.ones(n, dtype=bool)
        prev = np.zeros((n, 3))
        length = np.zeros(n, dtype=np.int64)
        for k in range(1, half_steps + 1):
            k1, ok = direction(pts, alive)
            if k > 1:
                ok &= np.einsum("ij,ij->i", k1, prev) >= cos_turn
            k2, ok2 = direction(pts + sign * 0.5 * h * k1, ok)
            k3, ok3 = direction(pts + sign * 0.5 * h * k2, ok2)
            k4, ok4 = direction(pts + sign * h * k3, ok3)
            step = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            nxt = pts + sign * h * step
            g = grid.coords(nxt)
            ok4 &= grid.inside(g)
            zi, yi, xi = grid.nearest(g)
            ok4 &= mask[zi, yi, xi]
            alive = ok4
            if not alive.any():
                break
            pts = np.where(alive[:, None], nxt, pts)
            prev = np.where(alive[:, None], k1, prev)
            path[k, alive] = pts[alive]
            length[alive] = k
        return path, length

    fwd, nf = integrate(1.0)
    bwd, nb = integrate(-1.0)

    # ---- greedy acceptance with an occupancy grid of half-spacing cells
    occ_step = 0.5 * spacing
    lo = grid.origin - 0.5 * grid.step
    occ_shape = np.maximum(1, np.ceil((grid.shape_xyz * grid.step) / occ_step).astype(np.int64)) + 1
    blocked = np.zeros(tuple(occ_shape[::-1]), dtype=bool)  # (z, y, x), dilated by one cell
    offsets_3 = np.stack(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"), axis=-1).reshape(-1, 3)

    def occ_index(pts: np.ndarray) -> np.ndarray:
        return np.clip(np.floor((pts - lo) / occ_step).astype(np.int64), 0, occ_shape - 1)

    lines: list[np.ndarray] = []
    min_vertices = max(3, int(math.ceil(min_len / h)) + 1)
    for c in range(len(seeds)):
        if len(lines) >= max_lines:
            break
        back = bwd[1 : nb[c] + 1, c][::-1]
        line = np.concatenate([back, fwd[: nf[c] + 1, c]], axis=0)
        if len(line) < min_vertices:
            continue
        s0 = len(back)
        oi = occ_index(line)
        hit = blocked[oi[:, 2], oi[:, 1], oi[:, 0]]
        if hit[s0]:
            continue
        # Grow from the seed until the line meets another's neighbourhood.
        a = s0
        while a > 0 and not hit[a - 1]:
            a -= 1
        b = s0
        while b < len(line) - 1 and not hit[b + 1]:
            b += 1
        if b - a + 1 < min_vertices:
            continue
        kept = line[a : b + 1]
        lines.append(kept)
        cells_kept = np.unique(occ_index(kept), axis=0)
        near = (cells_kept[:, None, :] + offsets_3[None, :, :]).reshape(-1, 3)
        near = np.clip(near, 0, occ_shape - 1)
        blocked[near[:, 2], near[:, 1], near[:, 0]] = True

    if not lines:
        raise ValueError("no streamline survived the length and spacing tests")

    counts = np.array([len(l) for l in lines], dtype=np.int64)
    offsets = np.zeros(len(lines) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(counts)
    positions = np.concatenate(lines, axis=0)
    g = grid.coords(positions)
    speed = np.linalg.norm(grid.trilinear(vel, g), axis=1)
    trust_v = np.ones(len(positions)) if trust_field is None else np.clip(grid.trilinear(trust_field, g), 0, 1)
    tau = np.zeros(len(positions))
    edge = np.zeros(len(positions))
    for a, b in zip(offsets[:-1], offsets[1:]):
        seg = np.linalg.norm(np.diff(positions[a:b], axis=0), axis=1)
        mean_speed = 0.5 * (speed[a : b - 1] + speed[a + 1 : b])
        dt = seg / np.maximum(mean_speed * KMS_TO_PC_PER_MYR, 1e-9)
        tau[a + 1 : b] = np.cumsum(dt)
        s = np.concatenate([[0.0], np.cumsum(seg)])
        edge[a:b] = np.clip(np.minimum(s, s[-1] - s) / max(fade, 1e-9), 0, 1)

    lengths = (counts - 1) * h
    summary = {
        "method": "RK4 on normalised velocity, trilinear interpolation, greedy evenly spaced acceptance",
        "lines": int(len(lines)),
        "vertices": int(len(positions)),
        "candidates": int(n_cand),
        "spacing_pc": spacing,
        "step_pc": h,
        "max_length_pc": float(max_length_pc),
        "min_length_pc": min_len,
        "min_speed_kms": float(min_speed_kms),
        "max_turn_deg": float(max_turn_deg),
        "edge_fade_pc": fade,
        "seed": int(seed),
        "length_pc_p10_p50_p90": [float(v) for v in np.percentile(lengths, [10, 50, 90])],
        "speed_kms_p05_p50_p98": [float(v) for v in np.percentile(speed, [5, 50, 98])],
        "tau_myr_p50_p90_max_per_line": [
            float(v) for v in np.percentile(tau[offsets[1:] - 1], [50, 90, 100])
        ],
    }
    return Streamlines(
        positions=positions.astype(np.float32),
        offsets=offsets,
        tau_myr=tau.astype(np.float32),
        speed_kms=speed.astype(np.float32),
        trust=trust_v.astype(np.float32),
        edge=edge.astype(np.float32),
        summary=summary,
    )


def colormap_lut(colormap: str | Sequence[tuple[float, str]], width: int = 256) -> np.ndarray:
    """An RGBA uint8 lookup table from a built-in ramp, stop list or Matplotlib name."""

    from .colors import parse_color

    stops = FLOW_COLORMAPS.get(colormap) if isinstance(colormap, str) else list(colormap)
    x = np.linspace(0.0, 1.0, width)
    if stops is None:
        import matplotlib

        cmap = matplotlib.colormaps[str(colormap)]
        return np.round(np.asarray(cmap(x)) * 255).astype(np.uint8)
    pos = np.array([float(s[0]) for s in stops])
    rgb = np.array([parse_color(s[1])[:3] for s in stops], dtype=float)
    out = np.empty((width, 4))
    for c in range(3):
        out[:, c] = np.interp(x, pos, rgb[:, c])
    out[:, 3] = 1.0
    return np.round(out * 255).astype(np.uint8)


def flow_layer(
    streamlines: Streamlines,
    *,
    key: str = "flow",
    name: str = "Flow",
    colormap: str | Sequence[tuple[float, str]] = "flow-ice",
    speed_range: tuple[float, float] | None = None,
    values: np.ndarray | None = None,
    value_range: tuple[float, float] | None = None,
    value_label: str = "speed",
    color: str | None = None,
    width_px: float = 1.6,
    opacity: float = 0.9,
    rail_opacity: float = 0.16,
    trust_floor: float = 0.35,
    rate_myr_per_s: float = 1.2,
    period_myr: float = 14.0,
    tail_myr: float = 3.5,
    playing: bool = True,
    visible: bool = True,
    show_in_legend: bool = True,
    present_day_only: bool = False,
    hide_in_sky: bool = True,
    unit: str = "km/s",
    description: str = "",
) -> dict[str, Any]:
    """Scene-spec entry for an animated flow layer.

    Light pulses run along every line: a new pulse every ``period_myr`` of
    integration time, each with a tail ``tail_myr`` long, advancing
    ``rate_myr_per_s`` of integration time per second of playback. Because
    the clock is integration time, fast gas carries its pulses further per
    second and draws longer tails. ``rail_opacity`` keeps a faint trace of
    every line visible when paused. Colour follows speed through
    ``colormap`` over ``speed_range`` (default: the 2nd to 98th percentile),
    unless a fixed ``color`` is given. ``values`` (one per vertex, e.g. a
    signed line-of-sight velocity) colour the lines instead, over
    ``value_range``, and ``value_label`` names them. Vertices below full
    trust fade toward ``trust_floor``. ``present_day_only`` fades the layer
    out away from t = 0 (a map of today's gas); ``hide_in_sky`` leaves it
    out of Sky view. Readers can change the colours, range, width and
    pulses in the figure's layers panel.
    """

    colour_values = streamlines.speed_kms if values is None else np.asarray(values, dtype=np.float32).ravel()
    if colour_values.size != len(streamlines.positions):
        raise ValueError("values must give one number per streamline vertex")
    rng = value_range if values is not None else speed_range
    lo, hi = rng or tuple(float(v) for v in np.percentile(colour_values, [2, 98]))
    if not hi > lo:
        hi = lo + 1.0
    lut = colormap_lut(colormap)
    cmap_name = colormap if isinstance(colormap, str) else f"{key}-colormap"
    legend = color or "#{:02x}{:02x}{:02x}".format(*lut[int(0.7 * (len(lut) - 1)), :3])
    return {
        "key": str(key),
        "name": str(name),
        "description": str(description),
        "positions": np.ascontiguousarray(streamlines.positions, dtype=np.float32),
        "offsets": np.ascontiguousarray(streamlines.offsets, dtype=np.int64),
        "tau_myr": np.ascontiguousarray(streamlines.tau_myr, dtype=np.float32),
        "speed_kms": np.ascontiguousarray(streamlines.speed_kms, dtype=np.float32),
        "trust": np.ascontiguousarray(streamlines.trust, dtype=np.float32),
        "edge": np.ascontiguousarray(streamlines.edge, dtype=np.float32),
        "colormap": {
            "name": str(cmap_name),
            "label": str(cmap_name),
            "lut_b64": base64.b64encode(lut.tobytes()).decode("ascii"),
            "legend_color": legend,
        },
        "speed_range": [float(lo), float(hi)],
        "color_values": np.ascontiguousarray(colour_values, dtype=np.float32),
        "color_label": str(value_label),
        "unit": str(unit),
        "color": color,
        "legend_color": legend,
        "width_px": float(width_px),
        "opacity": float(opacity),
        "rail_opacity": float(rail_opacity),
        "trust_floor": float(trust_floor),
        "animation": {
            "rate_myr_per_s": float(rate_myr_per_s),
            "period_myr": float(period_myr),
            "tail_myr": float(tail_myr),
            "playing": bool(playing),
        },
        "visible": bool(visible),
        "show_in_legend": bool(show_in_legend),
        "present_day_only": bool(present_day_only),
        "hide_in_sky": bool(hide_in_sky),
        "summary": dict(streamlines.summary),
    }
