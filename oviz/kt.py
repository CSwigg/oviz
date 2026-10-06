"""Kinetic tomography (KT) maps: 3D gas density with its line-of-sight velocity.

A KT map gives, on a heliocentric Cartesian grid, the density of some
tracer of the interstellar medium (dust, a diffuse interstellar band, H I,
CO...) and the mean line-of-sight velocity of that gas, often also the
residual velocity left after subtracting a model of Galactic rotation. Oviz
draws it two ways, together or apart:

* a **volume** whose opacity follows density and whose colour shows a
  velocity (``RdBu_r``: receding gas red, approaching gas blue) or the
  density; readers switch between the velocity, the residual and the
  density in the layers panel;
* animated **flow lines**. A KT map measures only the line-of-sight part of
  the motion, so the lines run along sight lines through the gas: pulses
  travel away from the Sun where the gas recedes and towards it where it
  approaches, at the gas's own speed. Trace them from the full velocity or
  from the residual (the gas's motion relative to Galactic rotation).

::

    from oviz.kt import read_kt_map

    kt = read_kt_map("cubes_full_mean_std.h5")      # DIB density, v_los, residual
    figure = scene.make_plot(
        time=time_myr,
        volumes=[kt.volume()],
        flows=[kt.flows(), kt.flows(field="residual", visible=False)],
    )

:func:`read_kt_map` reads HDF5 files that store ``(x, y, z)``-ordered cubes
with their voxel sizes in ``distances`` on a Sun-centred grid (the layout of
the DIB KT maps); other layouts go straight into :class:`KTMap`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np

__all__ = ["KTMap", "read_kt_map"]

_FIELDS = ("velocity", "residual")


@dataclass
class KTMap:
    """A KT map on a regular heliocentric grid.

    ``density``, ``velocity`` and the optional ``residual`` are
    ``(nz, ny, nx)`` arrays; ``axes`` holds the ``(x, y, z)`` cell centres
    in pc (Galactic Cartesian: x towards the Galactic centre, y towards
    l = 90°, z to the north Galactic pole). Velocities are line-of-sight
    km/s, positive receding; ``residual`` is what remains after subtracting
    a model of Galactic rotation (``residual_note`` says which, in the
    figure). The optional ``*_std`` arrays (same shape) are 1σ
    uncertainties.
    """

    density: np.ndarray
    velocity: np.ndarray
    axes: tuple[np.ndarray, np.ndarray, np.ndarray]
    density_std: np.ndarray | None = None
    velocity_std: np.ndarray | None = None
    residual: np.ndarray | None = None
    residual_std: np.ndarray | None = None
    name: str = "KT map"
    density_label: str = "Density"
    density_unit: str = ""
    velocity_label: str = "vₗₒₛ"
    residual_label: str = "Residual"
    residual_title: str = "vₗₒₛ − rotation"
    residual_note: str = ""
    velocity_unit: str = "km/s"
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.density = np.asarray(self.density, dtype=np.float32)
        if self.density.ndim != 3:
            raise ValueError("density must be a (nz, ny, nx) array")
        for name in ("velocity", "density_std", "velocity_std", "residual", "residual_std"):
            arr = getattr(self, name)
            if arr is None:
                continue
            arr = np.asarray(arr, dtype=np.float32)
            if arr.shape != self.density.shape:
                raise ValueError(f"{name} must have the density's shape {self.density.shape}")
            setattr(self, name, arr)
        self.axes = tuple(np.asarray(a, dtype=np.float64) for a in self.axes)
        if len(self.axes) != 3 or tuple(a.size for a in self.axes[::-1]) != self.density.shape:
            raise ValueError(f"axes (x, y, z) must have {self.density.shape[::-1]} cell centres")
        for a in self.axes:
            d = np.diff(a)
            if a.size < 2 or not np.all(d > 0) or not np.allclose(d, d[0], rtol=1e-4):
                raise ValueError("each axis must be increasing and uniformly spaced")

    # ------------------------------------------------------------------ grid

    @property
    def cell_pc(self) -> np.ndarray:
        """Cell size along x, y, z in pc."""
        return np.array([a[1] - a[0] for a in self.axes])

    @property
    def bounds(self) -> dict[str, list[float]]:
        """Outer edges of the grid (cell centres ± half a cell), in pc."""
        return {k: [float(a[0] - 0.5 * (a[1] - a[0])), float(a[-1] + 0.5 * (a[1] - a[0]))]
                for k, a in zip("xyz", self.axes)}

    @property
    def fields(self) -> list[str]:
        """The velocity fields this map has: ``"velocity"`` and, if given, ``"residual"``."""
        return [f for f in _FIELDS if getattr(self, f) is not None]

    def _field(self, name: str) -> tuple[np.ndarray, np.ndarray | None, str, str]:
        """(values, σ, short label, title) of a velocity field."""
        key = str(name).lower()
        if key not in _FIELDS:
            raise ValueError(f"field must be one of {_FIELDS}")
        values = getattr(self, key)
        if values is None:
            raise ValueError(f"this KT map has no {key} field")
        if key == "velocity":
            return values, self.velocity_std, self.velocity_label, self.velocity_label
        return values, self.residual_std, self.residual_label, self.residual_title

    def coarsened(self, factor: int) -> "KTMap":
        """The map on cells ``factor`` times larger along every axis.

        Density is the block mean and each velocity the density-weighted
        mean (the velocity of the gas in the block). The σ of a velocity is
        the density-weighted mean σ of its cells, so it never claims the
        precision of independent measurements. A remainder of cells that
        does not fill a block is dropped at the high end of each axis.
        """

        factor = int(factor)
        if factor <= 1:
            return self
        nz, ny, nx = (n // factor * factor for n in self.density.shape)
        if min(nz, ny, nx) < 2 * factor:
            raise ValueError("the map is too small to coarsen that much")

        def blocks(a: np.ndarray) -> np.ndarray:
            a = a[:nz, :ny, :nx]
            return a.reshape(nz // factor, factor, ny // factor, factor, nx // factor, factor).sum(axis=(1, 3, 5))

        rho = np.nan_to_num(np.clip(self.density, 0, None)).astype(np.float64)
        mass = blocks(rho)

        def weighted(a: np.ndarray | None, *, empty: float | None = None) -> np.ndarray | None:
            """Density-weighted block means; blocks without gas get ``empty``, or their plain mean."""
            if a is None:
                return None
            ok = np.isfinite(a)
            w = np.where(ok, rho, 0.0)
            num, den = blocks(w * np.where(ok, a, 0.0)), blocks(w)
            if empty is None:
                fill = blocks(np.where(ok, a, 0.0)) / np.maximum(blocks(ok.astype(np.float64)), 1.0)
            else:
                fill = np.full(den.shape, empty)
            return np.where(den > 0, num / np.where(den > 0, den, 1.0), fill).astype(np.float32)

        axes = tuple(a[:n].reshape(-1, factor).mean(axis=1) for a, n in zip(self.axes, (nx, ny, nz)))
        std = None if self.density_std is None else (blocks(np.nan_to_num(self.density_std)) / factor**3).astype(np.float32)
        return replace(
            self,
            density=(mass / factor**3).astype(np.float32),
            velocity=weighted(self.velocity),
            residual=weighted(self.residual),
            density_std=std,
            velocity_std=weighted(self.velocity_std, empty=np.inf),
            residual_std=weighted(self.residual_std, empty=np.inf),
            axes=axes,
        )

    def _gas(self) -> np.ndarray:
        d = self.density[np.isfinite(self.density) & (self.density > 0)]
        return d if d.size else np.ones(1, dtype=np.float32)

    def velocity_range(self, field: str = "velocity", quantile: float = 0.9) -> tuple[float, float]:
        """A symmetric range of ``field`` covering ``quantile`` of the gas (density-weighted)."""
        values = self._field(field)[0]
        w = np.nan_to_num(np.clip(self.density, 0, None)).ravel()
        v = np.abs(np.nan_to_num(values, posinf=0.0, neginf=0.0)).ravel()
        top = float(v.max()) if v.size else 0.0
        if not (w.sum() > 0 and top > 0):
            return (-1.0, 1.0)
        # A weighted quantile from a fine histogram (a sort of every cell is slow).
        hist, edges = np.histogram(v, bins=8192, range=(0.0, top), weights=w)
        cum = np.cumsum(hist)
        m = float(edges[1:][min(np.searchsorted(cum, quantile * cum[-1]), hist.size - 1)])
        m = _nice(m) if m > 0 else 1.0
        return (-m, m)

    # ---------------------------------------------------------------- volume

    def volume(
        self,
        *,
        name: str | None = None,
        key: str | None = None,
        color_by: str = "velocity",
        colormap: str = "magma",
        velocity_colormap: str = "RdBu_r",
        velocity_range: Sequence[float] | None = None,
        residual_range: Sequence[float] | None = None,
        density_range: Sequence[float] | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        max_resolution: int = 256,
        opacity: float = 1.0,
        alpha_coef: float = 10.0,
        samples: int = 256,
        stretch: str = "asinh",
        visible: bool = True,
        **extra: Any,
    ) -> dict[str, Any]:
        """A ``make_plot(volumes=[...])`` entry for this map.

        Opacity follows the density. ``color_by`` picks what colours it when
        the figure opens: ``"velocity"``, ``"residual"`` (each through
        ``velocity_colormap`` over its range, by default symmetric about
        zero and covering 90% of the gas) or ``"density"`` (through
        ``colormap``); readers can switch. The map is block-averaged to at
        most ``max_resolution`` cells along its longest axis (velocities
        density-weighted). Density is stored as 8-bit levels over
        ``density_range`` (default 0 to its 99.5th percentile, so the
        densest cores saturate) and opens with the window ``vmin``–``vmax``
        (default: the 85th percentile of the gas to the top of the range).
        Other keys (``only_at_t0``, ``lighting_mode``...) pass through.
        """

        mode = str(color_by).lower()
        if mode not in ("density", "value") and mode not in self.fields:
            raise ValueError(f"color_by must be 'density' or one of {self.fields}")
        factor = max(1, math.ceil(max(self.density.shape) / int(max_resolution)))
        kt = self.coarsened(factor)
        gas = kt._gas()
        lo, hi = density_range if density_range is not None else (0.0, float(np.percentile(gas, 99.5)))
        ranges = {"velocity": velocity_range, "residual": residual_range}
        fields = []
        for f in kt.fields:
            values, _, label, title = kt._field(f)
            rng = ranges[f] if ranges[f] is not None else self.velocity_range(f)
            fields.append({"key": f, "data": values, "label": label, "title": title, "unit": self.velocity_unit,
                           "range": [float(rng[0]), float(rng[1])], "colormap": velocity_colormap})
        cfg = {
            "name": name or self.name,
            "data": kt.density,
            "bounds": kt.bounds,
            "data_range": [float(lo), float(hi)],
            "vmin": float(vmin) if vmin is not None else float(np.percentile(gas, 85)),
            "vmax": float(vmax) if vmax is not None else float(hi),
            "colormap": colormap,
            "unit_label": self.density_unit,
            "opacity": float(opacity),
            "alpha_coef": float(alpha_coef),
            "samples": int(samples),
            "stretch": stretch,
            "visible": bool(visible),
            "color_fields": fields,
            "color_by": "density" if mode in ("density", "value") else mode,
        }
        if self.residual is not None and self.residual_note:
            cfg["display_note"] = f"{self.residual_label}: {self.residual_note}."
        
        if key:
            cfg["key"] = key
        cfg.update(extra)
        return cfg

    # ----------------------------------------------------------------- flows

    def flows(
        self,
        *,
        field: str = "velocity",
        name: str | None = None,
        key: str | None = None,
        cell_pc: float = 60.0,
        spacing_pc: float | None = None,
        max_lines: int = 1500,
        max_length_pc: float = 800.0,
        min_speed_kms: float = 2.0,
        density_quantile: float = 0.6,
        min_snr: float | None = 1.0,
        velocity_range: Sequence[float] | None = None,
        colormap: str = "RdBu_r",
        color: str | None = "#eef3ff",
        width_px: float = 1.1,
        opacity: float = 0.75,
        rail_opacity: float = 0.06,
        rate_myr_per_s: float | None = None,
        period_myr: float | None = None,
        tail_myr: float | None = None,
        present_day_only: bool = True,
        visible: bool = True,
        seed: int = 0,
    ) -> dict[str, Any]:
        """A ``make_plot(flows=[...])`` entry: animated lines along the sight lines.

        ``field`` is ``"velocity"`` (the measured line-of-sight velocity) or
        ``"residual"`` (what is left after subtracting Galactic rotation).
        The map is averaged onto cells about ``cell_pc`` wide (velocity
        density-weighted). Lines run where the density is above its
        ``density_quantile`` (of cells with gas), start preferentially in
        dense gas, keep about ``spacing_pc`` apart (default 1.5 cells) and
        stop where the velocity changes sign or drops below
        ``min_speed_kms``. With uncertainties, gas whose velocity is
        measured below ``min_snr`` σ carries no lines, and lines fade where
        it is poorly constrained. The lines are pale ``color`` (as wind
        particles over a coloured map) and also carry the signed velocity,
        so readers can colour them by it (``colormap`` over
        ``velocity_range``); ``color=None`` opens them that way. Pulses
        move ``rate_myr_per_s`` of gas travel time per second, one every
        ``period_myr``, trailing ``tail_myr``. By default these follow the
        lines' median speed, so a typical pulse crosses about 65 pc of the
        figure per second, 260 pc behind the one before, trailing 70 pc
        (the residual's slower gas then reads as clearly as the full
        velocity's); within a layer, faster gas still moves faster.
        Readers can change all of these in the figure.
        """

        from .viewer.flow import KMS_TO_PC_PER_MYR, _Grid, flow_layer, trace_streamlines

        title = self._field(field)[3]
        factor = max(1, int(round(float(cell_pc) / float(self.cell_pc.min()))))
        kt = self.coarsened(factor)
        velocity, vstd, _, _ = kt._field(field)
        density = kt.density
        x, y, z = kt.axes
        r = np.sqrt(x[None, None, :] ** 2 + y[None, :, None] ** 2 + z[:, None, None] ** 2)
        r = np.where(r > 0, r, np.inf)
        vx = velocity * (x[None, None, :] / r)
        vy = velocity * (y[None, :, None] / r)
        vz = velocity * (z[:, None, None] / r)
        gas = density[np.isfinite(density) & (density > 0)]
        if gas.size == 0:
            raise ValueError("the map has no gas to draw flow lines through")
        support = density > float(np.quantile(gas, density_quantile))
        trust = None
        if vstd is not None:
            snr = np.abs(velocity) / np.maximum(vstd, 1e-6)
            if min_snr is not None:
                support &= snr >= float(min_snr)
            # Fade lines through gas whose velocity is under ~3σ.
            trust = np.clip((snr - 1.0) / 2.0, 0.0, 1.0)
        cell = float(min(a[1] - a[0] for a in kt.axes))
        weight = np.where(support, np.clip(density, 0, float(np.percentile(gas, 99))), 0.0)
        lines = trace_streamlines(
            (vx, vy, vz), kt.axes,
            support=support, seed_weight=weight, trust=trust,
            spacing_pc=float(spacing_pc or 1.5 * cell),
            max_length_pc=float(max_length_pc), min_speed_kms=float(min_speed_kms),
            max_turn_deg=60.0, max_lines=int(max_lines), seed=int(seed),
        )
        grid = _Grid(kt.axes)
        signed = grid.trilinear(velocity.astype(np.float64), grid.coords(lines.positions.astype(np.float64)))
        rng = tuple(velocity_range) if velocity_range is not None else self.velocity_range(field)
        # Pulse timing in Myr from lengths on screen, at the lines' median speed.
        pc_per_myr = max(float(np.median(lines.speed_kms)), 1.0) * KMS_TO_PC_PER_MYR
        rate = float(rate_myr_per_s) if rate_myr_per_s is not None else 65.0 / pc_per_myr
        period = float(period_myr) if period_myr is not None else 260.0 / pc_per_myr
        tail = float(tail_myr) if tail_myr is not None else 70.0 / pc_per_myr
        residual = str(field).lower() == "residual"
        what = "relative to Galactic rotation" if residual else "along the line of sight"
        entry = flow_layer(
            lines, key=key or ("kt-flow-residual" if residual else "kt-flow"),
            name=name or f"{self.name} {'residual flow' if residual else 'flow'}", colormap=colormap,
            values=signed.astype(np.float32), value_range=(float(rng[0]), float(rng[1])),
            value_label=title, unit=self.velocity_unit, color=color,
            width_px=width_px, opacity=opacity, rail_opacity=rail_opacity,
            rate_myr_per_s=rate, period_myr=period, tail_myr=tail,
            visible=visible, present_day_only=present_day_only,
            description=(f"Gas motion {what} in {self.name}: pulses move away from the Sun where the gas "
                         "recedes and towards it where it approaches."
                         + (f" {self.residual_label}: {self.residual_note}." if residual and self.residual_note else "")),
        )
        entry["summary"] = {**entry.get("summary", {}), "field": f"{title} along r̂ (from the Sun)",
                            "cell_pc": cell, "density_quantile": float(density_quantile),
                            **({"min_snr": float(min_snr)} if vstd is not None and min_snr is not None else {})}
        return entry


def _nice(x: float) -> float:
    """``x`` rounded up to two significant figures, the second a multiple of 5 above 20 (84 → 85, 108 → 110)."""
    step = 10.0 ** (np.floor(np.log10(x)) - 1)
    m = x / step * (1 - 1e-12)
    m = np.ceil(m / 5) * 5 if m >= 20 else np.ceil(m)
    return float(round(m * step, 12))


def read_kt_map(
    path: str | Path,
    *,
    density: str = "mean/ew_cube",
    velocity: str = "mean/v_cube",
    residual: str | None = "mean/v_residual_cube",
    density_std: str | None = "std/ew_cube",
    velocity_std: str | None = "std/v_cube",
    residual_std: str | None = "std/v_residual_cube",
    name: str | None = None,
    density_label: str = "DIB density",
    density_unit: str = "Å/pc",
    velocity_label: str = "vₗₒₛ",
    residual_label: str = "Residual",
    residual_title: str = "vₗₒₛ − rotation",
    residual_note: str | None = None,
    cell_pc: Sequence[float] | None = None,
) -> KTMap:
    """Read a KT map from HDF5.

    The cubes are ``(x, y, z)``-ordered (as written by NumPy from an
    ``[ix, iy, iz]`` grid) on a grid centred on the Sun, with the voxel
    sizes in pc in the ``distances`` dataset (or ``cell_pc``). Defaults
    follow the DIB KT maps (``cubes_*_mean_std.h5``): equivalent-width
    density ``mean/ew_cube``, line-of-sight velocity ``mean/v_cube`` and
    the residual ``mean/v_residual_cube``, which subtracts a flat rotation
    curve with the map's own fitted circular speed (``mean/v0mean``; use
    ``mean/v_residual_reid_cube`` for 236 km/s), with uncertainties under
    ``std/``. Missing residual and uncertainty datasets are skipped.
    """

    import h5py

    path = Path(path)
    with h5py.File(path, "r") as f:
        def cube(key: str | None, required: bool = True) -> np.ndarray | None:
            if key is None or (not required and key not in f):
                return None
            if key not in f:
                raise KeyError(f"{path.name} has no dataset {key!r}")
            arr = np.asarray(f[key][...], dtype=np.float32)
            if arr.ndim != 3:
                raise ValueError(f"{key!r} is not a 3D cube")
            return np.ascontiguousarray(arr.transpose(2, 1, 0))  # (x, y, z) -> (z, y, x)

        rho = cube(density)
        vel = cube(velocity)
        res = cube(residual, required=False)
        rho_std = cube(density_std, required=False)
        vel_std = cube(velocity_std, required=False)
        res_std = cube(residual_std, required=False) if res is not None else None
        if cell_pc is None:
            if "distances" not in f:
                raise KeyError(f"{path.name} has no 'distances' (voxel sizes); pass cell_pc=(dx, dy, dz)")
            cell_pc = np.asarray(f["distances"][...], dtype=float).ravel()
        v0 = None
        if residual and residual.endswith("v_residual_cube") and "mean/v0mean" in f:
            v0 = float(np.asarray(f["mean/v0mean"][()]))
        elif residual and "reid" in residual and "vreid" in f.attrs:
            v0 = float(np.asarray(f.attrs["vreid"]).ravel()[-1])
        attrs = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in f.attrs.items()}
    if residual_note is None:
        residual_note = f"vₗₒₛ minus a flat rotation curve at {v0:.1f} km/s" if v0 else "vₗₒₛ minus Galactic rotation"
    dx, dy, dz = (float(c) for c in cell_pc)
    nz, ny, nx = rho.shape
    axes = tuple((np.arange(n) - 0.5 * (n - 1)) * d for n, d in ((nx, dx), (ny, dy), (nz, dz)))
    return KTMap(
        density=rho, velocity=vel, axes=axes, density_std=rho_std, velocity_std=vel_std,
        residual=res, residual_std=res_std,
        name=name or path.stem, density_label=density_label, density_unit=density_unit,
        velocity_label=velocity_label, residual_label=residual_label, residual_title=residual_title,
        residual_note=residual_note,
        attrs={**attrs, "source": path.name, "velocity_dataset": velocity, "residual_dataset": residual if res is not None else None},
    )
