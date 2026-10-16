"""Kinetic tomography maps: reading, volumes coloured by velocity, and flow lines."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oviz import Scene3D, Trace, TraceCollection
from oviz.kt import KTMap, _nice, read_kt_map
from oviz.viz import _inline_volume_color_fields

h5py = pytest.importorskip("h5py")

NX, NY, NZ, CELL = 24, 20, 10, 50.0


def _axes():
    return tuple((np.arange(n) - 0.5 * (n - 1)) * CELL for n in (NX, NY, NZ))


def _fields():
    """(z, y, x) density with two clouds, a velocity growing with x, and a residual."""
    x, y, z = _axes()
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    rho = (np.exp(-((xx - 300) ** 2 + yy**2 + zz**2) / (2 * 120.0**2))
           + 0.6 * np.exp(-((xx + 350) ** 2 + (yy - 150) ** 2 + zz**2) / (2 * 100.0**2)))
    rho[rho < 0.02] = 0.0
    vel = 0.05 * xx                        # km/s: receding at +x, approaching at -x
    res = np.where(xx > 0, 6.0, -4.0)      # what is left after "rotation"
    vel[rho == 0] = 0.0                    # like the DIB maps: no velocity without gas
    return rho.astype(np.float32), vel.astype(np.float32), res.astype(np.float32)


@pytest.fixture(scope="module")
def kt_file(tmp_path_factory):
    path = tmp_path_factory.mktemp("kt") / "cubes_test_mean_std.h5"
    rho, vel, res = _fields()
    to_xyz = lambda a: np.ascontiguousarray(a.transpose(2, 1, 0))  # noqa: E731  (z, y, x) -> (x, y, z)
    with h5py.File(path, "w") as f:
        f["distances"] = np.array([CELL, CELL, CELL])
        f["mean/ew_cube"] = to_xyz(rho)
        f["mean/v_cube"] = to_xyz(vel)
        f["mean/v_residual_cube"] = to_xyz(res)
        f["mean/v0mean"] = np.float32(219.8)
        f["std/ew_cube"] = to_xyz(0.1 * rho)
        f["std/v_cube"] = to_xyz(np.full_like(vel, 2.0))
        f["std/v_residual_cube"] = to_xyz(np.full_like(res, 2.0))
        f.attrs["vreid"] = np.array([10.4, 15.1, 8.2, 236.0])
    return path


def test_reader_transposes_and_centres_the_grid(kt_file):
    kt = read_kt_map(kt_file, name="Test KT")
    rho, vel, res = _fields()
    assert kt.density.shape == (NZ, NY, NX)
    np.testing.assert_array_equal(kt.density, rho)
    np.testing.assert_array_equal(kt.velocity, vel)
    np.testing.assert_array_equal(kt.residual, res)
    x, y, z = kt.axes
    assert x[0] == pytest.approx(-0.5 * (NX - 1) * CELL) and x[-1] == pytest.approx(0.5 * (NX - 1) * CELL)
    assert kt.bounds["z"] == pytest.approx([-0.5 * NZ * CELL, 0.5 * NZ * CELL])
    assert kt.fields == ["velocity", "residual"]
    assert "219.8" in kt.residual_note
    assert kt.attrs["source"] == kt_file.name


def test_reader_without_a_residual(kt_file):
    kt = read_kt_map(kt_file, residual=None, velocity_std=None)
    assert kt.fields == ["velocity"] and kt.residual is None and kt.velocity_std is None
    with pytest.raises(ValueError):
        kt.flows(field="residual")
    with pytest.raises(KeyError):
        read_kt_map(kt_file, velocity="mean/nope")


def test_coarsening_weights_velocity_by_density():
    rho, vel, res = _fields()
    kt = KTMap(density=rho, velocity=vel, residual=res, axes=_axes(), velocity_std=np.full_like(vel, 2.0))
    c = kt.coarsened(2)
    assert c.density.shape == (NZ // 2, NY // 2, NX // 2)
    block = (slice(2, 4), slice(4, 6), slice(14, 16))
    w = rho[block].astype(float)
    assert c.density[1, 2, 7] == pytest.approx(w.mean(), rel=1e-6)
    if w.sum() > 0:
        assert c.velocity[1, 2, 7] == pytest.approx((w * vel[block]).sum() / w.sum(), rel=1e-5)
    np.testing.assert_allclose(c.axes[0], _axes()[0].reshape(-1, 2).mean(axis=1))
    # Blocks without gas have no measured velocity.
    assert np.all(np.isinf(c.velocity_std[c.density == 0]))
    assert kt.coarsened(1) is kt


def test_velocity_ranges_are_symmetric_and_rounded():
    rho, vel, res = _fields()
    kt = KTMap(density=rho, velocity=vel, residual=res, axes=_axes())
    lo, hi = kt.velocity_range()
    assert lo == -hi and 0 < hi <= _nice(float(np.abs(vel).max()))
    assert kt.velocity_range("residual") == (-6.0, 6.0)
    assert [_nice(v) for v in (0.7, 3.2, 12, 70.78, 84.26, 108.1)] == [0.7, 3.5, 12.0, 75.0, 85.0, 110.0]


def test_volume_carries_both_velocities(kt_file):
    kt = read_kt_map(kt_file)
    cfg = kt.volume(max_resolution=12)
    assert cfg["data"].shape == (NZ // 2, NY // 2, NX // 2)
    assert cfg["bounds"] == kt.coarsened(2).bounds
    assert [f["key"] for f in cfg["color_fields"]] == ["velocity", "residual"]
    assert cfg["color_fields"][1]["range"] == [-6.0, 6.0]
    assert cfg["color_by"] == "velocity" and "Residual" in cfg["display_note"]
    assert kt.volume(color_by="density")["color_by"] == "density"
    assert kt.volume(color_by="residual")["color_by"] == "residual"
    with pytest.raises(ValueError):
        kt.volume(color_by="metallicity")


def test_flows_run_along_sight_lines(kt_file):
    kt = read_kt_map(kt_file)
    for field in ("velocity", "residual"):
        flow = kt.flows(field=field, cell_pc=CELL, max_lines=60, density_quantile=0.2, max_length_pc=400)
        pos = flow["positions"].astype(float)
        off = flow["offsets"]
        assert len(off) > 2
        for a, b in zip(off[:-1], off[1:]):
            seg = pos[a + 1:b] - pos[a:b - 1]
            radial = pos[a:b - 1] / np.linalg.norm(pos[a:b - 1], axis=1, keepdims=True)
            cos = np.einsum("ij,ij->i", seg / np.linalg.norm(seg, axis=1, keepdims=True), radial)
            # Along r̂ (to the few degrees trilinear sampling of r̂ allows):
            # outwards where the gas recedes, inwards where it approaches.
            assert np.all(np.abs(cos) > 0.98)
            assert np.all(np.sign(cos) == np.sign(flow["color_values"][a:b - 1]))
        assert flow["key"] == ("kt-flow" if field == "velocity" else "kt-flow-residual")
        assert flow["present_day_only"] and flow["hide_in_sky"]
        assert flow["color"] == "#eef3ff" and flow["color_label"] == (kt.velocity_label if field == "velocity" else kt.residual_title)
        # Pulse timing follows the median speed: ~65 pc of travel per second.
        median = float(np.median(flow["speed_kms"]))
        assert flow["animation"]["rate_myr_per_s"] * max(median, 1.0) * 1.0227121650537077 == pytest.approx(65.0, rel=1e-6)


def test_kt_figure_compiles(kt_file):
    kt = read_kt_map(kt_file, name="Test KT")
    df = pd.DataFrame({"name": ["a"], "x": [0.0], "y": [0.0], "z": [0.0], "U": [0.0], "V": [0.0], "W": [0.0], "age_myr": [5.0]})
    fig = Scene3D(TraceCollection([Trace(df, data_name="c")]), figure_theme="dark").make_plot(
        time=np.arange(0, -3, -1.0),
        volumes=[kt.volume(max_resolution=12)],
        flows=[kt.flows(cell_pc=CELL, max_lines=40, density_quantile=0.2),
               kt.flows(field="residual", cell_pc=CELL, max_lines=40, density_quantile=0.2, visible=False)],
    )
    bundle = fig.bundle
    m = bundle.manifest
    vol = m["volumes"][0]
    assert [f["key"] for f in vol["colorFields"]] == ["velocity", "residual"] and vol["colorDefault"] == "velocity"
    nx, ny, nz = vol["dims"]
    for f in vol["colorFields"]:
        assert f["colormap"] in m["colormaps"] and m["colormaps"][f["colormap"]]["label"] == "RdBu_r"
        assert bundle.decode(f["data"]["blob"]).size == nx * ny * nz
    # The residual field's bytes put +6 km/s near the top of its ±6 range.
    res = bundle.decode(vol["colorFields"][1]["data"]["blob"]).reshape(nz, ny, nx)
    dens = bundle.decode(vol["data"]["blob"]).reshape(nz, ny, nx)
    assert res[dens > 0].max() >= 250 and res[dens > 0].min() <= 50
    traces = {t["key"]: t for t in m["traces"]}
    assert traces["kt-flow"]["presentDayOnly"] and traces["kt-flow"]["skyHidden"]
    assert m["initialState"]["legend_state"]["kt-flow-residual"] is False
    assert traces["kt-flow-residual"]["flow"]["speedRange"] == [-6.0, 6.0]


def test_color_fields_in_array_volumes():
    rho, vel, _ = _fields()
    q = (np.clip(rho / rho.max(), 0, 1) * 255).astype(np.uint8)
    options = [{"name": "magma"}]
    out = _inline_volume_color_fields({"color_data": vel, "color_range": (-20, 20)}, rho, rho.shape, q, options, None)
    (field,) = out["color_fields"]
    assert field["key"] == "field" and field["range"] == [-20.0, 20.0] and out["color_default"] == "field"
    assert any(o["name"] == field["colormap"] for o in options)
    assert _inline_volume_color_fields({}, rho, rho.shape, q, options, None) is None
    with pytest.raises(ValueError):
        _inline_volume_color_fields({"color_data": vel[:-1]}, rho, rho.shape, q, options, None)
    with pytest.raises(ValueError):
        _inline_volume_color_fields({"color_fields": [{"data": vel}], "color_by": "nope"}, rho, rho.shape, q, options, None)
    # Downsampling averages the colour by density: empty voxels do not dilute it.
    small = tuple(n // 2 for n in rho.shape)
    out = _inline_volume_color_fields({"color_data": vel, "color_range": (-20, 20)}, rho, small, q[::2, ::2, ::2], options, None)
    assert out["color_fields"][0]["data_b64"]
