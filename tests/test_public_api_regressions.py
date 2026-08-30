import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from oviz import Animate3D, Trace, TraceCollection
from oviz import orbit_maker
from oviz.spiral_arms import SpiralArmTrace


def _trace():
    return Trace(
        pd.DataFrame(
            {
                "x": [10.0],
                "y": [20.0],
                "z": [30.0],
                "U": [1.0],
                "V": [2.0],
                "W": [3.0],
                "name": ["A"],
                "age_myr": [8.0],
            }
        ),
        data_name="Clusters",
    )


def _coordinate_groups(time):
    time = np.asarray(time, dtype=float)
    shape = (1, len(time))
    centered = tuple(np.full(shape, value) for value in (1.0, 2.0, 3.0))
    heliocentric = tuple(np.full(shape, value) for value in (10.0, 20.0, 30.0))
    galactocentric = tuple(np.full(shape, value) for value in (100.0, 200.0, 300.0))
    cylindrical = tuple(np.full(shape, value) for value in (400.0, 0.5, 500.0))
    return centered, heliocentric, galactocentric, cylindrical


@pytest.mark.parametrize(
    "time, message",
    [
        (None, "non-empty one-dimensional"),
        ([[0.0, -1.0]], "one-dimensional"),
        ([0.0, np.nan], "finite"),
        ([0.0, -1.0, -1.0], "duplicate"),
        ([-1.0, -2.0], "include 0"),
    ],
)
def test_time_grid_validation_fails_early(time, message):
    with pytest.raises(ValueError, match=message):
        orbit_maker.normalize_time_grid(time)


def test_time_grid_accepts_python_lists_and_preserves_order():
    normalized = orbit_maker.normalize_time_grid([2, -1, 0, 1])

    assert isinstance(normalized, np.ndarray)
    assert normalized.tolist() == [2.0, -1.0, 0.0, 1.0]


def test_orbit_integration_maps_permuted_times_back_to_requested_order():
    class _Cartesian:
        def __init__(self, values):
            self.x = type("Value", (), {"value": values})()

    class _Frame:
        def __init__(self, values):
            self.cartesian = _Cartesian(values)

    class _Coords:
        def __init__(self, values):
            self.galactic = _Frame(values)

    class _Orbit:
        integration_calls = []

        def __deepcopy__(self, _memo):
            return type(self)()

        def integrate(self, time, _potential):
            type(self).integration_calls.append(time.to_value(u.Myr).tolist())

        def SkyCoord(self, time):
            return _Coords(time.to_value(u.Myr))

    requested = np.array([2.0, -1.0, 0.0, 1.0])
    values = orbit_maker._integrate_orbit_in_requested_order(_Orbit(), requested, object())
    mapped = values(lambda coords: coords.galactic.cartesian.x.value)

    assert mapped.tolist() == requested.tolist()
    assert _Orbit.integration_calls == [[0.0, -1.0], [0.0, 1.0, 2.0]]


def test_cartesian_and_cylindrical_galactocentric_z_are_not_shadowed(monkeypatch):
    monkeypatch.setattr(orbit_maker, "coordFIX_to_coordROT", lambda frame, **_kwargs: frame)
    coords = _coordinate_groups([0.0, -1.0])

    trace = _trace()
    trace.cluster_int_coords = coords
    trace_frame = trace.create_integrated_dataframe([0.0, -1.0])

    arm = object.__new__(SpiralArmTrace)
    arm.cluster_int_coords = coords
    arm.df = pd.DataFrame({"name": ["arm-point"]})
    arm.arm_name = "Test Arm"
    arm_frame = arm._create_integrated_dataframe([0.0, -1.0])

    assert trace_frame["z_gc"].tolist() == [300.0, 300.0]
    assert trace_frame["z_gc_cyl"].tolist() == [500.0, 500.0]
    assert arm_frame["z_gc"].tolist() == [300.0, 300.0]
    assert arm_frame["z_gc_cyl"].tolist() == [500.0, 500.0]


def test_repeated_make_plot_reintegrates_changed_grid_and_refreshes_sizes(monkeypatch):
    integration_calls = []
    rotation_calls = []

    def fake_create_orbit(_coordinates, time, **kwargs):
        integration_calls.append((np.asarray(time, dtype=float).tolist(), kwargs["ro"], kwargs["vo"]))
        return _coordinate_groups(time)

    real_rotation = orbit_maker.coordFIX_to_coordROT

    def record_rotation(frame, *, r_sun=8.122, v_sun=236):
        rotation_calls.append((r_sun, v_sun))
        return real_rotation(frame, r_sun=r_sun, v_sun=v_sun)

    monkeypatch.setattr(orbit_maker, "create_orbit", fake_create_orbit)
    monkeypatch.setattr(orbit_maker, "coordFIX_to_coordROT", record_rotation)
    monkeypatch.setattr(
        orbit_maker,
        "get_center_orbit_coords",
        lambda time, *_args, **_kwargs: tuple(np.zeros(len(time)) for _ in range(3)),
    )

    collection = TraceCollection([_trace()])
    scene = Animate3D(collection, figure_theme="dark", ro=7.5, vo=210.0)
    scene.make_plot(time=[0.0, -1.0], show=False, show_gc_line=False)
    first_sizes = collection.get_cluster(0).df_int["size"].copy()
    scene.make_plot(time=[-2.0, 0.0, -1.0], show=False, show_gc_line=False)

    assert integration_calls == [([0.0, -1.0], 7.5, 210.0), ([-2.0, 0.0, -1.0], 7.5, 210.0)]
    assert rotation_calls == [(7.5, 210.0), (7.5, 210.0)]
    assert collection.time.tolist() == [-2.0, 0.0, -1.0]
    assert collection.get_cluster(0).df_int["time"].tolist() == [-2.0, 0.0, -1.0]
    assert "size" in collection.get_cluster(0).df_int
    assert len(first_sizes) == 2
    assert len(collection.get_cluster(0).df_int["size"]) == 3
