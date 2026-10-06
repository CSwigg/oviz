"""Published spiral arms and the Khalil et al. (2025) spiral potential."""

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from galpy.potential import MWPotential2014

from oviz import Animate3D, Trace, TraceCollection, orbit_maker
from oviz.spiral_models import (
    CASTRO_GINARD2021_ARMS,
    CASTRO_GINARD2021_ARM_TABLE,
    KHALIL2025_ARMS,
    KHALIL2025_SPIRAL_MODES,
    KMS_KPC_TO_RAD_MYR,
    khalil2025_potential,
    khalil2025_spiral_potentials,
    spiral_arm_coordinates,
)

RO, VO = 8.122, 236.0


def spiback_spiral(mode, radius, phi, time_myr):
    """SPIBACK's released planar spiral term (no cutoffs, fully grown), (km/s)^2."""
    m, pitch = mode["m"], mode["pitch_rad"]
    k = m / (radius * np.sin(pitch))
    kh = k * 0.130
    d = (1 + kh + 0.3 * kh * kh) / (1 + 0.3 * kh)
    phase = m * (phi - mode["phase_rad"] - mode["pattern_speed"] * KMS_KPC_TO_RAD_MYR * time_myr
                 + np.log(radius / 8.275) / np.tan(pitch))
    return -mode["amplitude_kms2"] / (8.275 * k * d) * np.cos(phase)


def potential_value(potential, radius, phi, time_myr):
    return float(potential(radius * u.kpc, 0 * u.kpc, phi=phi * u.rad, t=time_myr * u.Myr))


def test_khalil_spirals_match_the_released_spiback_potential_in_the_plane():
    for mode, spiral in zip(KHALIL2025_SPIRAL_MODES, khalil2025_spiral_potentials()):
        for time_myr in (0.0, -17.0, -60.0):
            for radius in (4.0, 6.5, RO, 11.0, 17.0):
                for phi in np.linspace(-np.pi, np.pi, 13):
                    expected = spiback_spiral(mode, radius, phi, time_myr) * np.exp(-(radius - RO) / 1000.0)
                    assert potential_value(spiral, radius, phi, time_myr) == pytest.approx(expected, rel=1e-9, abs=1e-9)
    # Troughs at the Sun: 167.6 and 70.7 (km/s)^2 deep (the m = 2 mode is the stronger).
    depths = [min(potential_value(s, RO, p, 0.0) for p in np.linspace(-np.pi, np.pi, 3601)) for s in khalil2025_spiral_potentials()]
    assert depths == pytest.approx([-167.63, -70.69], abs=0.05)


def test_khalil_potential_is_mwpotential2014_plus_the_two_spirals():
    potential = khalil2025_potential()
    assert potential[:len(MWPotential2014)] == list(MWPotential2014)
    assert [p.isNonAxi for p in potential[len(MWPotential2014):]] == [True, True]


def test_khalil_arms_are_the_troughs_of_the_potential_at_every_time():
    spirals = dict(zip((2, 3), khalil2025_spiral_potentials()))
    for time_myr in (0.0, -30.0, -60.0):
        curves = KHALIL2025_ARMS.curves(time_myr)
        assert len(curves) == 5  # two m = 2 arms, three m = 3 arms
        for name, radius, phi in curves:
            m = int(name[2])
            for i in range(0, len(radius), 50):
                offsets = np.radians(np.linspace(-4.0, 4.0, 161)) / m
                values = [potential_value(spirals[m], radius[i], phi[i] + d, time_myr) for d in offsets]
                assert abs(offsets[int(np.argmin(values))]) < np.radians(0.06) / m


def test_arms_trail_and_turn_at_their_published_pattern_speeds():
    for model, speeds in (
        (KHALIL2025_ARMS, [m["pattern_speed"] for m in KHALIL2025_SPIRAL_MODES for _ in range(m["m"])]),
        (CASTRO_GINARD2021_ARMS, [row[-1] for row in CASTRO_GINARD2021_ARM_TABLE]),
    ):
        now, past = model.curves(0.0), model.curves(-40.0)
        for (name, radius, phi), (_, radius_past, phi_past), speed in zip(now, past, speeds):
            # Trailing: along an arm, azimuth falls as radius grows (rotation is +phi).
            assert np.polyfit(radius, np.unwrap(phi), 1)[0] < 0, name
            np.testing.assert_allclose(radius_past, radius)
            np.testing.assert_allclose(phi_past - phi, -40.0 * speed * KMS_KPC_TO_RAD_MYR)
    # The Scutum segment turns fastest, Perseus slowest (Castro-Ginard et al., Table 2).
    assert CASTRO_GINARD2021_ARMS.pattern_speeds() == {"Perseus": 17.82, "Local": 33.76, "Sagittarius": 26.10, "Scutum": 49.81}


def test_castro_ginard_arms_sit_where_the_paper_put_them_today():
    helio, _ = spiral_arm_coordinates(CASTRO_GINARD2021_ARMS, 0.0)
    start = 0
    for arm, theta_ref, (t0, t1), r_ref, pitch, _speed in CASTRO_GINARD2021_ARM_TABLE:
        theta = np.radians(np.linspace(t0, t1, int(np.ceil((t1 - t0) / 0.5)) + 1))
        radius = r_ref * np.exp(-(theta - np.radians(theta_ref)) * np.tan(np.radians(pitch)))
        x, y, z = helio[:, start:start + theta.size]
        # Paper frame (Sun at 8.178 kpc): X towards the centre, Y along rotation.
        np.testing.assert_allclose(x, (8.178 - radius * np.cos(theta)) * 1000.0, atol=0.2)
        np.testing.assert_allclose(y, radius * np.sin(theta) * 1000.0, atol=0.2)
        assert np.all(np.abs(z) < 40.0), arm  # the Galactic plane, a few tens of pc below the Sun
        start += theta.size + 1
    assert start - 1 == helio.shape[1]
    first_arm = len(CASTRO_GINARD2021_ARMS.curves(0.0)[0][1])
    assert np.isnan(helio[:, first_arm]).all()  # a gap between arms


def test_arm_lines_keep_their_vertices_through_time():
    for model in (KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS):
        shapes = {spiral_arm_coordinates(model, t)[0].shape for t in (0.0, -1.0, -33.5, -60.0)}
        gaps = {tuple(np.flatnonzero(np.isnan(spiral_arm_coordinates(model, t)[0][0]))) for t in (0.0, -60.0)}
        assert len(shapes) == 1 and len(gaps) == 1


def test_the_lsr_frame_stays_circular_in_a_spiral_potential():
    time = np.arange(0.0, -61.0, -5.0)
    plain = orbit_maker.get_center_orbit_coords(time, None, MWPotential2014)
    spiral = orbit_maker.get_center_orbit_coords(time, None, khalil2025_potential())
    np.testing.assert_array_equal(np.asarray(spiral), np.asarray(plain))
    assert orbit_maker.reference_potential(khalil2025_potential(), None) == list(MWPotential2014)
    # A focus group's orbit is an object, so it feels the arms.
    custom = [0.0, 0.0, 0.0, -11.1, -12.24, -7.25]
    full = khalil2025_potential()
    assert orbit_maker.reference_potential(full, custom) is full
    assert orbit_maker.reference_potential(MWPotential2014, None) is MWPotential2014


def test_make_plot_adds_one_hidden_trace_per_arm_model():
    from oviz.viewer.compile import compile_scene_spec

    clusters = Trace(
        pd.DataFrame({
            "x": [100.0, -400.0], "y": [50.0, 300.0], "z": [10.0, -20.0],
            "U": [-10.0, -12.0], "V": [-15.0, -21.0], "W": [-7.0, -5.0],
            "name": ["A", "B"], "age_myr": [10.0, 30.0],
        }),
        data_name="Clusters",
    )
    time = np.arange(0.0, -11.0, -1.0)
    scene = Animate3D(data_collection=TraceCollection([clusters]), figure_theme="dark", potential=khalil2025_potential())
    figure = scene.make_plot(time=time, galactic_mode=True, spiral_arm_models=(KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS))
    spec = figure.to_dict()
    names = {KHALIL2025_ARMS.name, CASTRO_GINARD2021_ARMS.name}
    arm_traces = [[t for t in f["traces"] if t.get("name") in names] for f in spec["frames"]]
    assert all(len(traces) == 2 for traces in arm_traces)
    for name in names:
        counts = {len(t["segments"]) for traces in arm_traces for t in traces if t["name"] == name}
        assert len(counts) == 1  # the same vertices in every frame
    keys = {t["key"] for t in arm_traces[0]}
    for visibility in spec["group_visibility"].values():
        assert all(visibility[key] == "legendonly" for key in keys)
    manifest = compile_scene_spec(spec).manifest
    compiled = {t["name"]: t for t in manifest["traces"] if t["name"] in names}
    assert set(compiled) == names
    for trace in compiled.values():
        assert trace["showInLegend"] and trace["skyHidden"] and trace.get("role") != "guide"
        assert trace["lines"]["position"]["frames"] == len(time)  # they move with time
    assert compiled[CASTRO_GINARD2021_ARMS.name]["lines"]["dash"] == "dash"
