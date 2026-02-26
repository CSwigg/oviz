import numpy as np
import plotly.graph_objects as go
import copy

from oviz.app import (
    _extract_selected_point,
    _apply_footprint_to_figure,
    _build_aladin_catalog_payload,
    _filter_figure_by_age,
    create_dash_app,
    FOOTPRINT_TRACE_NAMES,
)


def _base_fig_dict():
    customdata = np.array([
        [5.0, 4.0, 120.0, -20.0, 150.0, 100.0, 50.0, 20.0],
        [15.0, 14.0, 125.0, -18.0, 190.0, 140.0, 40.0, 10.0],
    ])
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=[1.0, 2.0],
                y=[3.0, 4.0],
                z=[5.0, 6.0],
                mode="markers",
                customdata=customdata,
                meta={"trace_kind": "cluster"},
                name="Cluster A",
            )
        ]
    )
    return fig.to_dict()


def test_extract_selected_point_reads_customdata_schema():
    click_data = {
        "points": [
            {
                "customdata": [5.0, 4.0, 120.0, -20.0, 150.0, 100.0, 50.0, 20.0],
                "x": 1.0,
                "y": 3.0,
                "z": 5.0,
            }
        ]
    }
    selection = _extract_selected_point(click_data)
    assert selection is not None
    assert np.isclose(selection["l_deg"], 120.0)
    assert np.isclose(selection["b_deg"], -20.0)
    assert np.isclose(selection["x0"], 100.0)


def test_filter_figure_by_age_keeps_customdata_shape_consistent():
    fig_dict = _base_fig_dict()
    filtered = _filter_figure_by_age(fig_dict, age_min=0.0, age_max=10.0)
    trace = filtered["data"][0]
    assert len(trace["x"]) == 1
    assert len(trace["customdata"]) == 1
    assert np.isclose(trace["customdata"][0][0], 5.0)


def test_apply_footprint_adds_overlay_traces():
    fig_dict = _base_fig_dict()
    selection = {
        "l_deg": 120.0,
        "b_deg": -20.0,
        "dist_pc": 150.0,
        "x0": 100.0,
        "y0": 50.0,
        "z0": 20.0,
    }
    themed = {"footprint": "#0088ff", "axis": "#ffffff"}
    updated = _apply_footprint_to_figure(fig_dict, selection, sky_radius_deg=1.0, theme_colors=themed)
    names = {tr.get("name") for tr in updated["data"] if isinstance(tr, dict)}
    assert FOOTPRINT_TRACE_NAMES.issubset(names)
    assert "__oviz_sky_footprint_axis__" not in names


def test_apply_footprint_can_hide_overlay_with_visibility_flag():
    fig_dict = _base_fig_dict()
    trace0 = copy.deepcopy(fig_dict["data"][0])
    trace1 = copy.deepcopy(fig_dict["data"][0])
    cd0 = np.asarray(trace0["customdata"], dtype=float)
    cd1 = np.asarray(trace1["customdata"], dtype=float)
    cd0[:, 1] = cd0[:, 0] + 0.0
    cd1[:, 1] = cd1[:, 0] + 5.0
    trace0["customdata"] = cd0.tolist()
    trace1["customdata"] = cd1.tolist()

    fig_dict["frames"] = [
        {"name": "0.0", "data": [trace0]},
        {"name": "5.0", "data": [trace1]},
    ]
    fig_dict["layout"]["sliders"] = [{"steps": [{"label": "0.0"}, {"label": "5.0"}]}]

    selection = {
        "l_deg": 120.0,
        "b_deg": -20.0,
        "dist_pc": 150.0,
        "x0": 100.0,
        "y0": 50.0,
        "z0": 20.0,
    }
    themed = {"footprint": "#0088ff", "axis": "#ffffff"}
    updated = _apply_footprint_to_figure(
        fig_dict,
        selection,
        sky_radius_deg=1.0,
        theme_colors=themed,
    )
    frame0_cone = [tr for tr in updated["frames"][0]["data"] if tr.get("name") == "__oviz_sky_footprint_cone__"][0]
    frame1_cone = [tr for tr in updated["frames"][1]["data"] if tr.get("name") == "__oviz_sky_footprint_cone__"][0]
    assert len(frame0_cone.get("x", [])) > 0
    assert len(frame1_cone.get("x", [])) == 0


def test_create_dash_app_with_sky_panel():
    fig_dict = _base_fig_dict()
    fig = go.Figure(fig_dict)
    app = create_dash_app(
        figure=fig,
        enable_age_filter=True,
        enable_sky_panel=True,
        sky_radius_deg=1.0,
        sky_frame="galactic",
        sky_survey="P/DSS2/color",
    )
    assert app is not None


def test_catalog_payload_only_includes_selected_cluster_members():
    selection = {
        "l_deg": 120.0,
        "b_deg": -20.0,
        "dist_pc": 150.0,
        "x0": 100.0,
        "y0": 50.0,
        "z0": 20.0,
        "cluster_name": "Cluster A",
        "cluster_color": "#ff0000",
    }
    centers = {
        "Cluster A": {"l_deg": 120.0, "b_deg": -20.0, "color": "#00ff00"},
        "Cluster B": {"l_deg": 200.0, "b_deg": 40.0, "color": "#0000ff"},
    }
    members = {
        "Cluster A": {
            "l_deg": np.array([119.8, 120.1], dtype=float),
            "b_deg": np.array([-19.9, -20.2], dtype=float),
            "ra_deg": np.array([10.0, 10.1], dtype=float),
            "dec_deg": np.array([5.0, 5.1], dtype=float),
        },
        "Cluster B": {
            "l_deg": np.array([200.0], dtype=float),
            "b_deg": np.array([40.0], dtype=float),
            "ra_deg": np.array([20.0], dtype=float),
            "dec_deg": np.array([6.0], dtype=float),
        },
    }

    payload = _build_aladin_catalog_payload(
        selection=selection,
        sky_radius_deg=1.0,
        centers_by_name=centers,
        members_by_cluster=members,
    )

    assert len(payload) == 1
    assert payload[0]["name"] == "Cluster A"
    assert payload[0]["sourceSize"] == 7
    assert all(pt["label"] == "Cluster A" for pt in payload[0]["points"])
