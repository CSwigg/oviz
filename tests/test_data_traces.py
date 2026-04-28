import numpy as np
import pandas as pd

from oviz.traces import Layer, LayerCollection, Trace, TraceCollection


def _make_integrated_trace():
    df = pd.DataFrame(
        {
            "x": [10.0, 20.0],
            "y": [5.0, 15.0],
            "z": [1.0, 2.0],
            "U": [0.5, 1.5],
            "V": [2.5, 3.5],
            "W": [4.5, 5.5],
            "name": ["young", "old"],
            "age_myr": [5.0, 25.0],
        }
    )
    trace = Trace(df, data_name="Clusters")
    time = np.array([0.0, -1.0])

    centered = (
        np.array([[1.0, 1.5], [2.0, 2.5]]),
        np.array([[3.0, 3.5], [4.0, 4.5]]),
        np.array([[5.0, 5.5], [6.0, 6.5]]),
    )
    helio = (
        np.array([[10.0, 10.5], [20.0, 20.5]]),
        np.array([[30.0, 30.5], [40.0, 40.5]]),
        np.array([[50.0, 50.5], [60.0, 60.5]]),
    )
    galactocentric = (
        np.array([[100.0, 100.5], [200.0, 200.5]]),
        np.array([[300.0, 300.5], [400.0, 400.5]]),
        np.array([[500.0, 500.5], [600.0, 600.5]]),
    )
    cylindrical = (
        np.array([[700.0, 700.5], [800.0, 800.5]]),
        np.array([[0.1, 0.2], [0.3, 0.4]]),
        np.array([[900.0, 900.5], [1000.0, 1000.5]]),
    )

    trace.cluster_int_coords = (centered, helio, galactocentric, cylindrical)
    trace.integrated = True
    trace.df_int = trace.create_integrated_dataframe(time)
    return trace


def test_limit_cluster_age_preserves_integrated_coordinate_structure():
    trace = _make_integrated_trace()

    trace.limit_cluster_age(0.0, 10.0)

    assert trace.df["name"].tolist() == ["young"]
    assert set(trace.df_int["name"]) == {"young"}
    assert len(trace.cluster_int_coords) == 4
    assert trace.cluster_int_coords[0][0].shape == (1, 2)
    assert trace.cluster_int_coords[1][1].shape == (1, 2)
    assert trace.cluster_int_coords[2][2].shape == (1, 2)
    assert trace.cluster_int_coords[3][0].shape == (1, 2)


def test_limit_cluster_by_name_preserves_integrated_coordinate_structure():
    trace = _make_integrated_trace()

    trace.limit_cluster_by_name(["old"])

    assert trace.df["name"].tolist() == ["old"]
    assert set(trace.df_int["name"]) == {"old"}
    assert len(trace.cluster_int_coords) == 4
    assert trace.cluster_int_coords[0][0].shape == (1, 2)
    assert trace.cluster_int_coords[1][1].shape == (1, 2)
    assert trace.cluster_int_coords[2][2].shape == (1, 2)
    assert trace.cluster_int_coords[3][0].shape == (1, 2)


def test_trace_collection_uses_independent_default_cluster_lists():
    first = TraceCollection()
    second = TraceCollection()

    first.add_cluster(_make_integrated_trace())

    assert len(first.get_all_clusters()) == 1
    assert second.get_all_clusters() == []


def test_layer_can_create_static_integrated_state_without_velocities():
    df = pd.DataFrame(
        {
            "x": [1.0, 2.0],
            "y": [3.0, 4.0],
            "z": [5.0, 6.0],
        }
    )
    layer = Layer(df, layer_name="Markers", assume_stationary=True)

    layer.integrate_orbits(np.array([0.0, -1.0, -2.0]))

    assert layer.supports_orbit_tracing is False
    assert layer.integrated is True
    assert set(layer.df_int["name"]) == {"Markers"}
    assert set(layer.df_int["age_myr"]) == {0.0}
    assert np.allclose(layer.df_int["x"].to_numpy()[:3], np.array([1.0, 1.0, 1.0]))


def test_layer_collection_exposes_layer_alias_methods():
    df = pd.DataFrame(
        {
            "x": [0.0],
            "y": [0.0],
            "z": [0.0],
            "U": [0.0],
            "V": [0.0],
            "W": [0.0],
            "name": ["A"],
            "age_myr": [1.0],
        }
    )
    layer = Layer(df, layer_name="Layer A")
    collection = LayerCollection([layer])

    assert collection.layers == collection.get_all_clusters()
    assert collection.get_layer("Layer A") is layer
    assert collection.get_all_layers() == [layer]
