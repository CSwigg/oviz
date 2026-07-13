from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_chronos_july4_actions_audit.py")
SPEC = importlib.util.spec_from_file_location("oviz_july4_actions_audit", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
AUDIT_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT_MODULE)
build_audit_scene = AUDIT_MODULE.build_audit_scene
render_audit_html = AUDIT_MODULE.render_audit_html


def _audit_fixture_scene() -> dict:
    frames = []
    for index in range(121):
        trace_count = 23 if index == 120 else 19
        frames.append(
            {
                "time": float(index - 120),
                "traces": [
                    {
                        "key": f"trace-{trace_index}",
                        "name": "Clusters < 15 Myr" if trace_index == 0 else f"Trace {trace_index}",
                        "showlegend": True,
                        "points": [
                            {
                                "x": float(point_index),
                                "y": float(trace_index),
                                "z": 0.0,
                                "selection": {
                                    "cluster_name": f"Cluster {point_index}",
                                    "trace_name": f"Trace {trace_index}",
                                },
                            }
                            for point_index in range(16)
                        ] if trace_index == 0 else [],
                    }
                    for trace_index in range(trace_count)
                ],
            }
        )
    return {
        "width": 640,
        "height": 480,
        "center": {"x": 0.0, "y": 0.0, "z": 0.0},
        "max_span": 10000.0,
        "frames": frames,
        "group_order": ["Clusters"],
        "default_group": "Clusters",
        "group_visibility": {
            "Clusters": {f"trace-{index}": True for index in range(23)},
        },
        "initial_state": {
            "current_group": "Clusters",
            "current_frame_index": 120,
            "current_frame_value": 120.0,
            "global_controls": {},
        },
    }


def test_audit_scene_seeds_deterministic_states_actions_and_topology_boundary():
    scene = build_audit_scene(_audit_fixture_scene())

    assert [item["id"] for item in scene["states"]["items"]] == [
        "audit-lookback",
        "audit-topology-boundary",
        "audit-present",
    ]
    assert len(scene["actions"]["items"]) == 4
    assert len(scene["frames"][119]["traces"]) == 19
    assert len(scene["frames"][120]["traces"]) == 23
    concurrent = scene["actions"]["items"][-1]["steps"]
    assert concurrent[0]["start"] == "after_previous"
    assert concurrent[1]["start"] == "with_previous"
    lookback = scene["states"]["items"][0]["snapshot"]
    assert lookback["current_selection_mode"] == "lasso"
    assert len(lookback["selected_cluster_keys"]) == 12
    assert lookback["lasso_selection_mask"]["data_url"].startswith("data:image/png;base64,")
    assert scene["debug_transitions"] is True
    assert scene["initial_state"]["minimal_mode_enabled"] is True


def test_audit_html_contains_self_running_diagnostics_and_unified_runtime():
    html = render_audit_html(_audit_fixture_scene(), auto_run=True)

    assert 'id="oviz-actions-audit-panel"' in html
    assert "Run deterministic audit" in html
    assert "OvizTransitionRuntime.createFrameCoordinator" in html
    assert "updateOvizUnifiedTransitionSession(now)" in html
    assert 'button.dataset.actionKey = String(action.key || "")' in html
    assert "if (true) setTimeout(run, 0)" in html
    assert "await states.previous()" in html
    assert 'clickAction("Audit: State Present")' in html
    assert "State to Action did not cancel the State" in html
    assert "Rapid State 1 was not cancelled" in html
    assert "all-domains lasso membership" in html
    assert "renderedTraceFidelity" in html
    assert "actionTransitionMetrics" not in html
    panel_index = html.rfind('id="oviz-actions-audit-panel"')
    assert panel_index > html.rfind("</script>", 0, panel_index)
    assert panel_index < html.rfind("</body>")
