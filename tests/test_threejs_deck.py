from oviz.threejs_deck import DECK_SCHEMA_VERSION, normalize_deck_spec
from oviz.threejs_runtime_deck import THREEJS_DECK_RUNTIME_JS


def test_normalize_deck_spec_preserves_order_and_state_links():
    deck = normalize_deck_spec(
        {
            "enabled": True,
            "slides": [
                {
                    "id": "opening",
                    "name": "Opening",
                    "state_id": "original",
                    "blocks": [
                        {
                            "id": "title",
                            "type": "title",
                            "text": "Young clusters",
                            "x": -20,
                            "y": 130,
                            "width": 500,
                            "font_size": 500,
                            "font_weight": 2000,
                            "color": "#ABCDEF",
                            "align": "center",
                        }
                    ],
                },
                {"id": "science", "state_id": "state-sky", "blocks": []},
            ],
        }
    )

    assert deck["schema_version"] == DECK_SCHEMA_VERSION
    assert [slide["id"] for slide in deck["slides"]] == ["opening", "science"]
    assert deck["slides"][0]["state_id"] is None
    assert deck["slides"][1]["state_id"] == "state-sky"
    title = deck["slides"][0]["objects"][0]
    assert title["kind"] == "text"
    assert title["x"] == 0
    assert title["y"] == 880
    assert title["width"] == 1600
    assert title["font_size"] == 240
    assert title["font_weight"] == 900
    assert title["color"] == "#abcdef"


def test_normalize_deck_spec_repairs_duplicate_ids_and_invalid_styles():
    deck = normalize_deck_spec(
        {
            "slides": [
                {
                    "id": "slide",
                    "blocks": [
                        {"id": "block", "type": "unknown", "align": "diagonal"},
                        {
                            "id": "block",
                            "type": "subtitle",
                            "font_family": "serif",
                            "font_style": "italic",
                            "underline": True,
                            "line_height": 9,
                            "list_style": "bullet",
                        },
                    ],
                },
                {"id": "slide", "blocks": []},
            ]
        }
    )

    assert deck["enabled"] is True
    assert len({slide["id"] for slide in deck["slides"]}) == 2
    blocks = deck["slides"][0]["objects"]
    assert len({block["id"] for block in blocks}) == 2
    assert blocks[0]["type"] == "text"
    assert blocks[0]["align"] == "left"
    assert blocks[1]["font_size"] == 34
    assert blocks[1]["font_family"] == "serif"
    assert blocks[1]["font_style"] == "italic"
    assert blocks[1]["underline"] is True
    assert blocks[1]["line_height"] == 2.4
    assert blocks[1]["list_style"] == "bullet"


def test_empty_deck_is_disabled_and_backwards_compatible():
    deck = normalize_deck_spec(None)

    assert deck["enabled"] is False
    assert deck["slides"] == []
    assert deck["reveal"]["version"] == "5.2.1"
    assert deck["guides"] == {"smart": True, "grid": False, "grid_size": 20.0}
    assert deck["available"] is True


def test_deck_can_be_explicitly_unavailable():
    deck = normalize_deck_spec({"available": False, "enabled": False, "slides": []})

    assert deck["available"] is False
    assert deck["enabled"] is False


def test_v2_shape_object_and_guide_settings_round_trip():
    deck = normalize_deck_spec(
        {
            "schema_version": 2,
            "guides": {"smart": False, "grid": True, "grid_size": 40},
            "slides": [
                {
                    "id": "shapes",
                    "objects": [
                        {
                            "id": "shape-one",
                            "name": "Result callout",
                            "kind": "shape",
                            "shape_type": "rounded_rectangle",
                            "x": 100,
                            "y": 200,
                            "width": 420,
                            "height": 160,
                            "rotation": 12,
                            "opacity": 0.8,
                            "fill_color": "#123456",
                            "border_style": "dashed",
                            "corner_radius": 36,
                            "text": "Important result",
                        }
                    ],
                }
            ],
        }
    )

    shape = deck["slides"][0]["objects"][0]
    assert deck["schema_version"] == 2
    assert deck["guides"] == {"smart": False, "grid": True, "grid_size": 40.0}
    assert shape["kind"] == "shape"
    assert shape["shape_type"] == "rounded_rectangle"
    assert shape["width"] == 420
    assert shape["height"] == 160
    assert shape["rotation"] == 12
    assert shape["fill_color"] == "#123456"
    assert shape["border_style"] == "dashed"


def test_keynote_lite_runtime_exposes_object_tools_and_snap_modifiers():
    runtime = THREEJS_DECK_RUNTIME_JS

    assert "const OVIZ_DECK_VERSION = 2;" in runtime
    assert '"rectangle", "rounded_rectangle", "ellipse", "triangle", "diamond", "line", "arrow"' in runtime
    assert "function ovizDeckRenderSelectionFrame()" in runtime
    assert '["nw", "n", "ne", "e", "se", "s", "sw", "w"]' in runtime
    assert "const thresholdX = 6 / metrics.scaleX;" in runtime
    assert "const hysteresisX = 8 / metrics.scaleX;" in runtime
    assert "event.metaKey || event.ctrlKey" in runtime
    assert "event.altKey" in runtime
    assert "ovizDeckDuplicateObjects" in runtime
    assert "ovizDeckGroupObjects" in runtime
    assert "ovizDeckDistributeObjects" in runtime
    assert "addShape: ovizDeckAddShape" in runtime
    assert "deck-selection-changed" in runtime
    assert "deck-object-changed" in runtime
