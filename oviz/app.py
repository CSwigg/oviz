from __future__ import annotations

from typing import Optional
import copy

from dash import Dash, dcc, html, Input, Output, State
import plotly.graph_objects as go
import numpy as np


def _is_array_like(value):
    return isinstance(value, (list, tuple, np.ndarray))


def _extract_age_values_from_trace(trace: dict):
    customdata = trace.get("customdata")
    if customdata is None:
        return None

    arr = np.asarray(customdata)
    if arr.ndim != 2 or arr.shape[1] < 1:
        return None

    try:
        return arr[:, 0].astype(float)
    except Exception:
        return None


def _iter_cluster_traces(fig_dict: dict):
    for trace in fig_dict.get("data", []):
        if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
            yield trace

    for frame in fig_dict.get("frames", []):
        for trace in frame.get("data", []):
            if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
                yield trace


def _get_age_bounds(fig_dict: dict):
    all_ages = []
    for trace in _iter_cluster_traces(fig_dict):
        ages = _extract_age_values_from_trace(trace)
        if ages is not None and ages.size:
            all_ages.append(ages)

    if not all_ages:
        return None

    flat = np.concatenate(all_ages)
    if flat.size == 0:
        return None

    age_min = float(np.nanmin(flat))
    age_max = float(np.nanmax(flat))
    return age_min, age_max


def _theme_colors_from_figure(fig_dict: dict):
    scene = fig_dict.get("layout", {}).get("scene", {})
    xaxis = scene.get("xaxis", {}) if isinstance(scene, dict) else {}
    bg = str(xaxis.get("backgroundcolor", "")).lower()

    # Heuristic aligned with oviz themes: dark scenes use bright gray UI,
    # lighter scenes use dark UI elements.
    if "rgb(0" in bg or "rgba(0" in bg or "black" in bg:
        return {
            "text": "#b7b7b7",
            "panel_bg": "rgba(30,30,30,0.75)",
            "rail": "rgba(110,110,110,0.40)",
            "track": "#b7b7b7",
            "handle_border": "#d0d0d0",
            "handle_bg": "#6d6d6d",
        }

    return {
        "text": "#1f1f1f",
        "panel_bg": "rgba(255,255,255,0.82)",
        "rail": "rgba(25,25,25,0.25)",
        "track": "#2f2f2f",
        "handle_border": "#2f2f2f",
        "handle_bg": "#f2f2f2",
    }


def _apply_mask_to_trace(trace: dict, mask: np.ndarray):
    n = len(mask)

    for key in ("x", "y", "z", "text", "hovertext", "customdata"):
        val = trace.get(key)
        if _is_array_like(val) and len(val) == n:
            arr = np.asarray(val)
            trace[key] = arr[mask].tolist()

    marker = trace.get("marker")
    if isinstance(marker, dict):
        for mk in ("size", "color", "opacity", "symbol"):
            val = marker.get(mk)
            if _is_array_like(val) and len(val) == n:
                arr = np.asarray(val)
                marker[mk] = arr[mask].tolist()


def _filter_cluster_trace_by_age(trace: dict, age_min: float, age_max: float):
    ages = _extract_age_values_from_trace(trace)
    if ages is None:
        return

    mask = (ages >= float(age_min)) & (ages <= float(age_max))
    _apply_mask_to_trace(trace, mask)


def _filter_figure_by_age(fig_dict: dict, age_min: float, age_max: float):
    filtered = copy.deepcopy(fig_dict)

    for trace in filtered.get("data", []):
        if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
            _filter_cluster_trace_by_age(trace, age_min, age_max)

    for frame in filtered.get("frames", []):
        for trace in frame.get("data", []):
            if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
                _filter_cluster_trace_by_age(trace, age_min, age_max)

    return filtered


def create_dash_app(
    figure: go.Figure,
    title: str = "oviz",
    graph_id: str = "oviz-graph",
    enable_age_filter: bool = True,
) -> Dash:
    """Create a Dash app that renders an existing Plotly figure.

    Parameters
    ----------
    figure
        Plotly figure produced by Animate3D.make_plot(...).
    title
        Browser tab title.
    graph_id
        Component id for the main graph (useful for later callbacks).
    """
    app = Dash(__name__)
    app.title = title

    base_figure_dict = figure.to_dict()
    age_bounds = _get_age_bounds(base_figure_dict)
    theme_colors = _theme_colors_from_figure(base_figure_dict)

    controls = []
    if enable_age_filter and age_bounds is not None:
        age_min, age_max = age_bounds
        controls.append(
            html.Div(
                [
                    html.Div("Age (Myr)", style={"fontWeight": 600, "marginBottom": "2px", "fontSize": "12px"}),
                    dcc.RangeSlider(
                        id="oviz-age-range",
                        min=age_min,
                        max=age_max,
                        step=0.1,
                        value=[age_min, age_max],
                        allowCross=False,
                        marks={},
                        vertical=True,
                        verticalHeight=300,
                        tooltip={"placement": "bottom", "always_visible": False},
                        className="oviz-age-slider",
                    ),
                ],
                style={
                    "position": "absolute",
                    "right": "10px",
                    "top": "50%",
                    "transform": "translateY(-50%)",
                    "zIndex": 25,
                    "padding": "6px 6px",
                    "width": "62px",
                    "background": "transparent",
                    "borderRadius": "8px",
                    "display": "flex",
                    "flexDirection": "column",
                    "alignItems": "center",
                    "gap": "4px",
                    "color": theme_colors["text"],
                },
            )
        )

    slider_css = f"""
    .oviz-age-slider .rc-slider-rail {{
        background-color: {theme_colors['rail']};
        width: 4px;
    }}
    .oviz-age-slider .rc-slider-track {{
        background-color: {theme_colors['track']};
        width: 4px;
    }}
    .oviz-age-slider .rc-slider-handle {{
        border: 2px solid {theme_colors['handle_border']};
        background-color: {theme_colors['handle_bg']};
        box-shadow: none;
        width: 14px;
        height: 14px;
        margin-left: -5px;
    }}
    .oviz-age-slider .rc-slider-handle:hover,
    .oviz-age-slider .rc-slider-handle:active,
    .oviz-age-slider .rc-slider-handle:focus {{
        border-color: {theme_colors['handle_border']};
        box-shadow: none;
    }}
    .oviz-age-slider .rc-slider-mark,
    .oviz-age-slider .rc-slider-dot {{
        display: none;
    }}
    """

    # Dash compatibility: older versions do not expose html.Style component.
    # Inject custom CSS through index_string instead.
    app.index_string = f"""
    <!DOCTYPE html>
    <html>
        <head>
            {{%metas%}}
            <title>{{%title%}}</title>
            {{%favicon%}}
            {{%css%}}
            <style>{slider_css}</style>
        </head>
        <body>
            {{%app_entry%}}
            <footer>
                {{%config%}}
                {{%scripts%}}
                {{%renderer%}}
            </footer>
        </body>
    </html>
    """

    app.layout = html.Div(
        [
            dcc.Store(id="oviz-base-figure", data=base_figure_dict),
            dcc.Graph(
                id=graph_id,
                figure=figure,
                style={"height": "95vh"},
                config={"displaylogo": False, "responsive": True},
            ),
            *controls,
        ],
        style={"width": "100%", "height": "100vh", "margin": "0", "padding": "0", "position": "relative"},
    )

    if controls:
        @app.callback(
            Output(graph_id, "figure"),
            Input("oviz-age-range", "value"),
            State("oviz-base-figure", "data"),
            prevent_initial_call=False,
        )
        def _update_age_filtered_figure(age_range, base_figure):
            if not age_range or len(age_range) != 2:
                return copy.deepcopy(base_figure)

            age_min_val, age_max_val = age_range
            return _filter_figure_by_age(base_figure, age_min_val, age_max_val)

    return app


def run_dash_app(
    figure: go.Figure,
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
) -> Dash:
    """Run the figure as a standard Dash app and return the app instance."""
    app = create_dash_app(figure=figure, title=title)
    app.run(host=host, port=port, debug=debug)
    return app


def run_dash_app_in_notebook(
    figure: go.Figure,
    mode: str = "inline",
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
    height: int = 750,
) -> Dash:
    """Run a Dash app from a notebook.

    Uses Dash's native notebook support when available.

    mode
        One of: "inline", "external", "tab", "jupyterlab".
    """
    app = create_dash_app(figure=figure, title=title)

    # Dash >= 2.11 supports jupyter_mode directly on app.run.
    app.run(
        host=host,
        port=port,
        debug=debug,
        jupyter_mode=mode,
        jupyter_height=height,
    )
    return app


def launch_from_animate3d(
    animate3d,
    *,
    time,
    make_plot_kwargs: Optional[dict] = None,
    notebook: bool = True,
    mode: str = "inline",
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
    height: int = 750,
):
    """Build a figure via Animate3D.make_plot and launch a Dash app.

    This is a convenience helper for notebook workflows.
    """
    kwargs = dict(make_plot_kwargs or {})
    kwargs.setdefault("show", False)
    kwargs.setdefault("save_name", None)

    fig = animate3d.make_plot(time=time, **kwargs)
    if notebook:
        return run_dash_app_in_notebook(
            figure=fig,
            mode=mode,
            host=host,
            port=port,
            debug=debug,
            title=title,
            height=height,
        )

    return run_dash_app(
        figure=fig,
        host=host,
        port=port,
        debug=debug,
        title=title,
    )
