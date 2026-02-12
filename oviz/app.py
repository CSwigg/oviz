from __future__ import annotations

from typing import Optional

from dash import Dash, dcc, html
import plotly.graph_objects as go


def create_dash_app(
    figure: go.Figure,
    title: str = "oviz",
    graph_id: str = "oviz-graph",
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

    app.layout = html.Div(
        [
            dcc.Graph(
                id=graph_id,
                figure=figure,
                style={"height": "95vh"},
                config={"displaylogo": False, "responsive": True},
            )
        ],
        style={"width": "100%", "height": "100vh", "margin": "0", "padding": "0"},
    )
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
