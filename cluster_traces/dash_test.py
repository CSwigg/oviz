import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objects as go

def create_dash_app(collection, plotter):
    def make_plot(collection, plotter):
        plotter.generate_3d_plot(collection)
        fig = plotter.figure
        return fig

    def make_2d_plot(collection, ages):
        cluster_groups = collection.get_clusters_in_age_range(ages)
        scatter_list = []
        for cg in cluster_groups:
            df = cg.df
            x = df['x']
            y = df['y']
            scatter = go.Scatter(x=x, y=y, mode='markers', marker=dict(size=5, color=cg.color, opacity=cg.opacity))
            scatter_list.append(scatter)
        fig = go.Figure(data=scatter_list)
        return fig

    app = dash.Dash(__name__)

    min_age = collection.get_min_age()
    max_age = collection.get_max_age()

    app.layout = html.Div([
        html.Div(dcc.RangeSlider(id='age-slider', min=min_age, max=max_age, value=[min_age, max_age]), style={'width': '100%'}),
        html.Div(id='2d-plot', style={'width': '50%', 'display': 'inline-block'}),
        html.Div(id='3d-plot', style={'width': '50%', 'display': 'inline-block'}),
    ])

    @app.callback(
        [Output('2d-plot', 'children'), Output('3d-plot', 'children')],
        [Input('age-slider', 'value')]
    )
    def update_plots(ages):
        fig_2d = make_2d_plot(collection, ages)
        fig_3d = make_plot(collection, plotter)
        return dcc.Graph(figure=fig_2d), dcc.Graph(figure=fig_3d)

    return app