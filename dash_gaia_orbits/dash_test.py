import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import copy

def layout_2d_plot():
    """
    Create the layout for the 2D plot.

    Returns:
        go.Layout: The layout for the 2D plot.
    """
    layout = go.Layout(
        autosize=False,
        width=500,
        height=500,
        xaxis=dict(
            range=[-2000, 2000],
            showgrid=True,
            gridwidth=0.5,
            gridcolor='gray',
            zerolinecolor='gray',
            tickmode='linear',
            tick0=0,
            dtick=500
        ),
        yaxis=dict(
            range=[-2000, 2000],
            showgrid=True,
            gridwidth=0.5,
            gridcolor='gray',
            zerolinecolor='gray',
            tickmode='linear',
            tick0=0,
            dtick=500
        ),
        plot_bgcolor='black',
        paper_bgcolor='black',
        font=dict(color='white'),
        dragmode='select'
    )
    return layout


def create_dash_app(collection, plotter):
    """
    Create a Dash app with interactive plots.

    Args:
        collection (object): The collection object containing cluster data.
        plotter (object): The plotter object for generating 3D plots.

    Returns:
        dash.Dash: The Dash app.
    """
    def make_2d_plot(collection, age_min, age_max):
        """
        Create a 2D plot based on the given age range.

        Args:
            collection (object): The collection object containing cluster data.
            age_min (int): The minimum age for filtering clusters.
            age_max (int): The maximum age for filtering clusters.

        Returns:
            go.Figure: The 2D plot figure.
        """
        collection.limit_all_cluster_ages(age_min, age_max)
        cluster_groups = collection.get_all_clusters()

        scatter_list = []
        for cg in cluster_groups:
            df = cg.df
            x = df['x']
            y = df['y']
            scatter = go.Scatter(x=x, y=y, mode='markers', marker=dict(size=5, color=cg.color, opacity=cg.opacity), hovertext=cg.df['name'], name=cg.data_name)
            scatter_list.append(scatter)
        fig = go.Figure(data=scatter_list, layout=layout_2d_plot())
        return fig

    app = dash.Dash(__name__)
    # Create a copy of the original collection
    original_collection = copy.deepcopy(collection)

    app.layout = html.Div([
        html.Div(dcc.RangeSlider(id='age-slider', min=0, max=100, value=[0, 100]), style={'width': '100%'}),
        dcc.Graph(id='2d-plot', style={'width': '50%', 'display': 'inline-block'}),
        dcc.Graph(id='3d-plot', style={'width': '50%', 'display': 'inline-block'})
    ])

    @app.callback(
        [Output('2d-plot', 'figure'), Output('3d-plot', 'figure')],
        [Input('age-slider', 'value'), Input('2d-plot', 'selectedData')],
    )
    def update_plot(ages, selected_data):
        """
        Update the plots based on the selected age range and data points.

        Args:
            ages (list): The selected age range.
            selected_data (dict): The selected data points on the 2D plot.

        Returns:
            go.Figure: The updated 2D plot figure.
            go.Figure: The updated 3D plot figure.
        """
        age_min, age_max = ages
        collection = copy.deepcopy(original_collection)
        collection.limit_all_cluster_ages(age_min, age_max)

        if selected_data is not None:
            selected_points = selected_data['points']
            selected_names = [point['hovertext'] for point in selected_points]
            collection.limit_all_cluster_names(selected_names)

        fig_2d = make_2d_plot(collection, age_min, age_max)
        plotter.generate_3d_plot(collection, show=False)
        fig_3d = plotter.figure
        return fig_2d, fig_3d

    return app
