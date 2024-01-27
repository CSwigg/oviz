import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import pandas as pd

df = pd.read_csv('/Users/cam/Downloads/Hunt+2023 and Gagné+2018 (< 70 Myr)_integrated.csv')
df = df[['time', 'name', 'age_myr', 'x', 'y', 'z']]

app = dash.Dash(__name__)

app.layout = html.Div([
    dcc.Graph(id='2d-scatter', style={'width': '50%', 'display': 'inline-block'}),
    dcc.Graph(id='3d-scatter', style={'width': '50%', 'display': 'inline-block'}),
    dcc.RangeSlider(
        id='age-slider',
        min=df['age_myr'].min(),
        max=df['age_myr'].max(),
        value=[df['age_myr'].min(), df['age_myr'].max()],
        marks={str(age): str(age) for age in range(0, int(df['age_myr'].max()) + 1, 5)},
        step=None
    ),
])

def create_2d_scatter(df):
    fig = go.Figure(data=go.Scattergl(
        x = df[df['time'] == 0]['x'],
        y = df[df['time'] == 0]['y'],
        mode='markers'
    ))
    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        xaxis=dict(range=[-1000, 1000]),
        yaxis=dict(range=[-1000, 1000]),
        xaxis_title="X",
        yaxis_title="Y",
    )
    return fig

def create_3d_scatter(df):
    unique_times = sorted(df['time'].unique())
    data = []
    for time in unique_times:
        df_time = df[df['time'] == time]
        marker_size = [0 if age < abs(time) else 10 for age in df_time['age_myr']]

        scatter = go.Scatter3d(
            x=df_time['x'],
            y=df_time['y'],
            z=df_time['z'],
            mode='markers',
            marker=dict(
                size=marker_size,
                sizemode='diameter',
                line=dict(width=0)  # set line width to 0
            ),
            visible=(time == 0)
        )
        data.append(scatter)

    sliders = [dict(
        steps=[dict(
            method='restyle',
            args=['visible', [time == t for t in unique_times]],
            label=str(time)
        ) for time in unique_times],
        active=unique_times.index(0),
        pad={"t": 50}
    )]

    layout = go.Layout(
        sliders=sliders,
        scene=dict(
            xaxis=dict(range=[-1000, 1000]),
            yaxis=dict(range=[-1000, 1000]),
            zaxis=dict(range=[-300, 300]),
            aspectratio=dict(x=1, y=1, z=0.3),
            camera=dict(
                up=dict(x=0, y=1, z=0),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0., y=-.8, z=1.5),
                projection=dict(type='perspective')
            ),
            dragmode='turntable'
        )
    )

    fig = go.Figure(data=data, layout=layout)

    return fig

@app.callback(
    [Output('2d-scatter', 'figure'),
     Output('3d-scatter', 'figure')],
    [Input('age-slider', 'value'),
     Input('2d-scatter', 'selectedData')]
)
def update_plots(age_range, selectedData):
    filtered_df = df[(df['age_myr'] >= age_range[0]) & (df['age_myr'] <= age_range[1])]
    if selectedData:
        selected_points = [p['pointNumber'] for p in selectedData['points']]
        filtered_df = filtered_df.iloc[selected_points]

    return create_2d_scatter(filtered_df), create_3d_scatter(filtered_df)


if __name__ == '__main__':
    app.run_server(debug=True)