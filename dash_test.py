import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import pandas as pd

# ------------- Load Data -------------
df = pd.read_csv('/Users/cam/Downloads/hunt_final_sub.csv')
df = df.rename(columns = {
    'x_new' : 'x',
    'y_new' : 'y',
    'z_new' : 'z'
})
# limit to only important columns
df = df[['name', 'age_myr', 'x', 'y', 'z']]



# ------------- Create Dash App -------------

app = dash.Dash(__name__)
app.layout = html.Div([
    html.Div([
        dcc.Graph(id='3d-scatter', style={'width': '50%', 'display': 'inline-block'}),
        dcc.Graph(id='2d-scatter', style={'width': '50%', 'display': 'inline-block'})
    ]),
    html.Div(
        dcc.RangeSlider(
            id='age-slider',
            min=df['age_myr'].min(),
            max=df['age_myr'].max(),
            value=[df['age_myr'].min(), df['age_myr'].max()],
            marks={str(age): str(age) for age in range(int(df['age_myr'].min()), int(df['age_myr'].max())+1, 5)},
            step=None
        ),
        style={'width': '50%', 'margin': 'auto'}
    )
])
@app.callback(
    Output('2d-scatter', 'figure'),
    Input('age-slider', 'value')
)
def update_2d_scatter(selected_age_range):
    filtered_df = df[(df['age_myr'] >= selected_age_range[0]) & (df['age_myr'] <= selected_age_range[1])]
    return create_2d_scatter(filtered_df)

@app.callback(
    Output('3d-scatter', 'figure'),
    Input('age-slider', 'value'),
    Input('2d-scatter', 'selectedData')
)
def update_3d_scatter(selected_age_range, selectedData):
    filtered_df = df[(df['age_myr'] >= selected_age_range[0]) & (df['age_myr'] <= selected_age_range[1])]
    return create_3d_scatter(filtered_df, selectedData)

def create_3d_scatter(df, selectedData=None):
    selected_clusters = df['name'].isin([point['hovertext'] for point in selectedData['points']]) if selectedData else False

    # Create a scatter plot for all clusters
    fig = go.Figure(data=[go.Scatter3d(
        x=df['x'],
        y=df['y'],
        z=df['z'],
        mode='markers',
        marker=dict(
            size=6,
            color='gray',  # color all clusters as gray
            opacity=0.8
        ),
        hovertext=df['name']
    )])

    # Add another scatter plot for the selected clusters
    if selectedData:
        fig.add_trace(go.Scatter3d(
            x=df[selected_clusters]['x'],
            y=df[selected_clusters]['y'],
            z=df[selected_clusters]['z'],
            mode='markers',
            marker=dict(
                size=6,
                color='red',  # color selected clusters as red
                opacity=0.8
            ),
            hovertext=df[selected_clusters]['name']
        ))

    # Set the x, y, z ranges and the camera
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-1000, 1000]),
            yaxis=dict(range=[-1000, 1000]),
            zaxis=dict(range=[-500, 500]),
            camera=dict(
                up=dict(x=0, y=1, z=0),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0., y=-.8, z=1.5),
                projection=dict(type='perspective')
            )
        ),
        dragmode='turntable',
        template='plotly_dark'
    )

    return fig
    
def create_2d_scatter(df):
    fig = go.Figure(data=go.Scatter(
        x=df['x'],
        y=df['y'],
        mode='markers',
        marker=dict(
            size=6,
            color=df['age_myr'],
            colorscale='Viridis',
            opacity=0.8
        ),
        hovertext=df['name']
    ))

    # Set the x, y ranges
    fig.update_xaxes(range=[-1000, 1000])
    fig.update_yaxes(range=[-1000, 1000])

    # Set the x, y ranges and the figure size
    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        xaxis=dict(range=[-1000, 1000]),
        yaxis=dict(range=[-1000, 1000])
    )

    return fig

if __name__ == '__main__':
    app.run_server(debug=True)



    