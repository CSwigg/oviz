import numpy as np
import plotly.graph_objects as go

class StarClusters3DPlotter:
    """
    Class for generating 3D plots of star clusters.

    Methods:
    - generate_3d_plot(self, collection): Generates a 3D plot for a StarClusterCollection.
    """
    def default_figure_layout(self):
        """
        Returns the default layout for the 3D plot figure.

        Returns:
        layout (plotly.graph_objs.Layout): The default layout for the 3D plot figure.
        """
        
        # plotly scene, layout, and camera
        camera = dict(
            up=dict(x=0, y=1, z=0),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=0., y=-.8, z=1.5),
            projection = dict(type = 'perspective')
        ) 

        x_width = 1500
        y_width = 1500
        z_width = 600
        z_aspect = 1.*(z_width/x_width)
        scene = dict(
            camera = camera,
            aspectmode = 'manual',
            aspectratio = dict(x=1, y=1, z=z_aspect),
            xaxis=dict(
                title='X\' (pc)', 
                range = [-x_width, x_width],
                showgrid=True,
                zeroline=False,
                showline = True,
                mirror = True,
                showspikes = False,
                nticks = 5,
                linecolor = 'white',
                linewidth = 2.,
            ), 
            yaxis=dict(
                title='Y\' (pc)', 
                range = [-y_width, y_width],
                showgrid=True,
                zeroline=False,
                showline = True,
                mirror = True,
                showspikes = False,
                nticks = 5,
                linecolor = 'white',
                linewidth = 2.,
            ), 
            zaxis=dict(
                title='Z\' (pc)', 
                range = [-z_width, z_width],
                showgrid=True,
                zeroline=False,
                showline = True,
                showspikes = False,
                mirror = True,
                nticks = 5,
                linecolor = 'white',
                linewidth = 2.,
            )
        )
        
        legend = dict(
            x = 0,
            y = 1,
            title = dict(
                text = 'Click to toggle traces on/off',
                font = dict(
                    size = 18,  
                    color = 'white'
                ),  
                side = 'top'
            ),
            font = dict(
                size = 14,
                family = 'helvetica',
                color = 'white'
            ),
            itemsizing = 'constant',
            bgcolor = 'rgba(0,0,0,0)' # transparent legend
        )
        layout = go.Layout(
            scene = scene,
            template = 'plotly_dark', 
            dragmode = 'turntable',
            width = 500,
            height = 500
        )
        #self.layout = layout
        return layout
    
    def generate_slider(self):
        """
        Generates a slider for controlling the animation of the plot.

        Returns:
        list: A list containing the slider configuration.
        """
        if self.figure_layout_dict['template'] == 'plotly_dark':
            slider_color = 'white'
        else:
            slider_color = 'black'

        sliders = [
            dict(
                active = len(self.time) - 1,
                xanchor ="left",
                yanchor="top",
                transition = {"duration": 300, "easing": "bounce-in"},
                borderwidth = 0.,
                bordercolor = slider_color,
                bgcolor = slider_color,
                pad = {"b": 0, "t": 0, "l": 0, "r": 0},
                len = 0.9,
                x = 0.1,
                y = 0.,  
                currentvalue = {
                    "font": {"size": 18, "color" : slider_color, 'family' : 'helvetica'},
                    #'prefix' : 'Time: ',
                    'visible' : True,
                    'xanchor':'center',
                    "offset" : 20},
                steps = [
                    dict(
                        args = [
                            [str(t)],
                            dict(
                                frame = dict(
                                    duration = 5,
                                    easing = 'linear',
                                    redraw = True
                                ),
                                transition = dict(
                                    duration = 0,
                                    easing = 'linear'
                                )
                            )
                        ],
                        label = 'Time: {} Myr'.format(t),
                        method = 'animate'
                    ) for t in np.flip(self.time)
                ],
                tickcolor = slider_color,
                ticklen = 10,
                font = dict(color = 'rgba(0,0,0,0)', size = 1)
            )
        ]
        return sliders
    
    def generate_3d_plot(self, collection, figure_layout=None, show=False, save_name=None, static_traces = None, static_traces_times = None, reference_frame_center = None):
        """
        Generates a 3D plot of star clusters over time.

        Parameters:
        - collection: StarClusterCollection object containing the star clusters to be plotted.
        - figure_layout: Optional parameter specifying the layout of the plot. If not provided, a default layout will be used.
        - show: Boolean indicating whether to display the plot.
        - save: Boolean indicating whether to save the plot as an HTML file.
        - save_name: Optional parameter specifying the name of the HTML file to save.

        Returns:
        None
        """
        
        if collection.time is None: # handles if star clusters haven't been integrated yet
            print('Star clusters have not yet been integrated \n Integrating backwards 60 Myr with 1 Myr time steps...')
            collection.integrate_clusters(np.arange(0,-60,-1))
        collection.set_all_cluster_sizes()
        self.time = collection.time 
        if figure_layout is None:
            self.figure_layout = self.default_figure_layout()
        else:
            self.figure_layout = go.Layout(figure_layout)
            self.figure_layout_dict = figure_layout
        
        '''
        The following nested loops are strucutred such that for each step in time,
        every star cluster across the entire collection is plotted.
        '''
        cluster_groups = collection.get_all_clusters()
        fig = {}
        frames = []
        for t in self.time:
            scatter_list = []
            for cluster_group in cluster_groups:
                assert cluster_group.integrated == True # make sure that the integrated DataFrame is populated
                frame = {'data': [], 'name': str(t)}
                df_int = cluster_group.df_int
                if len(df_int) == 0:
                    continue
                
                df_t = df_int[df_int['time'] == t]
                scatter_cluster_group_t = go.Scatter3d(
                    x=df_t['x'],
                    y=df_t['y'],
                    z=df_t['z'],
                    mode='markers',
                    marker = dict(
                        size = df_t['size'],
                        color = cluster_group.color,
                        opacity = cluster_group.opacity,
                        line = dict(
                            color = 'black',
                            width = 0.0 if self.figure_layout_dict['template'] == 'plotly_dark' else 10.
                        )
                    ),
                    hovertext = df_t['name'],
                    name=cluster_group.data_name
                ) 
                scatter_list.append(scatter_cluster_group_t)
            frame['data'] = scatter_list 

            if static_traces is not None and static_traces_times is not None:
                for i, st in enumerate(static_traces):
                    if reference_frame_center is not None:
                        st['x'] = st['x'] - reference_frame_center[0]
                        st['y'] = st['y'] - reference_frame_center[1]
                        st['z'] = st['z'] - reference_frame_center[2]
                    if t in static_traces_times[i]:
                        visible = True
                    else:
                        visible = False
                        # place holder
                        st = go.Scatter3d()
                    st['visible'] = visible
                    frame['data'].append(st)
            frames.append(go.Frame(frame))
            
        fig['data'] = frames[0]['data']
        fig['layout'] = self.figure_layout
        # Generate slider
        slider = self.generate_slider()
        fig['layout']['sliders'] = slider
        fig['frames'] = frames
        self.fig_dict = fig
        self.figure = go.Figure(fig)
        if show:
            self.figure.show()
        if save_name is not None:
            self.figure.update_layout(width = None, height = None)
            self.figure.write_html(save_name, auto_play = False)