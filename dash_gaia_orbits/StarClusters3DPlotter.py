import numpy as np
import yaml
import plotly.graph_objects as go

class StarClusters3DPlotter:
    """
    Class for generating 3D plots of star clusters.

    Methods:
    - generate_3d_plot(self, collection): Generates a 3D plot for a StarClusterCollection.
    """


    def __init__(self, data_collection, xyz_widths = (1000, 1000, 300), figure_theme = None, plotly_light_template = None, figure_layout = None, figure_layout_dict = None):
        """
        Initializes a StarClusters3DPlotter object.

        Parameters:
        - data_collection (StarClusterCollection): The collection of star clusters to be plotted.
        - figure_theme (str): The theme of the plot. Options are 'light' and 'dark'.
        - figure_layout (plotly.graph_objs.Layout): The layout of the plot.
        - figure_layout_dict (dict): The layout of the plot as a dictionary.
        """
        self.data_collection = data_collection
        self.figure_theme = figure_theme
        self.figure_layout = figure_layout
        self.figure_layout_dict = figure_layout_dict
        self.time = None
        self.fig_dict = None
        self.figure = None
        self.xyz_widths = xyz_widths
        self.figure_theme = figure_theme
        self.plot_light_template = plotly_light_template

        self.figure_layout_dict = read_theme(self)


    
    def generate_slider(self):
        """
        Generates a slider for controlling the animation of the plot.

        Returns:
        list: A list containing the slider configuration.
        """
        if self.figure_theme == 'dark':
            slider_color = 'gray'
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
                len = 0.5,
                x = 0.27,
                y = 0.,  
                currentvalue = {
                    "font": {"size": 18, "color" : slider_color, 'family' : 'helvetica'},
                    'prefix' : 'Time (Myr ago): ',
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
                        #label = 'Time: {} Myr'.format(t),
                        label = -1*t,
                        method = 'animate'
                    ) for t in np.flip(self.time)
                ],
                tickcolor = slider_color,
                ticklen = 10,
                font = dict(color = 'rgba(0,0,0,0)', size = 8, family = 'helvetica')
                #font = dict(color = 'rgba(0,0,0,0)', size = 1)
            )
        ]
        return sliders
    
    
    def set_focus(self, focus_group):
        if focus_group is not None:
            focus_group_data = self.data_collection.get_cluster(focus_group).df
            focus_group_coords = focus_group_data[['x', 'y', 'z', 'U', 'V', 'W']].mean().values
            return focus_group_coords
        else:
            return None





    def generate_3d_plot(self, time = None, figure_layout=None, show=False, save_name=None, static_traces = None, static_traces_times = None, reference_frame_center = None, focus_group = None, fade_in_time = 5):
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
        
        if reference_frame_center is None:
            reference_frame_center = self.set_focus(focus_group)
        if self.data_collection.time is None: # handles if star clusters haven't been integrated yet
            self.data_collection.integrate_all_orbits(time, reference_frame_center = reference_frame_center)
        self.data_collection.set_all_cluster_sizes(fade_in_time)

        self.time = self.data_collection.time 
        self.figure_layout = go.Layout(self.figure_layout_dict)

        

        '''
        The following nested loops are strucutred such that for each step in time,
        every star cluster across the entire collection is plotted.
        '''
        cluster_groups = self.data_collection.get_all_clusters()
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
                        symbol = cluster_group.marker_style,
                        line = dict(
                            color = 'black',
                            width = 0.0
                            #width = 0.0 if self.figure_theme == 'dark' else 3.
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





def read_theme(plot : StarClusters3DPlotter):

    with open('../themes/{}.yaml'.format(plot.figure_theme), 'r') as file:
        layout = yaml.safe_load(file)

    if plot.figure_theme == 'light':
            layout['template'] = plot.plot_light_template if plot.plot_light_template is not None else 'ggplot2'

    x_width, y_width, z_width = plot.xyz_widths
    xy_aspect = 1.*(y_width/x_width)
    z_aspect = 1.*(z_width/x_width)

    layout['scene']['xaxis']['range'] = [-x_width, x_width]
    layout['scene']['yaxis']['range'] = [-y_width, y_width]
    layout['scene']['zaxis']['range'] = [-z_width, z_width]
    layout['scene']['aspectratio']['x'] = 1
    layout['scene']['aspectratio']['y'] = xy_aspect
    layout['scene']['aspectratio']['z'] = z_aspect


    return layout
