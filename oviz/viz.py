import numpy as np
import yaml
import plotly.graph_objects as go
from astropy.coordinates import SkyCoord
import astropy.units as u
from . import orbit_maker
import copy

class Animate3D:
    """
    Class for generating 3D plots of star clusters.

    Methods:
    - generate_3d_plot(self, collection): Generates a 3D plot for a StarClusterCollection.
    """

    def __init__(self, data_collection, xyz_widths=(1000, 1000, 300), xyz_ranges=None,  figure_title=None, figure_theme=None, plotly_light_template=None, figure_layout=None, figure_layout_dict=None, trace_grouping_dict=None):
        """
        Initializes a StarClusters3DPlotter object.

        Parameters:
        - data_collection (StarClusterCollection): The collection of star clusters to be plotted.
        - xyz_widths (tuple): Widths of the plot in x, y, and z directions.
        - xyz_ranges (tuple): Ranges of the plot in x, y, and z directions.
        - figure_title (str): The title of the plot.
        - figure_theme (str): The theme of the plot. Options are 'light' and 'dark'.
        - plotly_light_template (str): The Plotly template for light theme.
        - figure_layout (plotly.graph_objs.Layout): The layout of the plot.
        - figure_layout_dict (dict): The layout of the plot as a dictionary.
        - trace_grouping_dict (dict): Dictionary for grouping traces.
        """
        self.data_collection = data_collection
        self.figure_theme = figure_theme
        self.figure_layout = figure_layout
        self.figure_layout_dict = figure_layout_dict
        self.time = None
        self.fig_dict = None
        self.figure = None
        self.xyz_widths = xyz_widths
        self.xyz_ranges = xyz_ranges
        self.plot_light_template = plotly_light_template
        self.trace_grouping_dict = trace_grouping_dict or {}
        self.figure_title = figure_title

        self.figure_layout_dict = read_theme(self)

        if self.figure_title:
            if self.figure_theme == 'dark':
                self.font_color = 'white'
            else:
                self.font_color = 'black'
            self.figure_layout_dict['title'] = dict(text=self.figure_title, x=0.5, font=dict(family='Helvetica', size=20, color=self.font_color))

        self.trace_grouping_dict['All'] = [cluster.data_name for cluster in self.data_collection.get_all_clusters()]


    def rotating_gc_line(self, x_sub, y_sub):
        """
        Generates a rotating galactic center line.

        Parameters:
        - x_sub (float): X-coordinate subtraction value.
        - y_sub (float): Y-coordinate subtraction value.

        Returns:
        go.Scatter3d: A Plotly Scatter3d object representing the galactic center line.
        """
        n_marks = 1000
        R = -8.122 * np.ones(n_marks)
        phi = np.linspace(-180, 180, n_marks)
        z = np.zeros(n_marks)
        gc_line = SkyCoord(
            rho=R * u.kpc,
            phi=phi * u.deg,
            z=[0.]*len(R) * u.pc,
            frame='galactocentric',
            representation_type='cylindrical'
        )
        x_corr = gc_line.galactic.cartesian.x.value * 1000 - x_sub
        y_corr = gc_line.galactic.cartesian.y.value * 1000 - y_sub
        z_corr = [0.]*len(x_corr)
        
        line_color = 'gray' if self.figure_theme == 'dark' else 'black'

        scatter_gc_line = go.Scatter3d(
            x=x_corr,
            y=y_corr,
            z=z_corr,
            mode='lines',
            line=dict(color=line_color, width=3., dash='solid'),
            visible=True,
            name='R = 8.12 kpc',
            hovertext='R = 8.12 kpc'
        )
        return scatter_gc_line

    def generate_slider(self):
        """
        Generates a slider for controlling the animation of the plot.

        Returns:
        list: A list containing the slider configuration.
        """
        slider_color = 'gray' if self.figure_theme == 'dark' else 'black'

        sliders = [
            dict(
                active=len(self.time) - 1,
                xanchor="left",
                yanchor="top",
                transition={"duration": 300, "easing": "bounce-in"},
                borderwidth=0.,
                bordercolor=slider_color,
                bgcolor=slider_color,
                pad={"b": 0, "t": 0, "l": 0, "r": 0},
                len=0.5,
                x=0.27,
                y=0.,  
                currentvalue={
                    "font": {"size": 18, "color": slider_color, 'family': 'helvetica'},
                    'prefix': 'Time (Myr ago): ',
                    'visible': True,
                    'xanchor': 'center',
                    "offset": 20
                },
                steps=[
                    dict(
                        args=[
                            [str(t)],
                            dict(
                                frame=dict(
                                    duration=5,
                                    easing='linear',
                                    redraw=True
                                ),
                                transition=dict(
                                    duration=0,
                                    easing='linear'
                                )
                            )
                        ],
                        label=-1*t if t != 0 else 'Present Day',
                        method='animate'
                    ) for t in np.flip(self.time)
                ],
                tickcolor=slider_color,
                ticklen=10,
                font=dict(color='rgba(0,0,0,0)', size=8, family='helvetica')
            )
        ]
        return sliders

    def set_focus(self, focus_group):
        """
        Sets the focus on a specific group of star clusters.

        Parameters:
        - focus_group (str): The name of the focus group.

        Returns:
        np.ndarray: The mean coordinates of the focus group.
        """
        if focus_group:
            focus_group_data = self.data_collection.get_cluster(focus_group).df
            focus_group_coords = focus_group_data[['x', 'y', 'z', 'U', 'V', 'W']].median().values
            return focus_group_coords
        else:
            return None

    def get_visibility(self, traces_list):
        """
        Determines the visibility of traces in the plot.

        Parameters:
        - traces_list (list): List of trace names to be visible.

        Returns:
        list: A list of boolean values indicating the visibility of each trace.
        """
        visibility = []
        for trace in self.fig['data']:
            if trace['name'] in traces_list:
                visibility.append(True)
            elif trace['name'] in self.static_trace_names:
                visibility.append(False)
            elif trace['name'] == 'R = 8.12 kpc':
                visibility.append(True)
            else:
                visibility.append(False)
        return visibility

    def dropdown_menu(self):
        """
        Creates a dropdown menu for selecting trace groups in a 3D plot.

        Returns:
        list: A list containing a dictionary that defines the dropdown menu configuration for Plotly.
        """
        buttons = []
        for key, traces_list in self.trace_grouping_dict.items():
            visibility = self.get_visibility(traces_list)
            buttons.append(dict(
                args=[{'visible': visibility, 'frame': 0.}],
                label=key,
                method='restyle'
            ))
        bg_color = self.figure_layout_dict['scene']['xaxis']['backgroundcolor']
        text_color = 'black' if self.figure_theme != 'dark' else 'white'
        dropdown = [dict(
            buttons=buttons,
            direction='down',
            pad={'r': 10, 't': 10},
            showactive=True,
            x=0,
            xanchor='left',
            y=1.1,
            bgcolor=bg_color,
            font=dict(color=text_color, family='helvetica', size=14),
            yanchor='top',
            active=0
        )]
        return dropdown

    def make_plot(self, time=None, figure_layout=None, show=False, save_name=None, static_traces=None, static_traces_times=None, static_traces_legendonly=False, reference_frame_center=None, focus_group=None, fade_in_time=5, fade_in_and_out=False, show_gc_line=True):
        """
        Generates a 3D plot of star clusters over time.

        Parameters:
        - time (list): List of time points for the animation.
        - figure_layout (plotly.graph_objs.Layout): Optional layout for the plot.
        - show (bool): Whether to display the plot.
        - save_name (str): Optional name for saving the plot as an HTML file.
        - static_traces (list): List of static traces to be included in the plot.
        - static_traces_times (list): List of times for the static traces.
        - static_traces_legendonly (bool): Whether to show static traces only in the legend.
        - reference_frame_center (tuple): Center of the reference frame.
        - focus_group (str): Name of the focus group.
        - fade_in_time (int): Time for fading in the traces.
        - fade_in_and_out (bool): Whether to fade in and out the traces.
        - show_gc_line (bool): Whether to show the galactic center line.

        Returns:
        go.Figure: The generated 3D plot.
        """
        if reference_frame_center is None:
            reference_frame_center = self.set_focus(focus_group)
        if self.data_collection.time is None:
            self.data_collection.integrate_all_orbits(time, reference_frame_center=reference_frame_center)
        self.data_collection.set_all_cluster_sizes(fade_in_time, fade_in_and_out)

        self.time = np.array(self.data_collection.time, dtype=np.float64)
        self.figure_layout = go.Layout(self.figure_layout_dict)
        x_rf_int, y_rf_int, z_rf_int = orbit_maker.get_center_orbit_coords(time, reference_frame_center)
        self.coords_center_int = (x_rf_int, y_rf_int, z_rf_int)
        self.focus_group = focus_group

        cluster_groups = self.data_collection.get_all_clusters()
        self.fig = {}
        frames = []

        for i, t in enumerate(self.time):
            scatter_list = self._generate_scatter_list(cluster_groups, t, x_rf_int[i], y_rf_int[i], show_gc_line)
            frame = {'data': scatter_list, 'name': str(t)}
            self._add_static_traces(frame, static_traces, static_traces_times, reference_frame_center, t, static_traces_legendonly)
            frames.append(go.Frame(frame))

        self._initialize_figure(frames, static_traces_legendonly)
        self._add_slider_and_dropdown()

        if show:
            self.figure.show()
        if save_name:
            self.figure.update_layout(width=None, height=None)
            self.figure.write_html(save_name, auto_play=False)
        return self.figure

    def _generate_scatter_list(self, cluster_groups, t, x_rf, y_rf, show_gc_line):
        scatter_list = []
        for cluster_group in cluster_groups:
            assert cluster_group.integrated
            df_int = cluster_group.df_int
            if len(df_int) == 0:
                continue

            df_t = df_int[df_int['time'] == t]
            hovertext = df_t['name'].astype(str) + ': ' + df_t['age_myr'].round(1).astype(str) + ' Myr'
            if 'n_stars' in df_t.columns:
                hovertext += ', N = ' + df_t['n_stars'].astype(str)
                
            #df_t = orbit_maker.coordFIX_to_coordROT(df_t)
            x = df_t['x'].values
            y = df_t['y'].values
            z = df_t['z'].values
            scatter_cluster_group_t = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(
                    size=df_t['size'],
                    color=cluster_group.color,
                    opacity=cluster_group.opacity,
                    symbol=cluster_group.marker_style,
                    line=dict(color='black', width=0.0)
                ),
                hovertext=hovertext,
                name=cluster_group.data_name
            )
            scatter_list.append(scatter_cluster_group_t)

        if show_gc_line:
            gc_line_t = self.rotating_gc_line(x_rf, y_rf)
            scatter_list.append(gc_line_t)
        return scatter_list

    def _add_static_traces(self, frame, static_traces, static_traces_times, reference_frame_center, t, static_traces_legendonly):
        self.static_trace_names = []

        if static_traces is None:
            static_traces = []
            static_traces_times = []

        for sc in self.data_collection.get_all_clusters():
            if sc.show_tracks:
                track = plot_trace_tracks(sc)
                static_traces.append(track)
                static_traces_times.append([0])

        if static_traces and static_traces_times:
            for i, st in enumerate(static_traces):
                st_copy = copy.deepcopy(st)
                self.static_trace_names.append(st_copy['name'])
                if self.focus_group is not None and not st_copy['name'].endswith('Track'):
                    if st_copy['x'] is not None:
                        st_copy['x'] = np.array(st_copy['x']) - reference_frame_center[0]
                    if st_copy['y'] is not None:
                        st_copy['y'] = np.array(st_copy['y']) - reference_frame_center[1]
                    if st_copy['z'] is not None:
                        st_copy['z'] = np.array(st_copy['z']) - reference_frame_center[2]
                if t == 0 and not static_traces_legendonly:
                    visible = True
                elif t == 0 and static_traces_legendonly:
                    visible = 'legendonly'
                elif np.any(np.array(static_traces_times[i]) == t):
                    visible = True
                else:
                    visible = False
                    st_copy = go.Scatter3d()
                st_copy['visible'] = visible
                frame['data'].append(st_copy)

    def _initialize_figure(self, frames, static_traces_legendonly):
        grouping_0 = list(self.trace_grouping_dict.values())[0]
        starting_frame = copy.deepcopy(frames[0])
        data_updated = []
        for trace in starting_frame['data']:
            trace_name = trace['name']
            visible = True if trace_name in grouping_0 or trace_name in self.static_trace_names else False
            visible = 'legendonly' if trace_name in self.static_trace_names and static_traces_legendonly else visible
            visible = True if trace_name == 'R = 8.12 kpc' else visible
            trace['visible'] = visible
            data_updated.append(trace)

        starting_frame['data'] = data_updated

        self.fig['data'] = starting_frame['data']
        self.fig['layout'] = self.figure_layout
        self.fig['frames'] = frames

    def _add_slider_and_dropdown(self):
        slider = self.generate_slider()
        self.fig['layout']['sliders'] = slider

        if len(self.trace_grouping_dict) > 1:
            dropdown = self.dropdown_menu()
            self.fig['layout']['updatemenus'] = dropdown

        self.fig_dict = self.fig
        self.figure = go.Figure(self.fig)

def plot_trace_tracks(sc):
    """
    Plots the tracks of a star cluster.

    Parameters:
    - sc (StarCluster): The star cluster object.

    Returns:
    go.Scatter3d: A Plotly Scatter3d object representing the tracks of the star cluster.
    """
    df_int = sc.df_int
    df_int = df_int.loc[df_int['time'] > -1*df_int['age_myr']]
    tracks = go.Scatter3d(
        x=df_int['x'].iloc[::1],
        y=df_int['y'].iloc[::1],
        z=df_int['z'].iloc[::1],
        mode='markers',
        marker=dict(
            size=(sc.max_size/3)*(1 - np.abs(df_int['time'])/np.abs(df_int['time']).max()),
            #size=(sc.max_size/3)*(1 - np.abs(df_int['time'])/df_int['age_myr']),
            color=sc.color,
            opacity=sc.opacity/1.5,
            line=dict(width=0)
        ),
        name=sc.data_name + ' Track'
    )
    return tracks

def read_theme(plot):
    """
    Reads the theme configuration from a YAML file.

    Parameters:
    - plot (StarClusters3DPlotter): The plotter object.

    Returns:
    dict: The layout configuration.
    """
    theme_dir = './themes'
    with open(theme_dir + '/{}.yaml'.format(plot.figure_theme), 'r') as file:
        layout = yaml.safe_load(file)

    if plot.figure_theme == 'light':
        layout['template'] = plot.plot_light_template if plot.plot_light_template else 'ggplot2'

    if plot.xyz_ranges:
        x_ranges, y_ranges, z_ranges = plot.xyz_ranges
        x_low, x_high = x_ranges
        y_low, y_high = y_ranges
        z_low, z_high = z_ranges
        x_width = x_high - x_low
        y_width = y_high - y_low
        z_width = z_high - z_low
        xy_aspect = y_width / x_width
        z_aspect = z_width / x_width
    else:
        x_width, y_width, z_width = plot.xyz_widths
        x_low, x_high, y_low, y_high, z_low, z_high = -x_width, x_width, -y_width, y_width, -z_width, z_width
        xy_aspect = y_width / x_width
        z_aspect = z_width / x_width

    layout['scene']['xaxis']['range'] = [x_low, x_high]
    layout['scene']['yaxis']['range'] = [y_low, y_high]
    layout['scene']['zaxis']['range'] = [z_low, z_high]
    layout['scene']['aspectratio']['x'] = 1
    layout['scene']['aspectratio']['y'] = xy_aspect
    layout['scene']['aspectratio']['z'] = z_aspect

    return layout
