import numpy as np
import yaml
import plotly.graph_objects as go
from astropy.coordinates import SkyCoord
import astropy.units as u
import importlib.resources
from . import orbit_maker
import copy

class Animate3D:
    """
    Class for generating 3D plots of star clusters.
    """

    def __init__(
        self,
        data_collection,
        xyz_widths=(1000, 1000, 300),
        xyz_ranges=None,
        figure_title=None,
        figure_theme=None,
        plotly_light_template=None,
        figure_layout=None,
        figure_layout_dict=None,
        trace_grouping_dict=None
    ):
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

        # Read in the layout from the theme
        self.figure_layout_dict = read_theme(self)

        # Configure the main figure title if provided
        if self.figure_title:
            self.font_color = 'white' if self.figure_theme == 'dark' else 'black'
            self.figure_layout_dict['title'] = dict(
                text=self.figure_title,
                x=0.5,
                font=dict(family='Helvetica', size=20, color=self.font_color)
            )

        # Always include a default grouping "All"
        self.trace_grouping_dict['All'] = [
            cluster.data_name for cluster in self.data_collection.get_all_clusters()
        ]

    def make_plot(
        self,
        time=None,
        figure_layout=None,
        show=False,
        save_name=None,
        static_traces=None,
        static_traces_times=None,
        static_traces_legendonly=False,
        reference_frame_center=None,
        focus_group=None,
        fade_in_time=5,
        fade_in_and_out=False,
        show_gc_line=True
    ):
        """
        Main method to generate a 3D animated plot. Integrates orbits if necessary and
        constructs frames for each time step. Optionally saves the figure to HTML.
        """
        # Ensure t=0 is present
        if 0 not in time:
            raise ValueError("The time array must include 0 for the present day.")

        # Determine the reference frame center
        if reference_frame_center is None:
            reference_frame_center = self.set_focus(focus_group)

        # Integrate orbits if not already integrated
        if self.data_collection.time is None:
            self.data_collection.integrate_all_orbits(
                time, reference_frame_center=reference_frame_center
            )
        # Set cluster sizes for fade effects
        self.fade_in_time = fade_in_time
        self.data_collection.set_all_cluster_sizes(self.fade_in_time, fade_in_and_out)

        # Prepare time arrays and figure layout
        self.time = np.array(self.data_collection.time, dtype=np.float64)
        self.figure_layout = go.Layout(self.figure_layout_dict)
        self.focus_group = focus_group
        self.static_traces_legendonly = static_traces_legendonly

        x_rf_int, y_rf_int, z_rf_int = orbit_maker.get_center_orbit_coords(
            time, reference_frame_center
        )
        self.coords_center_int = (x_rf_int, y_rf_int, z_rf_int)

        # Initialize static traces if none given
        if static_traces is None:
            static_traces = []
        if static_traces_times is None:
            static_traces_times = []

        # Record static traces – for example, add orbit tracks as static traces
        for sc in self.data_collection.get_all_clusters():
            if sc.show_tracks:
                track = plot_trace_tracks(sc, self.fade_in_time)
                static_traces.append(track)
                static_traces_times.append([0])  # Show only at t=0

        # Save the names of all static traces (if they have a name)
        self.base_static_trace_names = [
            st.name for st in static_traces if hasattr(st, 'name') and st.name
        ]

        # Build frames for each time step
        cluster_groups = self.data_collection.get_all_clusters()
        self.fig = {}
        frames = []

        for i, t_i in enumerate(self.time):
            # Generate scatter traces for each cluster at time t_i
            scatter_list = self._generate_scatter_list(
                cluster_groups, t_i, x_rf_int[i], y_rf_int[i], show_gc_line
            )
            # Remove any 'visible' property so frames won't override grouping
            for sc_tr in scatter_list:
                sc_tr.pop('visible', None)

            # Create the frame dict
            frame_data = {'data': scatter_list, 'name': str(t_i)}

            # Add static traces if they should appear at this time.
            self._add_static_traces(
                frame_data, static_traces, static_traces_times,
                reference_frame_center, t_i
            )

            # Remove 'visible' from static traces too, for the same reason
            for item in frame_data['data']:
                item.pop('visible', None)

            frames.append(go.Frame(frame_data))

        # Initialize the figure at t=0 (base data)
        self._initialize_figure(frames)

        # Add slider & dropdown, and finalize
        self._add_slider_and_dropdown()

        # Show or save
        if show:
            self.figure.show()
        if save_name:
            self.figure.update_layout(width=None, height=None)
            self.figure.write_html(save_name, auto_play=False)

        return self.figure

    def rotating_gc_line(self, x_sub, y_sub):
        """
        Generates a rotating galactic center line, displayed as a 3D line in the figure.
        """
        n_marks = 1000
        R = -8.122 * np.ones(n_marks)
        phi = np.linspace(-180, 180, n_marks)
        z = np.zeros(n_marks)
        gc_line = SkyCoord(
            rho=R * u.kpc,
            phi=phi * u.deg,
            z=[0.] * len(R) * u.pc,
            frame='galactocentric',
            representation_type='cylindrical'
        )

        x_corr = gc_line.galactic.cartesian.x.value * 1000 - x_sub
        y_corr = gc_line.galactic.cartesian.y.value * 1000 - y_sub
        z_corr = [0.] * len(x_corr)

        line_color = 'gray' if self.figure_theme == 'dark' else 'black'
        return go.Scatter3d(
            x=x_corr,
            y=y_corr,
            z=z_corr,
            mode='lines',
            line=dict(color=line_color, width=3., dash='solid'),
            visible=True,
            name='R = 8.12 kpc',
            hovertext='R = 8.12 kpc'
        )

    def generate_play_pause(self):
        """
        Generates play/pause buttons for the plot, used to control frame animation.
        """
        if self.figure_theme == 'dark':
            button_color, text_color = 'gray', 'gray'
        elif self.figure_theme == 'light':
            button_color, text_color = 'black', 'black'
        else:  # 'gray'
            button_color, text_color = 'black', 'black'

        button_color = 'gray' if self.figure_theme == 'dark' else 'black'
        text_color = 'black'

        return [
            dict(
                type='buttons',
                showactive=False,
                x=0.2,
                y=-0.03,
                xanchor='left',
                yanchor='top',
                direction='left',
                pad={'r': 50, 't': 20, 'b': 20, 'l': 20},
                bgcolor='rgba(0,0,0,0)',
                bordercolor=button_color,
                font=dict(color=text_color, size=20, family='helvetica'),
                buttons=[
                    dict(
                        label='▶',
                        method='animate',
                        args=[
                            None,
                            dict(frame=dict(duration=500, redraw=True), fromcurrent=True)
                        ]
                    ),
                    dict(
                        label='⏸',
                        method='animate',
                        args=[[None], dict(frame=dict(duration=0, redraw=True), mode='immediate')]
                    )
                ]
            )
        ]

    def generate_slider(self):
        """
        Generates a slider to control the animation across different time steps.
        """
        slider_color = 'gray' if self.figure_theme == 'dark' else 'black'

        time_neg = self.time[self.time < 0]
        time_pos = self.time[self.time >= 0]

        # Sort out an order in which to present the slider times
        if (len(time_neg) > 0) and (len(time_pos) > 1):
            time_slider = np.append(time_neg, time_pos)
        elif (len(time_neg) > 0) and (len(time_pos) == 1):
            time_slider = np.flip(self.time)
        else:
            time_slider = self.time

        # Locate the index where time is zero
        zero_idx = np.where(time_slider == 0)[0][0]

        return [
            dict(
                active=zero_idx,
                xanchor="center",
                yanchor="top",
                transition={"duration": 300, "easing": "bounce-in"},
                borderwidth=0.0,
                bordercolor=slider_color,
                bgcolor=slider_color,
                pad={"b": 0, "t": 0, "l": 0, "r": 0},
                len=0.5,     # Shrink slider to half its default length
                x=0.5,       # Center slider on screen
                y=0.0,       # Keep it at the bottom
                currentvalue={
                    "font": {"size": 18, "color": slider_color, "family": "helvetica"},
                    "prefix": "Time (Myr): ",
                    "visible": True,
                    "xanchor": "center",
                    "offset": 20
                },
                steps=[
                    dict(
                        args=[
                            [str(t)],
                            dict(
                                frame=dict(duration=5, easing="linear", redraw=True),
                                transition=dict(duration=0, easing="linear")
                            )
                        ],
                        label=str(t),
                        method="animate"
                    ) for t in time_slider
                ],
                tickcolor=slider_color,
                ticklen=10,
                font=dict(color="rgba(0,0,0,0)", size=8, family="helvetica")
            )
        ]

    def set_focus(self, focus_group):
        """
        Returns median coordinates of a focus group if provided; otherwise returns None.
        """
        if not focus_group:
            return None

        focus_group_data = self.data_collection.get_cluster(focus_group).df
        coords = focus_group_data[['x', 'y', 'z', 'U', 'V', 'W']].median().values
        return coords

    def get_visibility(self, trace_name: str, grouping: list):
        """
        Determines if a given trace (by name) should be visible under the current grouping.
        
        - For static non‑track traces (i.e. those in self.base_static_trace_names that do not end with ' Track'):
          if static_traces_legendonly is True, return "legendonly", otherwise return True.
        - For static track traces (names ending with ' Track'):
          check if the corresponding base name is in grouping; if so, return "legendonly" when
          static_traces_legendonly is True, otherwise return True; if not, return False.
        - The galactic center line always returns True.
        - For non‑static traces, return True if the trace name is in grouping.
        """
        # For static non-track traces:
        if trace_name in self.base_static_trace_names and not trace_name.endswith(" Track"):
            if self.static_traces_legendonly:
                return "legendonly"
            else:
                return True

        # For static track traces:
        if trace_name.endswith(" Track"):
            base_name = trace_name.replace(" Track", "")
            if base_name in grouping:
                if self.static_traces_legendonly:
                    return "legendonly"
                else:
                    return True
            else:
                return False

        # Special case for the galactic center line.
        if trace_name == 'R = 8.12 kpc':
            return True
        # For non-static traces:
        elif trace_name in grouping:
            return True

        return False

    def _generate_scatter_list(self, cluster_groups, t, x_rf, y_rf, show_gc_line):
        """
        Generate the main cluster Scatter3d traces for a given time.
        Optionally includes the rotating galactic center line if show_gc_line is True.
        Non-static traces need not be “marked” since the dropdown update will control them.
        """
        scatter_list = []
        for cluster_group in cluster_groups:
            assert cluster_group.integrated
            df_int = cluster_group.df_int
            if df_int.empty:
                continue

            df_t = df_int[df_int['time'] == t]
            hovertext = (
                '<b style="font-size:16px;">' + df_t['name'].str.replace('_', ' ').astype(str) + '</b>' + '<br>'  # Bold cluster name
                + cluster_group.data_name + '<br>'  # Group name
                + 'Age = ' + df_t['age_myr'].round(1).astype(str) + ' Myr' + '<br>'  # Cluster age
            )

            if 'n_stars' in df_t.columns:
                hovertext += 'N = ' + df_t['n_stars'].astype(str) + ' stars <br>'  # Number of stars

            hovertext += (
                '(x,y,z) = (' +
                df_t['x'].round(1).astype(str) + ', ' +
                df_t['y'].round(1).astype(str) + ', ' +
                df_t['z'].round(1).astype(str) + ')'
            )

            scatter_list.append(
                go.Scatter3d(
                    x=df_t['x'].values,
                    y=df_t['y'].values,
                    z=df_t['z'].values,
                    mode='markers',
                    marker=dict(
                        size=df_t['size'],
                        color=cluster_group.color,
                        opacity=cluster_group.opacity,
                        symbol=cluster_group.marker_style,
                        line=dict(color='black', width=0.0)
                    ),
                    hovertext=hovertext,
                    hoverinfo='text',  # This removes default x, y, z
                    hovertemplate='%{hovertext}<extra></extra>',  # This ensures only custom hovertext is shown
                    name=cluster_group.data_name
                )
            )

        if show_gc_line:
            gc_line_t = self.rotating_gc_line(x_rf, y_rf)
            scatter_list.append(gc_line_t)

        return scatter_list

    def _add_static_traces(self, frame, static_traces, static_traces_times, reference_frame_center, t):
        """
        Adds static traces to the frame if the current time t is in static_traces_times[i].
        If not, an empty Scatter3d is appended to override any previously shown trace.
        Also re-centers if a focus group is set.
        Here we “mark” these traces as static by adding a meta field.
        """
        for i, st in enumerate(static_traces):
            st_copy = copy.deepcopy(st)
            st_copy['meta'] = {'static': True}

            # Re-center if focusing on a group (except for tracks)
            if (self.focus_group is not None) and (st_copy.name) and not st_copy.name.endswith('Track'):
                if st_copy.x is not None:
                    st_copy.x = np.array(st_copy.x) - reference_frame_center[0]
                if st_copy.y is not None:
                    st_copy.y = np.array(st_copy.y) - reference_frame_center[1]
                if st_copy.z is not None:
                    st_copy.z = np.array(st_copy.z) - reference_frame_center[2]

            if t in static_traces_times[i]:
                frame['data'].append(st_copy)
            else:
                frame['data'].append(go.Scatter3d(
                    x=[], y=[], z=[],
                    name=st_copy.name,
                    visible=False,
                    meta={'static': True}
                ))

    def _initialize_figure(self, frames):
        """
        Sets up the base figure at t=0 using the default grouping.
        """
        default_group_key = list(self.trace_grouping_dict.keys())[0]  # e.g. "All"
        grouping_0 = self.trace_grouping_dict[default_group_key]

        # Find the frame for t=0
        idx_zero = np.where(self.time == 0)[0][0]
        starting_frame = copy.deepcopy(frames[idx_zero])

        data_updated = []
        for trace in starting_frame['data']:
            visible_flag = self.get_visibility(trace['name'], grouping_0)
            trace['visible'] = visible_flag
            data_updated.append(trace)

        self.initial_data = copy.deepcopy(data_updated)
        self.fig = {
            'data': data_updated,
            'layout': self.figure_layout,
            'frames': frames
        }

    def _add_slider_and_dropdown(self):
        """
        Adds the slider, play/pause buttons, and (if applicable) the dropdown menu to the layout.
        """
        slider = self.generate_slider()
        play_pause = self.generate_play_pause()
        self.fig['layout']['sliders'] = slider

        self.figure = go.Figure(self.fig)

        if len(self.trace_grouping_dict) > 1:
            dropdown = self.dropdown_menu()
            self.figure.update_layout(updatemenus=dropdown)
            self.fig['layout']['updatemenus'] = dropdown

        self.fig_dict = self.fig



####################################################################################################
def read_theme(plot):
    """
    Reads the theme configuration from a YAML file and sets up figure layout based on
    the provided figure_theme (e.g., 'light', 'dark', etc.).
    """
    with importlib.resources.open_text("oviz.themes", f"{plot.figure_theme}.yaml") as file:
        layout = yaml.safe_load(file)

    if plot.figure_theme == 'light':
        layout['template'] = (
            plot.plot_light_template if plot.plot_light_template else 'ggplot2'
        )

    if plot.xyz_ranges:
        (x_low, x_high), (y_low, y_high), (z_low, z_high) = plot.xyz_ranges
        x_width = x_high - x_low
        y_width = y_high - y_low
        z_width = z_high - z_low
    else:
        x_width, y_width, z_width = plot.xyz_widths
        x_low, x_high = -x_width, x_width
        y_low, y_high = -y_width, y_width
        z_low, z_high = -z_width, z_width

    xy_aspect = y_width / x_width
    z_aspect = z_width / x_width

    layout['scene']['xaxis']['range'] = [x_low, x_high]
    layout['scene']['yaxis']['range'] = [y_low, y_high]
    layout['scene']['zaxis']['range'] = [z_low, z_high]
    layout['scene']['aspectratio']['x'] = 1
    layout['scene']['aspectratio']['y'] = xy_aspect
    layout['scene']['aspectratio']['z'] = z_aspect

    return layout


def plot_trace_tracks(sc, fade_in_time=0):
    """
    Plots the tracks of a star cluster over time < 0 in a Scatter3d trace.
    The size of markers changes with time in proportion to the cluster's 'max_size'.
    """
    df_int = sc.df_int

    max_size = sc.max_size/2
    min_size = sc.min_size/2
    df_int = df_int.loc[(df_int['time'] <= 0) & (df_int['time'] > -1 * df_int['age_myr'] - fade_in_time)]

    size_fade = min_size + (max_size - min_size) * (1 - np.abs(df_int['time']) / (df_int['age_myr'] + fade_in_time))

    tracks = go.Scatter3d(
        x=df_int['x'].iloc[::1],
        y=df_int['y'].iloc[::1],
        z=df_int['z'].iloc[::1],
        mode='markers',
        marker=dict(
            size=size_fade,
            color=sc.color,
            opacity=sc.opacity / 1.5,
            line=dict(width=0)
        ),
        hoverinfo='none',
        name=sc.data_name + ' Track'
    )
    return tracks