import numpy as np
import yaml
import plotly.graph_objects as go
import plotly.colors as pc
from astropy.coordinates import SkyCoord
import astropy.units as u
import importlib.resources
from . import orbit_maker
import copy

GALACTIC_GUIDE_TRACE_NAMES = {
    'Galactic Quadrants',
    'Galactic l Ticks',
    'Galactic l Labels',
    'Galactic Z Axis'
}
GALACTIC_RADIUS_CIRCLE_DEFS = (
    (4.0, 'R = 4 kpc'),
    (8.122, 'R = 8.12 kpc'),
    (12.0, 'R = 12 kpc'),
)
GALACTIC_RADIUS_TRACE_NAMES = (
    {name for _, name in GALACTIC_RADIUS_CIRCLE_DEFS}
    | {f'{name} Label' for _, name in GALACTIC_RADIUS_CIRCLE_DEFS}
)

SPIRAL_ARMS = {
    "Perseus": {
        "theta_ref_deg": -13.0,
        "theta_range_deg": (-20.9, 88.2),
        "Rref_kpc": 10.88,
        "psi_deg": 9.8,
        "Omega_p": 17.82
    },
    "Local": {
        "theta_ref_deg": -2.3,
        "theta_range_deg": (-26.9, 26.6),
        "Rref_kpc": 8.69,
        "psi_deg": 8.9,
        "Omega_p": 33.76
    },
    "Sagittarius": {
        "theta_ref_deg": 3.5,
        "theta_range_deg": (-39.3, 67.7),
        "Rref_kpc": 7.10,
        "psi_deg": 10.6,
        "Omega_p": 26.10
    },
    "Scutum": {
        "theta_ref_deg": -4.8,
        "theta_range_deg": (-32.7, 100.9),
        "Rref_kpc": 6.02,
        "psi_deg": 14.9,
        "Omega_p": 49.81
    },
}
SPIRAL_ARM_TRACE_NAMES = {f"Spiral Arm: {name}" for name in SPIRAL_ARMS}
SEC_PER_MYR = 1e6 * 365.25 * 24 * 3600.0
KPC_IN_KM = 3.085677581e16
KM_S_PER_KPC_TO_RAD_MYR = (1.0 / KPC_IN_KM) * SEC_PER_MYR

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
        trace_grouping_dict=None,
        potential=None,
        vo=236.,
        ro=8.122,
        zo=0.0208
    ):
        self.data_collection = data_collection
        self.potential = potential
        self.vo = vo
        self.ro = ro
        self.zo = zo
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
        galactic_mode=False,
        fade_in_time=5,
        fade_in_and_out=False,
        fade_in_and_disp=False,
        disp_time=0,
        show_gc_line=True,
        coord_system='centered',
        show_galactic_guides=True,
        show_galactic_center_circles=True,
        include_spiral_arms=False,
        camera_zoom_factor=1.0,
        galactic_reference_opacity=0.5
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
                time, 
                reference_frame_center=reference_frame_center,
                potential=self.potential,
                vo=self.vo, ro=self.ro, zo=self.zo
            )
        # Set cluster sizes for fade effects
        self.fade_in_time = fade_in_time
        self.data_collection.set_all_cluster_sizes(
            self.fade_in_time,
            fade_in_and_out,
            fade_in_and_disp,
            disp_time
        )

        # Prepare time arrays and figure layout
        self.time = np.array(self.data_collection.time, dtype=np.float64)

        if camera_zoom_factor <= 0:
            raise ValueError("camera_zoom_factor must be > 0.")
        if (galactic_reference_opacity < 0) or (galactic_reference_opacity > 1):
            raise ValueError("galactic_reference_opacity must be between 0 and 1.")

        self.galactic_reference_opacity = float(galactic_reference_opacity)

        # Build the figure layout, optionally overriding for galactic mode.
        layout_dict = copy.deepcopy(self.figure_layout_dict)

        # Apply initial camera zoom by scaling eye distance from the scene center.
        # Larger factor => zoom in (camera moves closer), smaller => zoom out.
        if 'scene' in layout_dict:
            camera_dict = layout_dict['scene'].setdefault('camera', {})
            eye_dict = camera_dict.get('eye')
            if isinstance(eye_dict, dict):
                for axis_key in ('x', 'y', 'z'):
                    if axis_key in eye_dict and eye_dict[axis_key] is not None:
                        eye_dict[axis_key] = float(eye_dict[axis_key]) / float(camera_zoom_factor)

        # Keep xyz axis lines stylistically aligned with each other.
        self._sync_scene_axis_style(layout_dict)

        if galactic_mode:
            # Force symmetric x/y ranges to +/- 10 kpc (in pc).
            xy_half_range = 10000
            layout_dict['scene']['xaxis']['range'] = [-xy_half_range, xy_half_range]
            layout_dict['scene']['yaxis']['range'] = [-xy_half_range, xy_half_range]

            # Recompute aspect ratio using the existing z range.
            z_low, z_high = layout_dict['scene']['zaxis'].get('range', [-300, 300])
            z_width = float(z_high) - float(z_low)
            x_width = 2.0 * float(xy_half_range)
            layout_dict['scene']['aspectratio']['x'] = 1
            layout_dict['scene']['aspectratio']['y'] = 1
            layout_dict['scene']['aspectratio']['z'] = z_width / x_width

            # Remove 3D axes lines/ticks/labels entirely.
            layout_dict['scene']['xaxis']['visible'] = False
            layout_dict['scene']['yaxis']['visible'] = False
            layout_dict['scene']['zaxis']['visible'] = False

        self.figure_layout = go.Layout(layout_dict)
        self.focus_group = focus_group
        self.static_traces_legendonly = static_traces_legendonly

        if coord_system not in ('centered', 'rot'):
            raise ValueError("coord_system must be either 'centered' or 'rot'")

        self.coord_system = coord_system

        x_rf_int, y_rf_int, z_rf_int = orbit_maker.get_center_orbit_coords(
            time, reference_frame_center, potential=self.potential, vo=self.vo, ro=self.ro, zo=self.zo
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
                track = plot_trace_tracks(sc, self.fade_in_time, coord_system=self.coord_system)
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
                cluster_groups,
                t_i,
                x_rf_int[i],
                y_rf_int[i],
                z_rf_int[i],
                show_gc_line,
                galactic_mode,
                show_galactic_guides=show_galactic_guides,
                show_galactic_center_circles=show_galactic_center_circles,
                include_spiral_arms=include_spiral_arms,
                coord_system=self.coord_system
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

    def _coordFIX_to_coordROT(self, x_gc_pc, y_gc_pc, z_gc_pc, time_myr):
        """Apply the same fixed->rotating transform used in orbit_maker.coordFIX_to_coordROT.

        This avoids importing pandas here and keeps the GC-line transform consistent.
        """
        w0 = self.vo / self.ro
        w1 = w0 / 10.0
        t1 = time_myr * 0.01022
        r_sun_pc = self.ro * 1000.0

        r = np.sqrt(x_gc_pc**2 + y_gc_pc**2)
        theta = np.arctan2(x_gc_pc, y_gc_pc)

        x_rot = r_sun_pc - r * np.cos(theta - w1 * t1 + np.pi / 2.0)
        y_rot = r * np.sin(theta - w1 * t1 + np.pi / 2.0)
        z_rot = z_gc_pc
        return x_rot, y_rot, z_rot

    def rotating_gc_line_rot(self, t_myr):
        """GC reference line in the same rotating coordinate system as x_rot/y_rot/z_rot."""
        return self._radius_circle_trace(
            radius_kpc=8.122,
            trace_name='R = 8.12 kpc',
            coord_system='rot',
            t_myr=float(t_myr)
        )

    def rotating_gc_line(self, x_sub, y_sub, z_sub=0.0):
        """
        Generates a rotating galactic center line, displayed as a 3D line in the figure.
        """
        return self._radius_circle_trace(
            radius_kpc=8.122,
            trace_name='R = 8.12 kpc',
            coord_system='centered',
            x_sub=float(x_sub),
            y_sub=float(y_sub),
            z_sub=float(z_sub)
        )

    def _radius_circle_trace(
        self,
        radius_kpc,
        trace_name,
        coord_system='centered',
        t_myr=0.0,
        x_sub=0.0,
        y_sub=0.0,
        z_sub=0.0
    ):
        """Create a galactocentric radius circle trace in the selected coordinate system."""
        n_marks = 1000
        R = -float(radius_kpc) * np.ones(n_marks)
        phi = np.linspace(-180, 180, n_marks)
        radius_circle = SkyCoord(
            rho=R * u.kpc,
            phi=phi * u.deg,
            z=[0.] * len(R) * u.pc,
            frame='galactocentric',
            representation_type='cylindrical'
        )

        if coord_system == 'rot':
            x_gc_pc = radius_circle.galactocentric.cartesian.x.value * 1000.0
            y_gc_pc = radius_circle.galactocentric.cartesian.y.value * 1000.0
            z_gc_pc = radius_circle.galactocentric.cartesian.z.value * 1000.0
            x_vals, y_vals, z_vals = self._coordFIX_to_coordROT(
                x_gc_pc, y_gc_pc, z_gc_pc, float(t_myr)
            )
        else:
            x_vals = radius_circle.galactic.cartesian.x.value * 1000.0 - float(x_sub)
            y_vals = radius_circle.galactic.cartesian.y.value * 1000.0 - float(y_sub)
            z_vals = radius_circle.galactic.cartesian.z.value * 1000.0 - float(z_sub)

        line_color = self._reference_line_color()
        return go.Scatter3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            mode='lines',
            line=dict(color=line_color, width=self._reference_line_width(), dash='solid'),
            opacity=self._reference_opacity(),
            visible=True,
            name=trace_name,
            showlegend=False,
            hovertext=trace_name
        )

    def _reference_line_color(self):
        """Common reference-line color used by GC circle and galactic guides."""
        return 'gray' if self.figure_theme == 'dark' else 'black'

    def _reference_line_width(self):
        """Common line width for galactic reference overlays."""
        return 2.0

    def _reference_opacity(self):
        """Common opacity for galactic reference overlays."""
        return float(getattr(self, 'galactic_reference_opacity', 0.5))

    def _galactic_center_position(self, t, x_rf, y_rf, z_rf, coord_system='centered'):
        """Galactic-center position in the active reference frame at time t."""
        if coord_system == 'rot':
            x_gc, y_gc, z_gc = self._coordFIX_to_coordROT(
                np.array([0.0]), np.array([0.0]), np.array([0.0]), float(t)
            )
            return float(x_gc[0]), float(y_gc[0]), float(z_gc[0])
        return (
            (self.ro * 1000.0) - float(x_rf),
            0.0 - float(y_rf),
            0.0 - float(z_rf)
        )

    def _radius_label_trace(
        self,
        radius_kpc,
        label_text,
        x_center,
        y_center,
        z_center,
        angle_deg=315.0
    ):
        """Place a radius label at a defined angle where 0 deg points +y away from GC."""
        angle_rad = np.deg2rad(float(angle_deg))
        radius_pc = float(radius_kpc) * 1000.0
        x_label = float(x_center) + radius_pc * np.sin(angle_rad)
        y_label = float(y_center) + radius_pc * np.cos(angle_rad)
        z_label = float(z_center)

        return go.Scatter3d(
            x=[x_label],
            y=[y_label],
            z=[z_label],
            mode='text',
            text=[label_text],
            textposition='middle left',
            textfont=dict(color=self._reference_line_color(), size=12, family='helvetica'),
            opacity=self._reference_opacity(),
            name=f'{label_text} Label',
            showlegend=False,
            hoverinfo='skip'
        )

    def _build_galactic_circles_with_labels(self, t, x_rf, y_rf, z_rf, coord_system='centered'):
        """Build R=4/8.12/12 kpc circles and their labels."""
        x_gc, y_gc, z_gc = self._galactic_center_position(
            t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
        )
        traces = []
        for radius_kpc, label_text in GALACTIC_RADIUS_CIRCLE_DEFS:
            traces.append(
                self._radius_circle_trace(
                    radius_kpc=radius_kpc,
                    trace_name=label_text,
                    coord_system=coord_system,
                    t_myr=float(t),
                    x_sub=float(x_rf),
                    y_sub=float(y_rf),
                    z_sub=float(z_rf)
                )
            )
            traces.append(
                self._radius_label_trace(
                    radius_kpc=radius_kpc,
                    label_text=label_text,
                    x_center=x_gc,
                    y_center=y_gc,
                    z_center=z_gc,
                    angle_deg=315.0
                )
            )
        return traces

    def _log_spiral_radius(self, theta_rad, rref_kpc, theta_ref_rad, psi_rad):
        """Log-spiral radius model: ln(R/Rref) = -(theta - theta_ref) * tan(psi)."""
        return float(rref_kpc) * np.exp(-(theta_rad - theta_ref_rad) * np.tan(psi_rad))

    def _spiral_arm_coords_at_time(self, arm_params, t_myr, npts=500):
        """Compute galactocentric and heliocentric-flipped spiral-arm coordinates at time t."""
        theta_ref = np.deg2rad(float(arm_params["theta_ref_deg"]))
        th0_deg, th1_deg = arm_params["theta_range_deg"]
        theta_min = np.deg2rad(float(th0_deg))
        theta_max = np.deg2rad(float(th1_deg))
        rref_kpc = float(arm_params["Rref_kpc"])
        psi_rad = np.deg2rad(float(arm_params["psi_deg"]))
        omega_p = float(arm_params["Omega_p"])  # km/s/kpc
        omega_rad_myr = omega_p * KM_S_PER_KPC_TO_RAD_MYR
        # OVIZ uses lookback time as negative t; convert to positive lookback
        # so the arm evolution matches the intended pattern-speed convention.
        lookback_myr = -float(t_myr)
        dtheta = omega_rad_myr * lookback_myr

        theta_ref_t = theta_ref - dtheta
        theta_t = np.linspace(theta_min - dtheta, theta_max - dtheta, int(npts))
        r_kpc_t = self._log_spiral_radius(theta_t, rref_kpc, theta_ref_t, psi_rad)
        r_pc_t = r_kpc_t * 1000.0

        # Galactocentric Cartesian in the paper's convention.
        x_gc = r_pc_t * np.cos(theta_t)
        y_gc = r_pc_t * np.sin(theta_t)
        z_gc = np.zeros_like(x_gc)

        # Convert to heliocentric XY with flipped x so GC lies at +R0 on x.
        r0_pc = float(self.ro) * 1000.0
        x_helio = -(x_gc - r0_pc)
        y_helio = y_gc
        return x_gc, y_gc, z_gc, x_helio, y_helio

    def _build_spiral_arm_traces(self, t, x_rf, y_rf, z_rf, coord_system='centered'):
        """Build moving spiral-arm line traces for the current time."""
        traces = []
        for arm_name, arm_params in SPIRAL_ARMS.items():
            x_gc, y_gc, z_gc, x_helio, y_helio = self._spiral_arm_coords_at_time(
                arm_params, t_myr=t, npts=500
            )

            if coord_system == 'rot':
                x_vals, y_vals, z_vals = self._coordFIX_to_coordROT(
                    x_gc, y_gc, z_gc, float(t)
                )
            else:
                x_vals = x_helio - float(x_rf)
                y_vals = y_helio - float(y_rf)
                z_vals = z_gc - float(z_rf)

            traces.append(
                go.Scatter3d(
                    x=x_vals,
                    y=y_vals,
                    z=z_vals,
                    mode='lines',
                    line=dict(color='cyan', width=14.0, dash='solid'),
                    opacity=0.3,
                    name=f'Spiral Arm: {arm_name}',
                    showlegend=True,
                    hoverinfo='skip'
                )
            )

        return traces

    def _sync_scene_axis_style(self, layout_dict):
        """Apply a shared x-axis line style to y/z so xyz axis lines stay visually consistent."""
        scene = layout_dict.get('scene', {})
        xaxis = scene.get('xaxis', {})
        linecolor = xaxis.get('linecolor')
        linewidth = xaxis.get('linewidth')
        showline = xaxis.get('showline')

        for axis_name in ('yaxis', 'zaxis'):
            axis = scene.setdefault(axis_name, {})
            if linecolor is not None:
                axis['linecolor'] = linecolor
            if linewidth is not None:
                axis['linewidth'] = linewidth
            if showline is not None:
                axis['showline'] = showline

    def _build_galactic_guide_traces(self, t, sun_x=0.0, sun_y=0.0):
        """Quadrant guides, l-ticks, and z-axis line shown only at t=0."""
        show_guides = np.isclose(float(t), 0.0, rtol=0.0, atol=1e-9)
        guide_color = self._reference_line_color()
        guide_width = self._reference_line_width()
        guide_opacity = self._reference_opacity()

        try:
            x_min, x_max = [float(v) for v in self.figure_layout['scene']['xaxis']['range']]
        except Exception:
            x_min, x_max = -10000.0, 10000.0
        try:
            y_min, y_max = [float(v) for v in self.figure_layout['scene']['yaxis']['range']]
        except Exception:
            y_min, y_max = -10000.0, 10000.0
        try:
            z_min, z_max = [float(v) for v in self.figure_layout['scene']['zaxis']['range']]
        except Exception:
            z_min, z_max = -300.0, 300.0

        if not show_guides:
            return [
                go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='lines',
                    name='Galactic Quadrants',
                    showlegend=False,
                    hoverinfo='skip'
                ),
                go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='lines',
                    name='Galactic l Ticks',
                    showlegend=False,
                    hoverinfo='skip'
                ),
                go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='text',
                    text=[],
                    name='Galactic l Labels',
                    showlegend=False,
                    hoverinfo='skip'
                ),
                go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='lines',
                    name='Galactic Z Axis',
                    showlegend=False,
                    hoverinfo='skip'
                ),
            ]

        # Cross-hair lines dividing the Galactic plane into quadrants around the Sun.
        quadrants = go.Scatter3d(
            x=[x_min, x_max, None, sun_x, sun_x],
            y=[sun_y, sun_y, None, y_min, y_max],
            z=[0.0, 0.0, None, 0.0, 0.0],
            mode='lines',
            line=dict(color=guide_color, width=guide_width, dash='solid'),
            opacity=guide_opacity,
            name='Galactic Quadrants',
            showlegend=False,
            hoverinfo='skip'
        )

        # Tick marks and labels for l = 0, 90, 180, 270 around the Sun.
        xy_span = min(x_max - x_min, y_max - y_min)
        tick_radius = 0.055 * xy_span
        tick_half = 0.007 * xy_span
        label_gap = 0.010 * xy_span
        label_z = (0.03 * (z_max - z_min)) + 50.0
        angles = [0.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0]
        labels = ['l=0', 'l=90', 'l=180', 'l=270']

        x_ticks, y_ticks, z_ticks = [], [], []
        x_labels, y_labels, z_labels = [], [], []

        for angle in angles:
            c = np.cos(angle)
            s = np.sin(angle)
            r0 = tick_radius - tick_half
            r1 = tick_radius + tick_half
            x_ticks.extend([sun_x + r0 * c, sun_x + r1 * c, None])
            y_ticks.extend([sun_y + r0 * s, sun_y + r1 * s, None])
            z_ticks.extend([0.0, 0.0, None])
            x_labels.append(sun_x + (tick_radius + tick_half + label_gap) * c)
            y_labels.append(sun_y + (tick_radius + tick_half + label_gap) * s)
            z_labels.append(label_z)

        l_ticks = go.Scatter3d(
            x=x_ticks,
            y=y_ticks,
            z=z_ticks,
            mode='lines',
            line=dict(color=guide_color, width=guide_width, dash='solid'),
            opacity=guide_opacity,
            name='Galactic l Ticks',
            showlegend=False,
            hoverinfo='skip'
        )

        l_labels = go.Scatter3d(
            x=x_labels,
            y=y_labels,
            z=z_labels,
            mode='text',
            text=labels,
            textposition='middle center',
            textfont=dict(color=guide_color, size=12, family='helvetica'),
            name='Galactic l Labels',
            showlegend=False,
            hoverinfo='skip'
        )

        z_axis = go.Scatter3d(
            x=[sun_x, sun_x],
            y=[sun_y, sun_y],
            z=[z_min, z_max],
            mode='lines',
            line=dict(color=guide_color, width=guide_width, dash='solid'),
            opacity=guide_opacity,
            name='Galactic Z Axis',
            showlegend=False,
            hoverinfo='skip'
        )

        return [quadrants, l_ticks, l_labels, z_axis]

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
            borderwidth=0.,
            bordercolor=slider_color,
            bgcolor=slider_color,
            pad={"b": 0, "t": 0, "l": 0, "r": 0},
            len=0.5,  # Adjusted the length to 0.5 for better centering
            x=0.5,
            y=0.,
            currentvalue={
                "font": {"size": 18, "color": slider_color, 'family': 'helvetica'},
                'prefix': 'Time (Myr): ',
                'visible': True,
                'xanchor': 'center',
                "offset": 20
            },
            steps=[
                dict(
                args=[
                    [str(t)],
                    dict(
                    frame=dict(duration=5, easing='linear', redraw=True),
                    transition=dict(duration=0, easing='linear')
                    )
                ],
                label=str(t),
                method='animate'
                ) for t in time_slider
            ],
            tickcolor=slider_color,
            ticklen=10,
            font=dict(color='rgba(0,0,0,0)', size=8, family='helvetica')
            )
        ]

    def dropdown_menu(self):
        """
        Creates a dropdown menu for selecting trace groups in the 3D plot.
        Uses 'method':'update' so the user’s grouping choice persists 
        even when the slider returns to t=0.
        In this update we ensure that if a static trace is a track (its name ends with ' Track'),
        its visibility is determined by whether its base name is in the grouping.
        """
        buttons = []
        for key, traces_list in self.trace_grouping_dict.items():
            visibility = []
            for trace in self.initial_data:
                # Use our get_visibility function so that static tracks (and non‑track statics)
                # obey the static_traces_legendonly flag.
                visibility.append(self.get_visibility(trace['name'], traces_list))
            buttons.append(
                dict(
                    label=key,
                    method='update',
                    args=[
                        {"visible": visibility},
                        {}
                    ]
                )
            )

        bg_color = self.figure_layout_dict['scene']['xaxis']['backgroundcolor']
        text_color = 'black' if self.figure_theme != 'dark' else 'white'
        for button in buttons:
            button['args'][1]['hoverlabel'] = dict(bgcolor='gray')

        return [
            dict(
                buttons=buttons,
                direction='down',
                pad={'r': 10, 't': 10},
                showactive=True,
                x=0,
                xanchor='left',
                y=1.1,
                yanchor='top',
                bgcolor=bg_color,
                font=dict(color=text_color, family='helvetica', size=14),
                active=0
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
        - Galactic reference overlays (GC, radius circles/labels, guide overlays) always return True.
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

        # In galactic mode, GC and galactic radius circles/labels stay visible.
        if trace_name == 'GC':
            return True
        if trace_name in GALACTIC_RADIUS_TRACE_NAMES:
            return True

        # Galactic-mode reference guides should not depend on grouping.
        if trace_name in GALACTIC_GUIDE_TRACE_NAMES:
            return True
        if trace_name in SPIRAL_ARM_TRACE_NAMES:
            return True
        # For non-static traces:
        elif trace_name in grouping:
            return True

        return False

    def _generate_scatter_list(
        self,
        cluster_groups,
        t,
        x_rf,
        y_rf,
        z_rf,
        show_gc_line,
        galactic_mode,
        show_galactic_guides=True,
        show_galactic_center_circles=True,
        include_spiral_arms=False,
        coord_system='centered'
    ):
        """
        Generate the main cluster Scatter3d traces for a given time.
        Optionally includes the non-galactic R=8.12 line when show_gc_line is True.
        In galactic mode, dedicated toggles control guide overlays and GC/circle overlays.
        Non-static traces need not be “marked” since the dropdown update will control them.
        """
        scatter_list = []
        if coord_system == 'rot':
            x_col, y_col, z_col = 'x_rot', 'y_rot', 'z_rot'
        else:
            x_col, y_col, z_col = 'x', 'y', 'z'

        sun_x = 0.0
        sun_y = 0.0

        for cluster_group in cluster_groups:
            assert cluster_group.integrated
            df_int = cluster_group.df_int
            if df_int.empty:
                continue

            # Use isclose to avoid float equality pitfalls (e.g. -10.0 vs -10).
            time_mask = np.isclose(df_int['time'].to_numpy(dtype=float), float(t), rtol=0.0, atol=1e-9)
            df_t = df_int[time_mask]
            if df_t.empty:
                continue

            if cluster_group.data_name.strip().lower() == 'sun':
                sun_x = float(np.nanmedian(df_t[x_col].to_numpy(dtype=float)))
                sun_y = float(np.nanmedian(df_t[y_col].to_numpy(dtype=float)))

            age_at_t = df_t['age_myr'] + t
            age_present = df_t['age_myr']
            hovertext = (
                '<b style="font-size:16px;">' + df_t['name'].str.replace('_', ' ').astype(str) + '</b>' + '<br>'  # Bold cluster name
                + cluster_group.data_name + '<br>'  # Group name
                + 'Age (now) = ' + age_present.round(1).astype(str) + ' Myr' + '<br>'
                + 'Age (t) = ' + age_at_t.round(1).astype(str) + ' Myr' + '<br>'  # Cluster age at time t
            )

            if 'n_stars' in df_t.columns:
                hovertext += 'N = ' + df_t['n_stars'].astype(str) + ' stars <br>'  # Number of stars

            hovertext += (
                f'({x_col},{y_col},{z_col}) = (' +
                df_t[x_col].round(1).astype(str) + ', ' +
                df_t[y_col].round(1).astype(str) + ', ' +
                df_t[z_col].round(1).astype(str) + ')'
            )

            marker_dict = dict(
                size=df_t['size'],
                symbol=cluster_group.marker_style,
                line=dict(color='black', width=0.0)
            )

            marker_dict.update(opacity=cluster_group.opacity)
            if cluster_group.colormap:
                age_full = df_int['age_myr'] + df_int['time']
                cmin = cluster_group.cmin if cluster_group.cmin is not None else float(age_full.min())
                cmax = cluster_group.cmax if cluster_group.cmax is not None else float(age_full.max())
                marker_dict.update(
                    color=age_at_t.values,
                    colorscale=cluster_group.colormap,
                    cmin=cmin,
                    cmax=cmax
                )
            else:
                marker_dict.update(color=cluster_group.color)

            scatter_list.append(
                go.Scatter3d(
                    x=df_t[x_col].values,
                    y=df_t[y_col].values,
                    z=df_t[z_col].values,
                    mode='markers',
                    marker=marker_dict,
                    customdata=np.column_stack((age_present.to_numpy(dtype=float), age_at_t.to_numpy(dtype=float))),
                    meta={'trace_kind': 'cluster'},
                    hovertext=hovertext,
                    hoverinfo='text',  # This removes default x, y, z
                    hovertemplate='%{hovertext}<extra></extra>',  # This ensures only custom hovertext is shown
                    name=cluster_group.data_name
                )
            )

        show_reference_lines = show_gc_line if not galactic_mode else show_galactic_center_circles

        if show_reference_lines:
            if galactic_mode:
                scatter_list.extend(
                    self._build_galactic_circles_with_labels(
                        t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
                    )
                )
            else:
                if coord_system == 'rot':
                    gc_line_t = self.rotating_gc_line_rot(t)
                else:
                    gc_line_t = self.rotating_gc_line(x_rf, y_rf, z_rf)
                scatter_list.append(gc_line_t)

        if galactic_mode and show_galactic_guides:
            scatter_list.extend(self._build_galactic_guide_traces(t, sun_x=sun_x, sun_y=sun_y))

        if galactic_mode and include_spiral_arms:
            scatter_list.extend(
                self._build_spiral_arm_traces(
                    t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
                )
            )

        if galactic_mode and show_galactic_center_circles:
            # Theme-aware styling for the GC marker/text.
            if self.figure_theme == 'dark':
                gc_marker_color = 'gray'
                gc_text_color = 'gray'
                gc_line_color = 'black'
            elif self.figure_theme in ('light', 'solarized_light'):
                gc_marker_color = 'black'
                gc_text_color = 'black'
                gc_line_color = 'white'
            else:
                # Fallback for other themes (e.g. 'gray')
                gc_marker_color = 'black'
                gc_text_color = 'black'
                gc_line_color = 'white'

            x_gc, y_gc, z_gc = self._galactic_center_position(
                t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
            )

            scatter_list.append(
                go.Scatter3d(
                    x=[x_gc],
                    y=[y_gc],
                    z=[z_gc],
                    mode='markers+text',
                    marker=dict(
                        size=12,
                        color=gc_marker_color,
                        symbol='circle',
                        opacity=self._reference_opacity(),
                        line=dict(color=gc_line_color, width=2)
                    ),
                    text=['GC'],
                    textposition='top center',
                    textfont=dict(color=gc_text_color, size=14, family='helvetica'),
                    name='GC',
                    showlegend=False,
                    hovertext='Galactic Center',
                    hoverinfo='text',
                    hovertemplate='%{hovertext}<extra></extra>'
                )
            )

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

    # Backward-compatible alias (in case external code used the old name)
    initialize_figure = _initialize_figure

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


def plot_trace_tracks(sc, fade_in_time=0, coord_system='centered'):
    """
    Plots the tracks of a star cluster over time < 0 in a Scatter3d trace.
    The size of markers changes with time in proportion to the cluster's 'max_size'.
    """
    df_int = sc.df_int

    max_size = sc.max_size/2
    min_size = sc.min_size/2
    df_int = df_int.loc[(df_int['time'] <= 0) & (df_int['time'] > -1 * df_int['age_myr'] - fade_in_time)]

    size_fade = min_size + (max_size - min_size) * (1 - np.abs(df_int['time']) / (df_int['age_myr'] + fade_in_time))

    if coord_system == 'rot':
        x_col, y_col, z_col = 'x_rot', 'y_rot', 'z_rot'
    else:
        x_col, y_col, z_col = 'x', 'y', 'z'

    tracks = go.Scatter3d(
        x=df_int[x_col].iloc[::1],
        y=df_int[y_col].iloc[::1],
        z=df_int[z_col].iloc[::1],
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
