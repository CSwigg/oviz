import base64
import json
import copy
import io

import numpy as np
import yaml
import plotly.graph_objects as go
import plotly.colors as pc
import webcolors
from astropy.coordinates import SkyCoord
import astropy.units as u
import importlib.resources
from . import orbit_maker
import pandas as pd
from .threejs_figure import ThreeJSFigure

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
KDE_TRACE_PREFIX = 'Age KDE: '
KDE_TIME_MARKER_TRACE_NAME = 'Age KDE Time Marker'
CUSTOMDATA_IDX_AGE_NOW = 0
CUSTOMDATA_IDX_AGE_AT_T = 1
CUSTOMDATA_IDX_L0_DEG = 2
CUSTOMDATA_IDX_B0_DEG = 3
CUSTOMDATA_IDX_DIST0_PC = 4
CUSTOMDATA_IDX_X0 = 5
CUSTOMDATA_IDX_Y0 = 6
CUSTOMDATA_IDX_Z0 = 7
CUSTOMDATA_IDX_CLUSTER_NAME = 8
CUSTOMDATA_IDX_CLUSTER_COLOR = 9
MAX_SELECTED_MEMBER_POINTS = 1200
DEFAULT_THREEJS_VOLUME_COLORMAPS = (
    'inferno',
    'magma',
    'plasma',
    'viridis',
    'cividis',
    'turbo',
    'Greys',
)

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
        show_age_kde_inset=False,
        age_kde_bandwidth_myr=3.0,
        camera_zoom_factor=1.0,
        galactic_reference_opacity=0.5,
        renderer='plotly',
        show_milky_way_model=False,
        enable_sky_panel=False,
        sky_radius_deg=1.0,
        sky_frame='galactic',
        sky_survey='P/DSS2/color',
        cluster_members_file=None,
        volumes=None,
    ):
        """
        Main method to generate a 3D animated plot. Integrates orbits if necessary and
        constructs frames for each time step. Optionally saves the figure to HTML.
        """
        renderer_name = _normalize_renderer_name(renderer)

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
        self.fade_in_and_out = bool(fade_in_and_out)
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
        self.show_age_kde_inset = bool(show_age_kde_inset)
        self.show_milky_way_model = bool(show_milky_way_model)
        self.enable_sky_panel = bool(enable_sky_panel)
        self.sky_radius_deg = float(sky_radius_deg)
        self.sky_frame = str(sky_frame)
        self.sky_survey = str(sky_survey)
        self.cluster_members_file = cluster_members_file
        self.sky_members_by_cluster = None
        self.volume_configs = _normalize_threejs_volume_configs(volumes)
        if age_kde_bandwidth_myr <= 0:
            raise ValueError("age_kde_bandwidth_myr must be > 0.")
        self.age_kde_bandwidth_myr = float(age_kde_bandwidth_myr)
        if renderer_name != 'threejs' and self.volume_configs:
            raise NotImplementedError(
                "volumes are currently only supported with renderer='threejs'."
            )
        if self.enable_sky_panel:
            if self.sky_radius_deg <= 0:
                raise ValueError("sky_radius_deg must be > 0.")
            if self.sky_frame not in ('galactic', 'icrs'):
                raise ValueError("sky_frame must be either 'galactic' or 'icrs'.")
            if self.sky_frame != 'galactic':
                raise ValueError("sky_frame='icrs' is not yet supported; use 'galactic'.")
            self.sky_members_by_cluster = _load_threejs_cluster_catalog(cluster_members_file)

        # Build the figure layout, optionally overriding for galactic mode.
        layout_dict = copy.deepcopy(self.figure_layout_dict)
        layout_dict['dragmode'] = 'turntable'
        layout_dict.setdefault('scene', {})['dragmode'] = 'turntable'

        if self.show_age_kde_inset:
            self._setup_age_kde_inset(layout_dict)

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
        if renderer_name == 'threejs' and self.volume_configs and self.coord_system != 'centered':
            raise NotImplementedError(
                "Three.js volume rendering currently supports coord_system='centered' only."
            )

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

        if renderer_name == 'plotly':
            # Add slider & dropdown, and finalize
            self._add_slider_and_dropdown()
        else:
            self._build_threejs_figure(frames)

        # Show or save
        if show:
            self.figure.show()
        if save_name:
            if renderer_name == 'plotly':
                self.figure.update_layout(width=None, height=None)
                post_script = self._legend_kde_sync_post_script() if self.show_age_kde_inset else None
                self.figure.write_html(save_name, auto_play=False, post_script=post_script)
            else:
                self.figure.write_html(save_name)

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

    def _legend_kde_sync_post_script(self):
        """JS hook for exported HTML: mirror 3D legend toggles to matching KDE traces."""
        script = """
(function() {
  const gd = document.getElementById('{plot_id}');
  if (!gd) return;

  const KDE_PREFIX = __KDE_PREFIX__;
  const KDE_MARKER_NAME = __KDE_MARKER_NAME__;
  let syncing = false;

  function findKdeIndexForTraceName(traceName) {
    const kdeName = KDE_PREFIX + traceName;
    for (let i = 0; i < gd.data.length; i++) {
      const tr = gd.data[i];
      if (tr && tr.name === kdeName) return i;
    }
    return -1;
  }

  gd.on('plotly_restyle', function(eventData) {
    if (syncing) return;
    if (!eventData || eventData.length < 2) return;

    const update = eventData[0] || {};
    const idxRaw = eventData[1] || [];
    if (!Object.prototype.hasOwnProperty.call(update, 'visible')) return;

    const idxs = Array.isArray(idxRaw) ? idxRaw : [idxRaw];
    const visRaw = update.visible;
    const visVals = Array.isArray(visRaw) ? visRaw : [visRaw];

    const kdeIdxs = [];
    const kdeVis = [];

    for (let j = 0; j < idxs.length; j++) {
      const iTrace = idxs[j];
      const tr = gd.data[iTrace];
      if (!tr || typeof tr.name !== 'string') continue;
      if (tr.name.startsWith(KDE_PREFIX)) continue;
      if (tr.name === KDE_MARKER_NAME) continue;

      const kdeIdx = findKdeIndexForTraceName(tr.name);
      if (kdeIdx < 0) continue;

      const vis = visVals[Math.min(j, visVals.length - 1)];
      kdeIdxs.push(kdeIdx);
      kdeVis.push(vis);
    }

    if (!kdeIdxs.length) return;

    syncing = true;
    Plotly.restyle(gd, { visible: kdeVis }, kdeIdxs)
      .then(function() { syncing = false; })
      .catch(function() { syncing = false; });
  });
})();
"""
        script = script.replace('__KDE_PREFIX__', json.dumps(KDE_TRACE_PREFIX))
        script = script.replace('__KDE_MARKER_NAME__', json.dumps(KDE_TIME_MARKER_TRACE_NAME))
        return script

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
        """Build R=4/8.12/12 kpc circles."""
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

    def _kde_trace_name(self, trace_name):
        return f'{KDE_TRACE_PREFIX}{trace_name}'

    def _kde_source_trace_name(self, trace_name):
        if isinstance(trace_name, str) and trace_name.startswith(KDE_TRACE_PREFIX):
            return trace_name[len(KDE_TRACE_PREFIX):]
        return None

    def _coerce_numeric(self, values):
        """Coerce iterable values to finite float array."""
        vals = pd.to_numeric(values, errors='coerce').to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        return vals

    def _kde_color_for_cluster(self, cluster_group):
        """Choose a representative color for a cluster KDE line."""
        color = getattr(cluster_group, 'color', None)
        if isinstance(color, str) and color:
            return color
        return 'white' if self.figure_theme == 'dark' else 'black'

    def _kde_opacity_for_cluster(self, cluster_group):
        """Use trace opacity for KDE line opacity."""
        opacity = getattr(cluster_group, 'opacity', 1.0)
        try:
            opacity = float(opacity)
        except (TypeError, ValueError):
            opacity = 1.0
        return max(0.0, min(1.0, opacity))

    def _color_with_alpha(self, color, alpha):
        """Convert common color forms to rgba(r,g,b,a) with requested alpha."""
        alpha = max(0.0, min(1.0, float(alpha)))

        if isinstance(color, (tuple, list)) and len(color) >= 3:
            r, g, b = int(color[0]), int(color[1]), int(color[2])
            return f'rgba({r},{g},{b},{alpha})'

        if not isinstance(color, str):
            if self.figure_theme == 'dark':
                return f'rgba(255,255,255,{alpha})'
            return f'rgba(0,0,0,{alpha})'

        c = color.strip()
        try:
            if c.startswith('rgba(') and c.endswith(')'):
                vals = [v.strip() for v in c[5:-1].split(',')]
                r, g, b = int(float(vals[0])), int(float(vals[1])), int(float(vals[2]))
                return f'rgba({r},{g},{b},{alpha})'
            if c.startswith('rgb(') and c.endswith(')'):
                vals = [v.strip() for v in c[4:-1].split(',')]
                r, g, b = int(float(vals[0])), int(float(vals[1])), int(float(vals[2]))
                return f'rgba({r},{g},{b},{alpha})'
            if c.startswith('#'):
                rgb = webcolors.hex_to_rgb(c)
                return f'rgba({rgb.red},{rgb.green},{rgb.blue},{alpha})'
            rgb = webcolors.name_to_rgb(c)
            return f'rgba({rgb.red},{rgb.green},{rgb.blue},{alpha})'
        except Exception:
            # Keep original color if parsing fails.
            return c

    def _collect_trace_lookback_times(self, max_lookback_myr=None):
        """Collect lookback times (-age_myr) per trace, optionally clipped by integration window."""
        lookback_by_trace = {}
        color_by_trace = {}
        opacity_by_trace = {}
        for cluster_group in self.data_collection.get_all_clusters():
            data_df = cluster_group.df_int if cluster_group.df_int is not None else cluster_group.df
            if data_df is None or ('age_myr' not in data_df.columns):
                continue
            vals = self._coerce_numeric(data_df['age_myr'])
            if max_lookback_myr is not None:
                vals = vals[vals <= float(max_lookback_myr) + 1e-9]
            if vals.size:
                lookback_by_trace[cluster_group.data_name] = -np.abs(vals)
                color_by_trace[cluster_group.data_name] = self._kde_color_for_cluster(cluster_group)
                opacity_by_trace[cluster_group.data_name] = self._kde_opacity_for_cluster(cluster_group)
        return lookback_by_trace, color_by_trace, opacity_by_trace

    def _gaussian_kde(self, values, x_grid, bandwidth_myr=None):
        """Lightweight Gaussian KDE (SciPy-free)."""
        vals = np.asarray(values, dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return np.zeros_like(x_grid, dtype=float)

        if bandwidth_myr is not None:
            bandwidth = float(bandwidth_myr)
        else:
            std = float(np.std(vals, ddof=1))
            iqr = float(np.subtract(*np.percentile(vals, [75.0, 25.0])))
            sigma = min(std, iqr / 1.34) if (std > 0 and iqr > 0) else std
            bandwidth = 0.9 * sigma * (vals.size ** (-1.0 / 5.0))
            if (not np.isfinite(bandwidth)) or (bandwidth <= 0):
                grid_span = float(np.max(x_grid) - np.min(x_grid))
                bandwidth = max(grid_span / 60.0, 1e-3)

        if (not np.isfinite(bandwidth)) or (bandwidth <= 0):
            bandwidth = 1.0

        u = (x_grid[:, None] - vals[None, :]) / bandwidth
        density = np.exp(-0.5 * (u ** 2)).sum(axis=1)
        density /= (vals.size * bandwidth * np.sqrt(2.0 * np.pi))
        return density

    def _setup_age_kde_inset(self, layout_dict):
        """Prepare KDE inset data/axes and precompute one KDE line per integrated cluster trace."""
        self.kde_trace_order = []
        self.kde_density_by_trace = {}
        self.kde_trace_name_by_trace = {}
        self.kde_color_by_trace = {}
        self.kde_opacity_by_trace = {}

        finite_time = self.time[np.isfinite(self.time)] if self.time.size else np.array([], dtype=float)
        non_positive_time = finite_time[finite_time <= 0]
        if non_positive_time.size:
            lookback_limit = abs(float(np.min(non_positive_time)))
        else:
            lookback_limit = None

        lookback_by_trace, color_by_trace, opacity_by_trace = self._collect_trace_lookback_times(
            max_lookback_myr=lookback_limit
        )
        self.kde_trace_order = list(lookback_by_trace.keys())
        self.kde_color_by_trace = color_by_trace
        self.kde_opacity_by_trace = opacity_by_trace

        if non_positive_time.size:
            time_min = float(np.min(non_positive_time))
        elif self.kde_trace_order:
            all_lookback_values = np.concatenate([lookback_by_trace[k] for k in self.kde_trace_order])
            time_min = float(np.min(all_lookback_values))
        else:
            time_min = 0.0

        # Star-formation-history x-axis should be non-positive lookback time only.
        x_min = min(time_min, 0.0)
        x_max = 0.0
        if x_min >= x_max:
            x_min = x_max - 1.0

        # Keep KDE x-range exactly matched to slider's time span for visual alignment.
        grid_min = x_min
        grid_max = x_max
        self.kde_x_grid = np.linspace(grid_min, grid_max, 300)
        self.kde_x_range = [float(grid_min), float(grid_max)]

        y_max = 0.0
        for trace_name in self.kde_trace_order:
            dens = self._gaussian_kde(
                lookback_by_trace[trace_name],
                self.kde_x_grid,
                bandwidth_myr=self.age_kde_bandwidth_myr
            )
            dens_max = float(np.max(dens)) if dens.size else 0.0
            if dens_max > 0:
                dens = dens / dens_max
            self.kde_density_by_trace[trace_name] = dens
            self.kde_trace_name_by_trace[trace_name] = self._kde_trace_name(trace_name)
            if dens.size:
                y_max = max(y_max, float(np.max(dens)))

        if y_max <= 0:
            y_max = 1.0
        self.kde_y_max = max(y_max * 1.05, 1.0)

        # Reserve a small lower strip so the inset cannot be occluded by the 3D scene.
        scene_domain = layout_dict.setdefault('scene', {}).setdefault('domain', {})
        scene_domain.setdefault('x', [0.0, 1.0])
        scene_domain['y'] = [0.19, 1.0]

        if self.figure_theme == 'dark':
            axis_color = 'gray'
        elif self.figure_theme == 'gray':
            axis_color = 'black'
        else:
            axis_color = 'black'

        self.kde_axis_color = axis_color

        # Center the inset horizontally; slider is tied to this domain.
        x2_domain = [0.30, 0.70]
        y2_domain = [0.03, 0.17]
        self.kde_x2_domain = list(x2_domain)
        self.kde_y2_domain = list(y2_domain)
        panel_bg = layout_dict.get('scene', {}).get('bgcolor', 'black')

        # Draw a compact background panel only under the inset area.
        layout_dict.setdefault('shapes', [])
        layout_dict['shapes'].append(
            dict(
                type='rect',
                xref='paper',
                yref='paper',
                x0=x2_domain[0],
                x1=x2_domain[1],
                y0=y2_domain[0],
                y1=y2_domain[1],
                fillcolor=panel_bg,
                line=dict(width=0),
                layer='below'
            )
        )

        layout_dict['xaxis2'] = dict(
            domain=x2_domain,
            anchor='y2',
            range=self.kde_x_range,
            title='',
            titlefont=dict(color=axis_color, size=10, family='helvetica'),
            tickfont=dict(color=axis_color, size=9, family='helvetica'),
            showgrid=False,
            zeroline=False,
            showline=True,
            linecolor=axis_color,
            linewidth=1.5,
            mirror=True,
            layer='above traces',
            fixedrange=True
        )
        layout_dict['yaxis2'] = dict(
            domain=y2_domain,
            anchor='x2',
            range=[0.0, self.kde_y_max],
            title='Relative SFH',
            titlefont=dict(color=axis_color, size=10, family='helvetica'),
            tickfont=dict(color=axis_color, size=9, family='helvetica'),
            showgrid=False,
            zeroline=False,
            showline=True,
            linecolor=axis_color,
            linewidth=1.5,
            mirror=True,
            layer='above traces',
            fixedrange=True
        )
        layout_dict.setdefault('annotations', [])

    def _build_kde_inset_traces(self):
        """Build one static KDE trace per integrated cluster trace."""
        traces = []
        for trace_name in self.kde_trace_order:
            base_color = self.kde_color_by_trace.get(trace_name, self.kde_axis_color)
            line_color = self._color_with_alpha(base_color, self.kde_opacity_by_trace.get(trace_name, 1.0))
            fill_color = self._color_with_alpha(base_color, 0.30)
            traces.append(
                go.Scatter(
                    x=self.kde_x_grid,
                    y=self.kde_density_by_trace[trace_name],
                    mode='lines',
                    line=dict(color=line_color, width=2),
                    fill='tozeroy',
                    fillcolor=fill_color,
                    xaxis='x2',
                    yaxis='y2',
                    name=self.kde_trace_name_by_trace[trace_name],
                    showlegend=False,
                    hovertemplate=(
                        f'{trace_name}<br>'
                        + 'Lookback = %{x:.1f} Myr<br>'
                        + 'Relative KDE = %{y:.2f}<extra></extra>'
                    )
                )
            )
        return traces

    def _build_kde_time_marker_trace(self, t):
        """Build moving vertical time marker for the KDE inset."""
        if not np.isfinite(float(t)):
            t = 0.0
        x_t = float(np.clip(float(t), self.kde_x_range[0], 0.0))
        return go.Scatter(
            x=[x_t, x_t],
            y=[0.0, self.kde_y_max],
            mode='lines',
            line=dict(color=self.kde_axis_color, width=2, dash='dash'),
            xaxis='x2',
            yaxis='y2',
            name=KDE_TIME_MARKER_TRACE_NAME,
            showlegend=False,
            hovertemplate='t = %{x:.1f} Myr<extra></extra>'
        )

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

        time_slider = self._ordered_slider_times()

        # Locate the index where time is zero
        zero_idx = np.where(time_slider == 0)[0][0]

        if self.show_age_kde_inset and hasattr(self, 'kde_x2_domain'):
            slider_x0 = float(self.kde_x2_domain[0])
            slider_x1 = float(self.kde_x2_domain[1])
            slider_len = max(0.0, slider_x1 - slider_x0)
            slider_x = slider_x0
            slider_xanchor = "left"
        else:
            slider_len = 0.5
            slider_x = 0.5
            slider_xanchor = "center"

        return [
            dict(
            active=zero_idx,
            xanchor=slider_xanchor,
            yanchor="top",
            transition={"duration": 300, "easing": "bounce-in"},
            borderwidth=0.,
            bordercolor=slider_color,
            bgcolor=slider_color,
            pad={"b": 0, "t": 0, "l": 0, "r": 0},
            len=slider_len,
            x=slider_x,
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

    def _ordered_slider_times(self):
        """Return times in the same order used by the Plotly slider."""
        time_neg = self.time[self.time < 0]
        time_pos = self.time[self.time >= 0]

        if (len(time_neg) > 0) and (len(time_pos) > 1):
            return np.append(time_neg, time_pos)
        if (len(time_neg) > 0) and (len(time_pos) == 1):
            return np.flip(self.time)
        return self.time

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
                trace_name = self._trace_name(trace)
                # Use our get_visibility function so that static tracks (and non‑track statics)
                # obey the static_traces_legendonly flag.
                visibility.append(self.get_visibility(trace_name, traces_list))
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

    def _trace_name(self, trace):
        """Return trace name for both graph_objects traces and dict-like traces."""
        if isinstance(trace, dict):
            return trace.get('name')
        return getattr(trace, 'name', None)

    def _set_trace_visible(self, trace, visible_flag):
        """Set trace visibility for both graph_objects traces and dict-like traces."""
        if isinstance(trace, dict):
            trace['visible'] = visible_flag
        else:
            trace.visible = visible_flag

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
        if trace_name is None:
            return False

        if trace_name == KDE_TIME_MARKER_TRACE_NAME:
            return True

        kde_source_trace = self._kde_source_trace_name(trace_name)
        if kde_source_trace is not None:
            return kde_source_trace in grouping

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

            # Build present-day (t=0) sky quantities for click->sky-panel callbacks.
            t0_mask = np.isclose(df_int['time'].to_numpy(dtype=float), 0.0, rtol=0.0, atol=1e-9)
            df_t0 = df_int[t0_mask]
            if len(df_t0) != len(df_t):
                # Fallback for any unexpected ordering/shape mismatch.
                df_t0 = df_t

            x0 = pd.to_numeric(df_t0[x_col], errors='coerce').to_numpy(dtype=float)
            y0 = pd.to_numeric(df_t0[y_col], errors='coerce').to_numpy(dtype=float)
            z0 = pd.to_numeric(df_t0[z_col], errors='coerce').to_numpy(dtype=float)

            x_helio0 = pd.to_numeric(df_t0['x_helio'], errors='coerce').to_numpy(dtype=float)
            y_helio0 = pd.to_numeric(df_t0['y_helio'], errors='coerce').to_numpy(dtype=float)
            z_helio0 = pd.to_numeric(df_t0['z_helio'], errors='coerce').to_numpy(dtype=float)
            dist0 = np.sqrt(x_helio0 ** 2 + y_helio0 ** 2 + z_helio0 ** 2)

            with np.errstate(invalid='ignore', divide='ignore'):
                l0 = np.rad2deg(np.arctan2(y_helio0, x_helio0))
                l0 = np.mod(l0, 360.0)
                b0 = np.rad2deg(np.arcsin(np.clip(z_helio0 / np.where(dist0 > 0, dist0, np.nan), -1.0, 1.0)))

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

            trace_color = cluster_group.color if isinstance(cluster_group.color, str) else 'white'
            cluster_names = df_t['name'].astype(str).to_numpy(dtype=object)
            trace_colors = np.repeat(trace_color, len(df_t)).astype(object)

            scatter_list.append(
                go.Scatter3d(
                    x=df_t[x_col].values,
                    y=df_t[y_col].values,
                    z=df_t[z_col].values,
                    mode='markers',
                    marker=marker_dict,
                    customdata=np.column_stack((
                        age_present.to_numpy(dtype=float),
                        age_at_t.to_numpy(dtype=float),
                        l0,
                        b0,
                        dist0,
                        x0,
                        y0,
                        z0,
                        cluster_names,
                        trace_colors
                    )),
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

        if self.show_age_kde_inset:
            scatter_list.append(self._build_kde_time_marker_trace(t))
            scatter_list.extend(self._build_kde_inset_traces())

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
            trace_name = self._trace_name(trace)
            visible_flag = self.get_visibility(trace_name, grouping_0)
            self._set_trace_visible(trace, visible_flag)
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

    def _build_threejs_figure(self, frames):
        """Build a standalone three.js figure wrapper from the current frame data."""
        scene_spec = self._build_threejs_scene_spec(frames)
        self.fig = scene_spec
        self.fig_dict = scene_spec
        self.figure = ThreeJSFigure(scene_spec)

    def _build_threejs_scene_spec(self, frames):
        """Serialize the current animation into a renderer-agnostic scene spec."""
        layout_json = self.figure_layout.to_plotly_json()
        scene_layout = layout_json.get('scene', {})
        x_range = _coerce_range(scene_layout.get('xaxis', {}).get('range'), [-1.0, 1.0])
        y_range = _coerce_range(scene_layout.get('yaxis', {}).get('range'), [-1.0, 1.0])
        z_range = _coerce_range(scene_layout.get('zaxis', {}).get('range'), [-1.0, 1.0])
        center = {
            'x': 0.5 * (x_range[0] + x_range[1]),
            'y': 0.5 * (y_range[0] + y_range[1]),
            'z': 0.5 * (z_range[0] + z_range[1]),
        }
        max_span = max(
            x_range[1] - x_range[0],
            y_range[1] - y_range[0],
            z_range[1] - z_range[0],
            1.0,
        )

        axis_x = scene_layout.get('xaxis', {})
        axis_y = scene_layout.get('yaxis', {})
        axis_z = scene_layout.get('zaxis', {})
        show_axes = not (
            axis_x.get('visible') is False
            and axis_y.get('visible') is False
            and axis_z.get('visible') is False
        )

        volume_layers = self._build_threejs_volume_layers()
        trace_keys = []
        legend_items = []
        for idx, trace in enumerate(self.initial_data):
            trace_json = _trace_to_plotly_json(trace)
            trace_key = f'trace-{idx}'
            trace_spec = self._threejs_trace_spec(trace_key, trace_json)
            if trace_spec is None:
                continue

            trace_keys.append((trace_key, trace_json.get('name')))
            if trace_spec.get('showlegend'):
                legend_items.append({
                    'key': trace_key,
                    'name': trace_spec.get('name') or trace_key,
                    'color': trace_spec.get('legend_color'),
                    'kind': 'trace',
                    'has_points': bool(trace_spec.get('points')),
                    'has_segments': bool(trace_spec.get('segments')),
                    'default_color': trace_spec.get('legend_color'),
                    'default_opacity': float(trace_spec.get('default_opacity', 1.0)),
                    'default_point_size': trace_spec.get('default_point_size'),
                })

        trace_key_by_name = {
            str(trace_name): trace_key
            for trace_key, trace_name in trace_keys
            if isinstance(trace_name, str) and trace_name
        }

        seen_volume_state_keys = set()
        for layer in volume_layers:
            state_key = str(layer.get('state_key') or layer.get('key'))
            if state_key in seen_volume_state_keys:
                continue
            seen_volume_state_keys.add(state_key)
            legend_items.append({
                'key': state_key,
                'name': str(layer.get('state_name') or layer.get('name') or state_key),
                'color': layer.get('legend_color'),
                'kind': 'volume',
            })

        group_visibility = {}
        for group_name, traces_list in self.trace_grouping_dict.items():
            group_visibility[group_name] = {}
            for trace_key, trace_name in trace_keys:
                visible = self.get_visibility(trace_name, traces_list)
                if visible == "legendonly":
                    group_visibility[group_name][trace_key] = "legendonly"
                else:
                    group_visibility[group_name][trace_key] = bool(visible)

        for group_name in group_visibility:
            for layer in volume_layers:
                state_key = str(layer.get('state_key') or layer.get('key'))
                group_visibility[group_name][state_key] = bool(layer.get('visible', True))

        frames_by_time = {}
        for time_val, frame in zip(self.time, frames):
            frames_by_time[round(float(time_val), 12)] = frame.to_plotly_json()

        frame_specs = []
        ordered_times = [float(t) for t in self._ordered_slider_times()]
        for time_val in ordered_times:
            frame_json = frames_by_time.get(round(float(time_val), 12))
            if frame_json is None:
                continue

            traces = []
            for idx, trace_json in enumerate(frame_json.get('data', [])):
                trace_spec = self._threejs_trace_spec(f'trace-{idx}', trace_json)
                if trace_spec is not None:
                    traces.append(trace_spec)

            frame_specs.append({
                'name': _format_time_label(time_val),
                'time': float(time_val),
                'traces': traces,
                'decorations': self._threejs_frame_decorations(
                    frame_json=frame_json,
                    time_value=float(time_val),
                    x_range=x_range,
                    y_range=y_range,
                    z_range=z_range,
                    fallback_center=center,
                    volume_layers=volume_layers,
                ),
            })

        default_group = list(self.trace_grouping_dict.keys())[0]
        initial_frame_index = 0
        for idx, frame_spec in enumerate(frame_specs):
            if np.isclose(float(frame_spec['time']), 0.0, atol=1e-9):
                initial_frame_index = idx
                break
        default_sky_catalog = {}
        if getattr(self, 'enable_sky_panel', False) and frame_specs:
            default_sky_catalog = _threejs_catalog_from_frame_spec(frame_specs[initial_frame_index])

        _annotate_threejs_point_motion_ranges(frame_specs)

        title_cfg = layout_json.get('title', {})
        if isinstance(title_cfg, dict):
            title_text = title_cfg.get('text') or ''
        else:
            title_text = str(title_cfg)

        return {
            'renderer': 'threejs',
            'title': title_text,
            'width': int(layout_json.get('width') or 900),
            'height': int(layout_json.get('height') or 700),
            'center': center,
            'max_span': float(max_span),
            'ranges': {
                'x': x_range,
                'y': y_range,
                'z': z_range,
            },
            'layout': layout_json,
            'axes': {
                'x': axis_x,
                'y': axis_y,
                'z': axis_z,
            },
            'theme': self._threejs_theme(layout_json),
            'frames': frame_specs,
            'initial_frame_index': int(initial_frame_index),
            'group_order': list(self.trace_grouping_dict.keys()),
            'default_group': default_group,
            'group_visibility': group_visibility,
            'legend': {
                'items': legend_items,
            },
            'show_axes': bool(show_axes),
            'playback_interval_ms': 500,
            'camera_up': {'x': 0.0, 'y': 0.0, 'z': 1.0},
            'sky_panel': self._build_threejs_sky_panel_spec(default_sky_catalog),
            'age_kde': self._build_threejs_age_kde_spec(trace_key_by_name),
            'cluster_filter': self._build_threejs_cluster_filter_spec(trace_key_by_name),
            'dendrogram': self._build_threejs_dendrogram_spec(trace_key_by_name),
            'volumes': {
                'enabled': bool(volume_layers),
                'layers': volume_layers,
            },
            'animation': {
                'fade_in_time_default': float(self.fade_in_time),
                'fade_in_and_out_default': bool(getattr(self, 'fade_in_and_out', False)),
                'focus_trace_key_default': trace_key_by_name.get(str(self.focus_group)) if self.focus_group else '',
                'focus_options': [
                    {
                        'key': trace_key,
                        'name': str(trace_name),
                    }
                    for trace_key, trace_name in trace_keys
                    if isinstance(trace_name, str)
                    and trace_name
                    and not str(trace_name).endswith(' Track')
                    and any(
                        trace.get('key') == trace_key and trace.get('points')
                        for frame_spec in frame_specs[:1]
                        for trace in frame_spec.get('traces', [])
                    )
                ],
            },
            'initial_state': {
                'click_selection_enabled': True,
            },
            'note': self._threejs_note_text(),
        }

    def _build_threejs_sky_panel_spec(self, default_catalog=None):
        if not getattr(self, 'enable_sky_panel', False):
            return {'enabled': False}
        merged_catalog = _merge_threejs_member_catalogs(default_catalog, self.sky_members_by_cluster)

        return {
            'enabled': True,
            'radius_deg': float(self.sky_radius_deg),
            'frame': self.sky_frame,
            'survey': self.sky_survey,
            'members_by_cluster': merged_catalog,
        }

    def _build_threejs_age_kde_spec(self, trace_key_by_name=None):
        if not getattr(self, 'show_age_kde_inset', False):
            return {'enabled': False}

        x_grid = np.asarray(getattr(self, 'kde_x_grid', np.array([], dtype=float)), dtype=float)
        if x_grid.size == 0:
            return {'enabled': False}

        trace_key_by_name = trace_key_by_name or {}
        traces = []
        for trace_name in getattr(self, 'kde_trace_order', []):
            density = np.asarray(
                self.kde_density_by_trace.get(trace_name, np.array([], dtype=float)),
                dtype=float,
            )
            if density.size != x_grid.size:
                continue
            traces.append({
                'trace_name': str(trace_name),
                'trace_key': trace_key_by_name.get(str(trace_name)),
                'color': str(self.kde_color_by_trace.get(trace_name, self.kde_axis_color)),
                'opacity': float(self.kde_opacity_by_trace.get(trace_name, 1.0)),
                'x': [float(value) for value in x_grid],
                'y': [float(value) for value in density],
            })

        cluster_points = []
        for cluster_group in self.data_collection.get_all_clusters():
            trace_name = str(getattr(cluster_group, 'data_name', '') or '')
            data_df = (
                cluster_group.df_int
                if getattr(cluster_group, 'df_int', None) is not None
                else getattr(cluster_group, 'df', None)
            )
            if data_df is None or ('age_myr' not in data_df.columns):
                continue

            df_present = data_df
            if 'time' in data_df.columns:
                time_values = pd.to_numeric(data_df['time'], errors='coerce').to_numpy(dtype=float)
                zero_mask = np.isclose(time_values, 0.0, rtol=0.0, atol=1e-9)
                if np.any(zero_mask):
                    df_present = data_df.loc[zero_mask]

            age_values = pd.to_numeric(df_present['age_myr'], errors='coerce').to_numpy(dtype=float)
            if 'name' in df_present.columns:
                cluster_names = df_present['name'].astype(str).to_numpy(dtype=object)
            else:
                cluster_names = np.array([trace_name] * len(df_present), dtype=object)

            trace_key = trace_key_by_name.get(trace_name)
            for cluster_name, age_now in zip(cluster_names, age_values):
                if not np.isfinite(age_now):
                    continue
                cluster_points.append({
                    'cluster_name': str(cluster_name),
                    'trace_name': trace_name,
                    'trace_key': trace_key,
                    'age_now_myr': float(age_now),
                })

        return {
            'enabled': True,
            'title': 'Relative SFH',
            'x_range': [
                float(value)
                for value in getattr(
                    self,
                    'kde_x_range',
                    [float(np.min(x_grid)), float(np.max(x_grid))],
                )
            ],
            'y_range': [0.0, float(getattr(self, 'kde_y_max', 1.0))],
            'axis_color': str(
                getattr(
                    self,
                    'kde_axis_color',
                    'white' if self.figure_theme == 'dark' else 'black',
                )
            ),
            'bandwidth_myr': float(self.age_kde_bandwidth_myr),
            'traces': traces,
            'cluster_points': cluster_points,
        }

    def _build_threejs_cluster_filter_spec(self, trace_key_by_name=None):
        trace_key_by_name = trace_key_by_name or {}
        entries = []
        seen_keys = set()

        for cluster_group in self.data_collection.get_all_clusters():
            trace_name = str(getattr(cluster_group, 'data_name', '') or '')
            trace_key = trace_key_by_name.get(trace_name)
            data_df = (
                cluster_group.df_int
                if getattr(cluster_group, 'df_int', None) is not None
                else getattr(cluster_group, 'df', None)
            )
            if data_df is None or data_df.empty or ('age_myr' not in data_df.columns):
                continue

            df_present = data_df
            if 'time' in data_df.columns:
                time_values = pd.to_numeric(data_df['time'], errors='coerce').to_numpy(dtype=float)
                zero_mask = np.isclose(time_values, 0.0, rtol=0.0, atol=1e-9)
                if np.any(zero_mask):
                    df_present = data_df.loc[zero_mask]

            if df_present.empty:
                continue

            age_values = pd.to_numeric(df_present['age_myr'], errors='coerce').to_numpy(dtype=float)
            if 'n_stars' in df_present.columns:
                n_stars_values = pd.to_numeric(df_present['n_stars'], errors='coerce').to_numpy(dtype=float)
            else:
                n_stars_values = np.full(len(df_present), np.nan, dtype=float)

            if 'name' in df_present.columns:
                cluster_names = df_present['name'].astype(str).to_numpy(dtype=object)
            else:
                cluster_names = np.array([trace_name] * len(df_present), dtype=object)

            for cluster_name, age_now, n_stars in zip(cluster_names, age_values, n_stars_values):
                selection_key = _selection_identity_key({
                    'cluster_name': str(cluster_name),
                    'trace_name': trace_name,
                })
                if not selection_key or selection_key in seen_keys:
                    continue
                seen_keys.add(selection_key)
                entries.append({
                    'selection_key': str(selection_key),
                    'cluster_name': str(cluster_name),
                    'trace_name': trace_name,
                    'trace_key': trace_key,
                    'age_now_myr': float(age_now) if np.isfinite(age_now) else np.nan,
                    'n_stars': float(n_stars) if np.isfinite(n_stars) else np.nan,
                })

        parameters = []
        for key, label, unit in (
            ('age_now_myr', 'Age', 'Myr'),
            ('n_stars', 'Stars', ''),
        ):
            values = np.asarray([entry.get(key, np.nan) for entry in entries], dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            parameters.append({
                'key': key,
                'label': label,
                'unit': unit,
                'min': float(np.nanmin(values)),
                'max': float(np.nanmax(values)),
            })

        return {
            'enabled': bool(entries and parameters),
            'default_parameter_key': parameters[0]['key'] if parameters else '',
            'parameters': parameters,
            'entries': entries,
        }

    def _build_threejs_dendrogram_spec(self, trace_key_by_name=None):
        trace_key_by_name = trace_key_by_name or {}
        entries = []
        trace_options = []

        if getattr(self, 'coord_system', 'centered') == 'rot':
            x_col, y_col, z_col = 'x_rot', 'y_rot', 'z_rot'
        else:
            x_col, y_col, z_col = 'x', 'y', 'z'

        for cluster_group in self.data_collection.get_all_clusters():
            trace_name = str(getattr(cluster_group, 'data_name', '') or '')
            trace_key = trace_key_by_name.get(trace_name)
            if not trace_name or not trace_key:
                continue

            data_df = (
                cluster_group.df_int
                if getattr(cluster_group, 'df_int', None) is not None
                else getattr(cluster_group, 'df', None)
            )
            if data_df is None or data_df.empty or ('age_myr' not in data_df.columns):
                continue

            required_cols = {x_col, y_col, z_col}
            if data_df.empty or not required_cols.issubset(data_df.columns):
                continue

            cluster_color = getattr(cluster_group, 'color', None)
            if not isinstance(cluster_color, str) or not cluster_color:
                cluster_color = '#ffffff'

            working_df = data_df.copy()
            if 'name' in working_df.columns:
                working_df['__cluster_name'] = working_df['name'].astype(str)
            else:
                working_df['__cluster_name'] = trace_name

            if 'time' in working_df.columns:
                working_df['__time_myr'] = pd.to_numeric(working_df['time'], errors='coerce')
            else:
                working_df['__time_myr'] = 0.0
            working_df['__age_now_myr'] = pd.to_numeric(working_df['age_myr'], errors='coerce')
            working_df['__x'] = pd.to_numeric(working_df[x_col], errors='coerce')
            working_df['__y'] = pd.to_numeric(working_df[y_col], errors='coerce')
            working_df['__z'] = pd.to_numeric(working_df[z_col], errors='coerce')

            trace_entry_count = 0
            trace_age_values = pd.to_numeric(working_df['__age_now_myr'], errors='coerce').to_numpy(dtype=float)
            trace_age_values = trace_age_values[np.isfinite(trace_age_values)]
            trace_max_age = float(np.nanmax(trace_age_values)) if trace_age_values.size else 0.0

            for cluster_name, cluster_df in working_df.groupby('__cluster_name', sort=False):
                age_values = pd.to_numeric(cluster_df['__age_now_myr'], errors='coerce').to_numpy(dtype=float)
                age_values = age_values[np.isfinite(age_values)]
                if age_values.size == 0:
                    continue
                age_now = float(np.nanmax(age_values))
                birth_time_myr = -age_now

                time_samples = pd.to_numeric(cluster_df['__time_myr'], errors='coerce').to_numpy(dtype=float)
                x_samples = pd.to_numeric(cluster_df['__x'], errors='coerce').to_numpy(dtype=float)
                y_samples = pd.to_numeric(cluster_df['__y'], errors='coerce').to_numpy(dtype=float)
                z_samples = pd.to_numeric(cluster_df['__z'], errors='coerce').to_numpy(dtype=float)
                valid_mask = (
                    np.isfinite(time_samples)
                    & np.isfinite(x_samples)
                    & np.isfinite(y_samples)
                    & np.isfinite(z_samples)
                )
                if not np.any(valid_mask):
                    continue

                time_samples = time_samples[valid_mask]
                x_samples = x_samples[valid_mask]
                y_samples = y_samples[valid_mask]
                z_samples = z_samples[valid_mask]

                order = np.argsort(time_samples)
                time_samples = time_samples[order]
                x_samples = x_samples[order]
                y_samples = y_samples[order]
                z_samples = z_samples[order]

                unique_times, unique_indices = np.unique(time_samples, return_index=True)
                time_samples = unique_times.astype(float)
                x_samples = x_samples[unique_indices].astype(float)
                y_samples = y_samples[unique_indices].astype(float)
                z_samples = z_samples[unique_indices].astype(float)

                if time_samples.size == 0:
                    continue

                birth_time_for_interp = float(np.clip(birth_time_myr, np.nanmin(time_samples), np.nanmax(time_samples)))
                birth_x = float(np.interp(birth_time_for_interp, time_samples, x_samples))
                birth_y = float(np.interp(birth_time_for_interp, time_samples, y_samples))
                birth_z = float(np.interp(birth_time_for_interp, time_samples, z_samples))
                if not (np.isfinite(birth_x) and np.isfinite(birth_y) and np.isfinite(birth_z)):
                    continue

                selection_key = _selection_identity_key({
                    'cluster_name': str(cluster_name),
                    'trace_name': trace_name,
                })
                if not selection_key:
                    continue
                entries.append({
                    'selection_key': str(selection_key),
                    'cluster_name': str(cluster_name),
                    'trace_name': trace_name,
                    'trace_key': str(trace_key),
                    'color': str(cluster_color),
                    'age_now_myr': float(age_now),
                    'birth_time_myr': float(birth_time_myr),
                    'x_birth': birth_x,
                    'y_birth': birth_y,
                    'z_birth': birth_z,
                    'time_samples': [float(value) for value in time_samples.tolist()],
                    'x_samples': [float(value) for value in x_samples.tolist()],
                    'y_samples': [float(value) for value in y_samples.tolist()],
                    'z_samples': [float(value) for value in z_samples.tolist()],
                })
                trace_entry_count += 1

            if trace_entry_count:
                trace_options.append({
                    'trace_name': trace_name,
                    'trace_key': str(trace_key),
                    'color': str(cluster_color),
                    'count': int(trace_entry_count),
                    'max_age_myr': float(trace_max_age),
                })

        age_values = np.asarray([entry['age_now_myr'] for entry in entries], dtype=float)
        max_age = float(np.nanmax(age_values)) if age_values.size else 0.0
        return {
            'enabled': bool(entries and trace_options),
            'title': 'Birth Tree',
            'default_trace_key': trace_options[0]['trace_key'] if trace_options else '',
            'default_threshold_pc': 100.0,
            'threshold_min_pc': 0.0,
            'threshold_max_pc': max(5000.0, max_age * 50.0 if np.isfinite(max_age) else 5000.0),
            'max_age_myr': max_age,
            'traces': trace_options,
            'entries': entries,
        }

    def _build_threejs_volume_layers(self):
        if not getattr(self, 'volume_configs', None):
            return []

        zero_matches = np.flatnonzero(np.isclose(np.asarray(self.time, dtype=float), 0.0, rtol=0.0, atol=1e-9))
        if zero_matches.size:
            zero_idx = int(zero_matches[0])
        else:
            zero_idx = 0

        x_rf = float(np.asarray(self.coords_center_int[0], dtype=float)[zero_idx])
        y_rf = float(np.asarray(self.coords_center_int[1], dtype=float)[zero_idx])
        z_rf = float(np.asarray(self.coords_center_int[2], dtype=float)[zero_idx])
        center_offset = {'x': x_rf, 'y': y_rf, 'z': z_rf}

        layers = []
        for idx, volume_cfg in enumerate(self.volume_configs):
            layer = _build_threejs_volume_layer_spec(
                volume_cfg,
                center_offset=center_offset,
                index=idx,
                include_sky_overlay=bool(getattr(self, 'enable_sky_panel', False)),
            )
            if layer is not None:
                layers.append(layer)
        return layers

    def _threejs_note_text(self):
        note = 'Three.js modules are loaded from a CDN in this renderer path.'
        if getattr(self, 'enable_sky_panel', False):
            note += ' Use the widget menu to open Aladin Lite.'
        if getattr(self, 'show_age_kde_inset', False):
            note += ' The widget menu also opens the age KDE view.'
        if getattr(self, 'volume_configs', None):
            if any(
                _coerce_threejs_volume_time_myr((cfg or {}).get('time_myr')) is not None
                for cfg in (self.volume_configs or [])
                if isinstance(cfg, dict)
            ):
                note += ' Time-resolved volumetric layers track the current animation frame and are controlled from the left toolbar.'
            else:
                note += ' Volumetric layers render at t=0 and are controlled from the left toolbar.'
        return note

    def _threejs_frame_decorations(
        self,
        frame_json,
        time_value,
        x_range,
        y_range,
        z_range,
        fallback_center,
        volume_layers=None,
    ):
        """Return frame-local decorative objects for the Three.js renderer."""
        decorations = []
        for layer in volume_layers or []:
            layer_time = _coerce_threejs_volume_time_myr(layer.get('time_myr'))
            if layer_time is not None:
                if not np.isclose(float(time_value), float(layer_time), atol=1e-9):
                    continue
            elif layer.get('only_at_t0', True) and (not np.isclose(float(time_value), 0.0, atol=1e-9)):
                continue
            decorations.append({
                'kind': 'volume_layer',
                'key': layer.get('key'),
                'state_key': layer.get('state_key') or layer.get('key'),
            })

        if (not getattr(self, 'show_milky_way_model', False)) or (not np.isclose(float(time_value), 0.0, atol=1e-9)):
            return decorations

        x_span = max(float(x_range[1]) - float(x_range[0]), 1.0)
        y_span = max(float(y_range[1]) - float(y_range[0]), 1.0)
        z_span = max(float(z_range[1]) - float(z_range[0]), 1.0)
        disc_radius_pc = float(np.clip(0.49 * min(x_span, y_span), 4500.0, 14000.0))
        center = self._threejs_milky_way_center(frame_json, fallback_center)

        decorations.append({
            'kind': 'milky_way_model',
            'center': center,
            'disc_radius_pc': disc_radius_pc,
            'disc_thickness_pc': max(180.0, 0.018 * disc_radius_pc),
            'bulge_radius_pc': 0.18 * disc_radius_pc,
            'halo_radius_pc': 1.16 * disc_radius_pc,
            'bulge_height_pc': min(max(0.75 * z_span, 700.0), 2400.0),
            'dust_inner_radius_pc': 0.14 * disc_radius_pc,
            'dust_outer_radius_pc': 0.92 * disc_radius_pc,
        })
        return decorations

    def _threejs_milky_way_center(self, frame_json, fallback_center):
        """Prefer the plotted Galactic-center marker when anchoring the Milky Way model."""
        for trace_json in frame_json.get('data', []):
            if trace_json.get('type', 'scatter3d') != 'scatter3d':
                continue
            if trace_json.get('name') != 'GC':
                continue
            x_vals = _as_object_list(trace_json.get('x'))
            y_vals = _as_object_list(trace_json.get('y'))
            z_vals = _as_object_list(trace_json.get('z'))
            if not x_vals or not y_vals or not z_vals:
                continue
            try:
                return {
                    'x': float(x_vals[0]),
                    'y': float(y_vals[0]),
                    'z': float(z_vals[0]),
                }
            except Exception:
                continue
        return {
            'x': float(fallback_center['x']),
            'y': float(fallback_center['y']),
            'z': float(fallback_center['z']),
        }

    def _threejs_theme(self, layout_json):
        """Translate the Plotly layout theme into a smaller scene theme."""
        scene_layout = layout_json.get('scene', {})
        axis_x = scene_layout.get('xaxis', {})
        axis_color = (
            axis_x.get('linecolor')
            or axis_x.get('tickfont', {}).get('color')
            or axis_x.get('title_font', {}).get('color')
            or ('gray' if self.figure_theme == 'dark' else 'black')
        )
        text_color = (
            layout_json.get('legend', {}).get('font', {}).get('color')
            or axis_x.get('tickfont', {}).get('color')
            or axis_x.get('title_font', {}).get('color')
            or ('white' if self.figure_theme == 'dark' else 'black')
        )
        if self.figure_theme == 'dark':
            panel_bg = 'rgba(0,0,0,0.48)'
            panel_border = 'rgba(110,110,110,0.60)'
            panel_solid = '#121212'
            footprint = '#6ec5ff'
        else:
            panel_bg = 'rgba(255,255,255,0.78)'
            panel_border = 'rgba(80,80,80,0.35)'
            panel_solid = '#f7f7f7'
            footprint = '#0b67b1'

        return {
            'paper_bgcolor': layout_json.get('paper_bgcolor', 'black'),
            'scene_bgcolor': scene_layout.get('bgcolor', layout_json.get('paper_bgcolor', 'black')),
            'text_color': text_color,
            'axis_color': axis_color,
            'panel_bg': panel_bg,
            'panel_border': panel_border,
            'panel_solid': panel_solid,
            'footprint': footprint,
        }

    def _threejs_trace_spec(self, trace_key, trace_json):
        """Convert a Plotly trace into a simpler primitive spec for three.js."""
        if trace_json.get('type', 'scatter3d') != 'scatter3d':
            return None

        mode = str(trace_json.get('mode', 'markers'))
        spec = {
            'key': trace_key,
            'name': trace_json.get('name') or trace_key,
            'showlegend': bool(trace_json.get('showlegend', True)),
        }
        legend_color = None

        if 'lines' in mode:
            segments = _line_segments_from_trace(trace_json)
            if segments:
                line_color, line_opacity = _color_to_css_and_opacity(
                    trace_json.get('line', {}).get('color'),
                    base_opacity=float(trace_json.get('opacity', 1.0)),
                )
                spec['segments'] = segments
                spec['line'] = {
                    'color': line_color,
                    'width': float(trace_json.get('line', {}).get('width', 1.0) or 1.0),
                    'dash': trace_json.get('line', {}).get('dash', 'solid'),
                }
                spec['opacity'] = line_opacity
                legend_color = line_color

        if 'markers' in mode:
            points = _points_from_trace(
                trace_json,
                default_opacity=float(trace_json.get('opacity', 1.0)),
                include_selection=bool(getattr(self, 'enable_sky_panel', False)),
            )
            if points:
                spec['points'] = points
                point_sizes = [
                    float(point.get('size'))
                    for point in points
                    if np.isfinite(point.get('size')) and float(point.get('size')) > 0.0
                ]
                if point_sizes:
                    spec['default_point_size'] = float(np.median(point_sizes))
                point_opacities = [
                    float(point.get('opacity'))
                    for point in points
                    if np.isfinite(point.get('opacity'))
                ]
                if point_opacities:
                    spec['default_opacity'] = float(np.median(point_opacities))
                if legend_color is None:
                    legend_color = points[0].get('color')

        if 'text' in mode:
            labels = _labels_from_trace(trace_json)
            if labels:
                spec['labels'] = labels
                if legend_color is None:
                    legend_color = labels[0].get('color')

        if legend_color is not None:
            spec['legend_color'] = legend_color

        if 'default_opacity' not in spec:
            spec['default_opacity'] = float(spec.get('opacity', 1.0))

        if any(key in spec for key in ('segments', 'points', 'labels')):
            return spec
        return None



####################################################################################################
def _normalize_renderer_name(renderer):
    renderer_name = str(renderer).strip().lower()
    if renderer_name in ('plotly',):
        return 'plotly'
    if renderer_name in ('three', 'threejs', 'three.js'):
        return 'threejs'
    raise ValueError("renderer must be one of: 'plotly', 'threejs'.")


def _trace_to_plotly_json(trace):
    if isinstance(trace, dict):
        return copy.deepcopy(trace)
    if hasattr(trace, 'to_plotly_json'):
        return trace.to_plotly_json()
    raise TypeError(f'Unsupported trace type for serialization: {type(trace)!r}')


def _coerce_range(values, default):
    if not isinstance(values, (list, tuple, np.ndarray)) or len(values) != 2:
        return [float(default[0]), float(default[1])]
    try:
        return [float(values[0]), float(values[1])]
    except Exception:
        return [float(default[0]), float(default[1])]


def _format_time_label(time_value):
    rounded = round(float(time_value), 10)
    if np.isclose(rounded, round(rounded), atol=1e-9):
        return str(int(round(rounded)))
    return f'{rounded:.1f}'.rstrip('0').rstrip('.')


def _is_sequence_value(value):
    return isinstance(value, (list, tuple, np.ndarray, pd.Series))


def _expand_value(value, length, default=None):
    if length <= 0:
        return []
    if value is None:
        return [default] * length
    if _is_sequence_value(value):
        values = _as_object_list(value)
        if len(values) == length:
            return values
        if len(values) == 1:
            return values * length
        if len(values) < length:
            values.extend([values[-1]] * (length - len(values)))
            return values
        return values[:length]
    return [value] * length


def _coerce_float(value, default=0.0):
    try:
        out = float(value)
        if np.isfinite(out):
            return out
    except Exception:
        pass
    return float(default)


def _as_object_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, pd.Series):
        return value.tolist()
    if isinstance(value, np.ndarray):
        return np.asarray(value, dtype=object).reshape(-1).tolist()
    return [value]


def _color_to_css_and_opacity(color, base_opacity=1.0):
    opacity = _coerce_float(base_opacity, 1.0)
    if color is None:
        return '#808080', opacity

    if isinstance(color, (tuple, list)):
        if len(color) == 4:
            r, g, b, a = color
            return (
                f'rgb({int(r)}, {int(g)}, {int(b)})',
                opacity * _coerce_float(a, 1.0),
            )
        if len(color) == 3:
            r, g, b = color
            return f'rgb({int(r)}, {int(g)}, {int(b)})', opacity

    color_text = str(color).strip()
    if not color_text:
        return '#808080', opacity

    lower = color_text.lower()
    if lower.startswith('rgba(') and lower.endswith(')'):
        inner = color_text[color_text.find('(') + 1: color_text.rfind(')')]
        parts = [part.strip() for part in inner.split(',')]
        if len(parts) == 4:
            try:
                r, g, b = [int(float(part)) for part in parts[:3]]
                a = float(parts[3])
                return f'rgb({r}, {g}, {b})', opacity * a
            except Exception:
                return color_text, opacity

    if lower.startswith('rgb(') and lower.endswith(')'):
        return color_text, opacity

    if color_text.startswith('#'):
        return color_text, opacity

    try:
        return webcolors.name_to_hex(color_text), opacity
    except Exception:
        return color_text, opacity


def _resolve_marker_color_values(marker, length, default_opacity=1.0):
    color_value = marker.get('color')
    if color_value is None:
        css, alpha = _color_to_css_and_opacity('#808080', default_opacity)
        return [css] * length, [alpha] * length

    if _is_sequence_value(color_value):
        color_values = list(np.asarray(color_value, dtype=object).tolist())
        if color_values and all(_looks_numeric(val) for val in color_values):
            colorscale = marker.get('colorscale', 'Viridis')
            cmin = _coerce_float(marker.get('cmin'), np.nan)
            cmax = _coerce_float(marker.get('cmax'), np.nan)
            numeric = np.asarray([float(val) for val in color_values], dtype=float)
            if not np.isfinite(cmin):
                cmin = float(np.nanmin(numeric)) if numeric.size else 0.0
            if not np.isfinite(cmax):
                cmax = float(np.nanmax(numeric)) if numeric.size else 1.0
            if np.isclose(cmax, cmin, atol=1e-12):
                scaled = np.zeros_like(numeric)
            else:
                scaled = np.clip((numeric - cmin) / (cmax - cmin), 0.0, 1.0)
            sampled = pc.sample_colorscale(colorscale, scaled.tolist(), colortype='rgb')
            css_values = []
            alpha_values = []
            for item in sampled:
                css, alpha = _color_to_css_and_opacity(item, default_opacity)
                css_values.append(css)
                alpha_values.append(alpha)
            return css_values, alpha_values

        expanded = _expand_value(color_values, length, '#808080')
        css_values = []
        alpha_values = []
        for item in expanded:
            css, alpha = _color_to_css_and_opacity(item, default_opacity)
            css_values.append(css)
            alpha_values.append(alpha)
        return css_values, alpha_values

    css, alpha = _color_to_css_and_opacity(color_value, default_opacity)
    return [css] * length, [alpha] * length


def _looks_numeric(value):
    try:
        return np.isfinite(float(value))
    except Exception:
        return False


def _line_segments_from_trace(trace_json):
    x_vals = _as_object_list(trace_json.get('x'))
    y_vals = _as_object_list(trace_json.get('y'))
    z_vals = _as_object_list(trace_json.get('z'))
    segments = []
    prev = None
    for x_val, y_val, z_val in zip(x_vals, y_vals, z_vals):
        if x_val is None or y_val is None or z_val is None:
            prev = None
            continue
        try:
            point = [float(x_val), float(y_val), float(z_val)]
        except Exception:
            prev = None
            continue
        if not all(np.isfinite(point)):
            prev = None
            continue
        if prev is not None:
            segments.append(prev + point)
        prev = point
    return segments


def _selection_from_customdata_row(row):
    if row is None:
        return None

    values = _as_object_list(row)
    if len(values) <= CUSTOMDATA_IDX_Z0:
        return None

    l_deg = _coerce_float(values[CUSTOMDATA_IDX_L0_DEG], np.nan)
    b_deg = _coerce_float(values[CUSTOMDATA_IDX_B0_DEG], np.nan)
    dist_pc = _coerce_float(values[CUSTOMDATA_IDX_DIST0_PC], np.nan)
    x0 = _coerce_float(values[CUSTOMDATA_IDX_X0], np.nan)
    y0 = _coerce_float(values[CUSTOMDATA_IDX_Y0], np.nan)
    z0 = _coerce_float(values[CUSTOMDATA_IDX_Z0], np.nan)
    if (not np.isfinite(l_deg)) or (not np.isfinite(b_deg)):
        return None

    age_now = _coerce_float(values[CUSTOMDATA_IDX_AGE_NOW], np.nan)
    age_at_t = _coerce_float(values[CUSTOMDATA_IDX_AGE_AT_T], np.nan)
    click_time_myr = np.nan
    if np.isfinite(age_now) and np.isfinite(age_at_t):
        click_time_myr = float(age_at_t - age_now)
    ra_deg = np.nan
    dec_deg = np.nan
    try:
        icrs = SkyCoord(l=[l_deg] * u.deg, b=[b_deg] * u.deg, frame='galactic').icrs
        ra_deg = float(icrs.ra.deg[0])
        dec_deg = float(icrs.dec.deg[0])
    except Exception:
        pass

    cluster_name = None
    if len(values) > CUSTOMDATA_IDX_CLUSTER_NAME and values[CUSTOMDATA_IDX_CLUSTER_NAME] not in (None, ''):
        cluster_name = str(values[CUSTOMDATA_IDX_CLUSTER_NAME])

    cluster_color = None
    if len(values) > CUSTOMDATA_IDX_CLUSTER_COLOR and values[CUSTOMDATA_IDX_CLUSTER_COLOR] not in (None, ''):
        cluster_color = str(values[CUSTOMDATA_IDX_CLUSTER_COLOR])

    return {
        'l_deg': float(l_deg),
        'b_deg': float(b_deg),
        'dist_pc': float(dist_pc) if np.isfinite(dist_pc) else np.nan,
        'x0': float(x0) if np.isfinite(x0) else np.nan,
        'y0': float(y0) if np.isfinite(y0) else np.nan,
        'z0': float(z0) if np.isfinite(z0) else np.nan,
        'age_now_myr': float(age_now) if np.isfinite(age_now) else np.nan,
        'age_at_t_myr': float(age_at_t) if np.isfinite(age_at_t) else np.nan,
        'click_time_myr': float(click_time_myr) if np.isfinite(click_time_myr) else np.nan,
        'ra_deg': float(ra_deg) if np.isfinite(ra_deg) else np.nan,
        'dec_deg': float(dec_deg) if np.isfinite(dec_deg) else np.nan,
        'cluster_name': cluster_name,
        'cluster_color': cluster_color,
    }


def _selection_identity_key(selection):
    if not isinstance(selection, dict):
        return ''

    cluster_name = str(selection.get('cluster_name') or '').strip()
    if cluster_name:
        return cluster_name

    x0 = _coerce_float(selection.get('x0'), np.nan)
    y0 = _coerce_float(selection.get('y0'), np.nan)
    z0 = _coerce_float(selection.get('z0'), np.nan)
    if np.isfinite(x0) and np.isfinite(y0) and np.isfinite(z0):
        return f'{x0:.6f}|{y0:.6f}|{z0:.6f}'

    ra_deg = _coerce_float(selection.get('ra_deg'), np.nan)
    dec_deg = _coerce_float(selection.get('dec_deg'), np.nan)
    if np.isfinite(ra_deg) and np.isfinite(dec_deg):
        return f'{ra_deg:.6f}|{dec_deg:.6f}'

    return str(selection.get('trace_name') or '').strip()


def _motion_from_selection(selection):
    if not isinstance(selection, dict):
        return None

    age_now_myr = _coerce_float(selection.get('age_now_myr'), np.nan)
    age_at_t_myr = _coerce_float(selection.get('age_at_t_myr'), np.nan)
    if (not np.isfinite(age_now_myr)) or (not np.isfinite(age_at_t_myr)):
        return None

    return {
        'key': _selection_identity_key(selection),
        'cluster_name': str(selection.get('cluster_name') or '').strip(),
        'trace_name': str(selection.get('trace_name') or '').strip(),
        'age_now_myr': float(age_now_myr),
        'age_at_t_myr': float(age_at_t_myr),
        'time_myr': float(age_at_t_myr - age_now_myr),
    }


def _points_from_trace(trace_json, default_opacity=1.0, include_selection=False):
    x_vals = _as_object_list(trace_json.get('x'))
    y_vals = _as_object_list(trace_json.get('y'))
    z_vals = _as_object_list(trace_json.get('z'))
    n_points = min(len(x_vals), len(y_vals), len(z_vals))
    if n_points == 0:
        return []

    marker = trace_json.get('marker', {})
    sizes = [_coerce_float(val, 4.0) for val in _expand_value(marker.get('size', 4.0), n_points, 4.0)]
    symbols = [str(val) for val in _expand_value(marker.get('symbol', 'circle'), n_points, 'circle')]
    hovertext = _expand_value(trace_json.get('hovertext'), n_points, '')

    base_marker_opacity = marker.get('opacity')
    if base_marker_opacity is None:
        marker_opacity = [float(default_opacity)] * n_points
    else:
        marker_opacity = [
            _coerce_float(val, default_opacity) * float(default_opacity)
            for val in _expand_value(base_marker_opacity, n_points, default_opacity)
        ]

    colors, color_opacity = _resolve_marker_color_values(marker, n_points, default_opacity=1.0)
    customdata = trace_json.get('customdata')
    custom_rows = None
    if customdata is not None:
        custom_arr = np.asarray(customdata, dtype=object)
        if custom_arr.ndim == 1:
            custom_arr = custom_arr.reshape(-1, 1)
        if custom_arr.ndim == 2 and custom_arr.shape[0] >= n_points:
            custom_rows = [custom_arr[idx, :].tolist() for idx in range(n_points)]

    points = []
    for idx in range(n_points):
        try:
            x_val = float(x_vals[idx])
            y_val = float(y_vals[idx])
            z_val = float(z_vals[idx])
        except Exception:
            continue
        if not np.isfinite(x_val) or not np.isfinite(y_val) or not np.isfinite(z_val):
            continue

        points.append({
            'x': x_val,
            'y': y_val,
            'z': z_val,
            'size': float(max(sizes[idx], 0.0)),
            'symbol': symbols[idx],
            'color': colors[idx],
            'opacity': float(np.clip(marker_opacity[idx] * color_opacity[idx], 0.0, 1.0)),
            'hovertext': str(hovertext[idx]) if hovertext[idx] is not None else '',
        })
        if custom_rows is not None:
            selection = _selection_from_customdata_row(custom_rows[idx])
            if selection is not None:
                selection['trace_name'] = trace_json.get('name')
                motion = _motion_from_selection(selection)
                if motion is not None:
                    points[-1]['motion'] = motion
            if include_selection:
                points[-1]['selection'] = selection

    return points


def _annotate_threejs_point_motion_ranges(frame_specs):
    motion_ranges = {}
    for frame_spec in frame_specs:
        for trace in frame_spec.get('traces', []):
            for point in trace.get('points', []):
                motion = point.get('motion')
                if not isinstance(motion, dict):
                    continue
                motion_key = str(motion.get('key') or '').strip()
                if not motion_key:
                    continue
                entry = motion_ranges.setdefault(
                    motion_key,
                    {
                        'size_min': np.inf,
                        'size_max': -np.inf,
                        'opacity_min': np.inf,
                        'opacity_max': -np.inf,
                    },
                )
                point_size = _coerce_float(point.get('size'), np.nan)
                point_opacity = _coerce_float(point.get('opacity'), np.nan)
                if np.isfinite(point_size):
                    entry['size_min'] = min(entry['size_min'], point_size)
                    entry['size_max'] = max(entry['size_max'], point_size)
                if np.isfinite(point_opacity):
                    entry['opacity_min'] = min(entry['opacity_min'], point_opacity)
                    entry['opacity_max'] = max(entry['opacity_max'], point_opacity)

    for frame_spec in frame_specs:
        for trace in frame_spec.get('traces', []):
            for point in trace.get('points', []):
                motion = point.get('motion')
                if not isinstance(motion, dict):
                    continue
                motion_key = str(motion.get('key') or '').strip()
                ranges = motion_ranges.get(motion_key)
                if not ranges:
                    continue
                size_min = ranges['size_min']
                size_max = ranges['size_max']
                opacity_min = ranges['opacity_min']
                opacity_max = ranges['opacity_max']
                motion['size_min'] = float(size_min) if np.isfinite(size_min) else float(_coerce_float(point.get('size'), 0.0))
                motion['size_max'] = float(size_max) if np.isfinite(size_max) else float(_coerce_float(point.get('size'), 0.0))
                motion['opacity_min'] = float(opacity_min) if np.isfinite(opacity_min) else float(_coerce_float(point.get('opacity'), 1.0))
                motion['opacity_max'] = float(opacity_max) if np.isfinite(opacity_max) else float(_coerce_float(point.get('opacity'), 1.0))


def _labels_from_trace(trace_json):
    x_vals = _as_object_list(trace_json.get('x'))
    y_vals = _as_object_list(trace_json.get('y'))
    z_vals = _as_object_list(trace_json.get('z'))
    n_labels = min(len(x_vals), len(y_vals), len(z_vals))
    if n_labels == 0:
        return []

    texts = _expand_value(trace_json.get('text'), n_labels, '')
    textfont = trace_json.get('textfont', {})
    color, _ = _color_to_css_and_opacity(textfont.get('color'), 1.0)
    font_size = _coerce_float(textfont.get('size'), 12.0)
    font_family = textfont.get('family', 'helvetica')
    trace_name = str(trace_json.get('name') or '')
    is_gc_label = trace_name == 'GC'
    is_galactic_radius_label = trace_name.startswith('R=') and trace_name.endswith(' Label')

    labels = []
    for idx in range(n_labels):
        text = texts[idx]
        if text in (None, ''):
            continue
        try:
            x_val = float(x_vals[idx])
            y_val = float(y_vals[idx])
            z_val = float(z_vals[idx])
        except Exception:
            continue
        if not np.isfinite(x_val) or not np.isfinite(y_val) or not np.isfinite(z_val):
            continue

        labels.append({
            'text': str(text),
            'x': x_val,
            'y': y_val,
            'z': z_val,
            'color': color,
            'size': font_size,
            'family': font_family,
            'screen_stable': bool(is_gc_label or is_galactic_radius_label),
            'screen_px': 10.0 if is_gc_label else (9.0 if is_galactic_radius_label else None),
        })

    return labels


def _load_threejs_cluster_catalog(cluster_members_file):
    if not cluster_members_file:
        return None

    try:
        df = pd.read_csv(cluster_members_file, usecols=['name', 'l', 'b'])
    except Exception:
        return None

    if df.empty:
        return None

    df = df.rename(columns={'name': 'cluster_name', 'l': 'l_deg', 'b': 'b_deg'})
    df['cluster_name'] = df['cluster_name'].astype(str)
    df['l_deg'] = pd.to_numeric(df['l_deg'], errors='coerce')
    df['b_deg'] = pd.to_numeric(df['b_deg'], errors='coerce')
    df = df.loc[df['l_deg'].notnull() & df['b_deg'].notnull()].copy()
    if df.empty:
        return None

    grouped = {}
    for cluster_name, grp in df.groupby('cluster_name', sort=False):
        l_vals = grp['l_deg'].to_numpy(dtype=np.float64)
        b_vals = grp['b_deg'].to_numpy(dtype=np.float64)
        coords = SkyCoord(l=l_vals * u.deg, b=b_vals * u.deg, frame='galactic').icrs
        idx = np.arange(l_vals.size, dtype=int)
        if idx.size > MAX_SELECTED_MEMBER_POINTS:
            choose = np.linspace(0, idx.size - 1, int(MAX_SELECTED_MEMBER_POINTS), dtype=int)
            idx = idx[choose]

        grouped[str(cluster_name)] = [
            {
                'l': float(l_vals[i]),
                'b': float(b_vals[i]),
                'ra': float(coords.ra.deg[i]),
                'dec': float(coords.dec.deg[i]),
                'label': str(cluster_name),
            }
            for i in idx
        ]

    return grouped or None


def _catalog_point_from_selection(selection, default_label=''):
    if not isinstance(selection, dict):
        return None

    l_deg = _coerce_float(selection.get('l_deg'), np.nan)
    b_deg = _coerce_float(selection.get('b_deg'), np.nan)
    ra_deg = _coerce_float(selection.get('ra_deg'), np.nan)
    dec_deg = _coerce_float(selection.get('dec_deg'), np.nan)
    if not (np.isfinite(l_deg) and np.isfinite(b_deg) and np.isfinite(ra_deg) and np.isfinite(dec_deg)):
        return None

    label = selection.get('cluster_name') or selection.get('trace_name') or default_label or 'Selection'
    return {
        'l': float(l_deg),
        'b': float(b_deg),
        'ra': float(ra_deg),
        'dec': float(dec_deg),
        'label': str(label),
    }


def _limit_catalog_points(points):
    if not points:
        return []
    if len(points) <= MAX_SELECTED_MEMBER_POINTS:
        return copy.deepcopy(points)

    idx = np.linspace(0, len(points) - 1, int(MAX_SELECTED_MEMBER_POINTS), dtype=int)
    return [copy.deepcopy(points[int(i)]) for i in idx]


def _threejs_catalog_from_frame_spec(frame_spec):
    if not isinstance(frame_spec, dict):
        return {}

    catalogs = {}
    for trace in frame_spec.get('traces', []):
        if not isinstance(trace, dict):
            continue
        trace_name = trace.get('name')
        trace_points = []
        grouped_points = {}

        for point in trace.get('points', []):
            if not isinstance(point, dict):
                continue
            selection = point.get('selection')
            catalog_point = _catalog_point_from_selection(selection, default_label=trace_name or '')
            if catalog_point is None:
                continue

            trace_points.append(catalog_point)

            cluster_name = selection.get('cluster_name') if isinstance(selection, dict) else None
            if cluster_name:
                grouped_points.setdefault(str(cluster_name), []).append(catalog_point)

        if trace_name and trace_points:
            catalogs[str(trace_name)] = _limit_catalog_points(trace_points)
        for cluster_name, points in grouped_points.items():
            catalogs.setdefault(cluster_name, _limit_catalog_points(points))

    return catalogs


def _merge_threejs_member_catalogs(*catalogs):
    merged = {}
    for catalog in catalogs:
        if not isinstance(catalog, dict):
            continue
        for key, points in catalog.items():
            if not key or not isinstance(points, list) or not points:
                continue
            merged[str(key)] = copy.deepcopy(points)
    return merged


def _normalize_threejs_volume_configs(volumes):
    if volumes in (None, False):
        return []

    if isinstance(volumes, (str, bytes)):
        return [{'path': str(volumes)}]

    if isinstance(volumes, dict):
        return [copy.deepcopy(volumes)]

    normalized = []
    for item in volumes:
        if item in (None, False):
            continue
        if isinstance(item, (str, bytes)):
            normalized.append({'path': str(item)})
        elif isinstance(item, dict):
            normalized.append(copy.deepcopy(item))
        else:
            raise TypeError(
                "volumes must be None, a path string, a dict, or a sequence of dict/path values."
            )
    return normalized


def _normalize_threejs_volume_stretch(stretch):
    value = str(stretch or 'linear').strip().lower()
    if value in {'linear', 'log10', 'asinh'}:
        return value
    return 'linear'


def _coerce_threejs_volume_time_myr(value):
    if value in (None, '', False):
        return None
    try:
        coerced = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(coerced):
        return None
    return coerced


def _threejs_volume_state_key(volume_cfg, fallback_key):
    candidate = (
        volume_cfg.get('state_key')
        or volume_cfg.get('legend_key')
        or volume_cfg.get('control_key')
        or volume_cfg.get('group_key')
        or fallback_key
    )
    candidate = str(candidate).strip()
    return candidate or str(fallback_key)


def _coerce_threejs_volume_bound_offset(bound_offset):
    if isinstance(bound_offset, dict):
        source = bound_offset
    else:
        source = {}
    return {
        'x': float(_coerce_float(source.get('x'), 0.0)),
        'y': float(_coerce_float(source.get('y'), 0.0)),
        'z': float(_coerce_float(source.get('z'), 0.0)),
    }


def _build_threejs_inline_volume_layer_spec(volume_cfg, center_offset=None, index=0):
    center_offset = center_offset or {'x': 0.0, 'y': 0.0, 'z': 0.0}
    data = volume_cfg.get('data')
    if data is None:
        raise ValueError("Each Three.js volume config must include a FITS 'path' or a 3D 'data' array.")

    data = np.asarray(data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(
            "Inline Three.js volume data must be a 3D array with axis order (z, y, x)."
        )

    data_shape_zyx = tuple(int(v) for v in data.shape)
    max_resolution = volume_cfg.get('max_resolution')
    if max_resolution is None:
        target_shape_zyx = data_shape_zyx
    else:
        max_resolution = int(np.clip(_coerce_float(max_resolution, max(data_shape_zyx)), 8, 512))
        target_shape_zyx = tuple(min(max_resolution, dim) for dim in data_shape_zyx)
    sampled = (
        _downsample_threejs_volume_with_zoom(data, target_shape_zyx)
        if target_shape_zyx != data_shape_zyx
        else np.ascontiguousarray(data)
    )

    finite_mask = np.isfinite(sampled)
    positive_mask = finite_mask & (sampled > 0)
    stats_values = sampled[positive_mask] if np.any(positive_mask) else sampled[finite_mask]
    if stats_values.size == 0:
        stats_values = np.array([0.0], dtype=float)

    data_range = volume_cfg.get('data_range')
    if data_range is not None:
        data_min = float(data_range[0])
        data_max = float(data_range[1])
    else:
        data_min = float(np.nanmin(stats_values))
        data_max = float(np.nanmax(stats_values))
    if not np.isfinite(data_min):
        data_min = 0.0
    if not np.isfinite(data_max) or not data_max > data_min:
        data_max = data_min + 1.0

    lower_quantile = float(np.clip(_coerce_float(volume_cfg.get('default_vmin_quantile'), 0.70), 0.0, 1.0))
    upper_quantile = float(np.clip(_coerce_float(volume_cfg.get('default_vmax_quantile'), 0.995), 0.0, 1.0))
    if upper_quantile < lower_quantile:
        upper_quantile = lower_quantile

    default_vmin = volume_cfg.get('vmin')
    default_vmax = volume_cfg.get('vmax')
    if default_vmin is None:
        default_vmin = float(np.nanquantile(stats_values, lower_quantile))
    else:
        default_vmin = float(default_vmin)
    if default_vmax is None:
        default_vmax = float(np.nanquantile(stats_values, upper_quantile))
    else:
        default_vmax = float(default_vmax)

    default_vmin = float(np.clip(default_vmin, data_min, data_max))
    default_vmax = float(np.clip(default_vmax, data_min, data_max))
    if not default_vmax > default_vmin:
        default_vmin = data_min
        default_vmax = data_max

    def _quantize_sampled_uint8(sampled_data):
        sampled_finite_mask = np.isfinite(sampled_data)
        if data_max > data_min:
            normalized = np.zeros_like(sampled_data, dtype=np.float32)
            normalized[sampled_finite_mask] = (sampled_data[sampled_finite_mask] - data_min) / (data_max - data_min)
            normalized = np.clip(normalized, 0.0, 1.0)
        else:
            normalized = np.zeros_like(sampled_data, dtype=np.float32)
            normalized[sampled_finite_mask] = 1.0
        quantized = np.zeros_like(sampled_data, dtype=np.uint8)
        quantized[sampled_finite_mask] = np.rint(normalized[sampled_finite_mask] * 255.0).astype(np.uint8)
        return np.ascontiguousarray(quantized)

    def _encode_sampled_uint8(sampled_data):
        quantized = _quantize_sampled_uint8(sampled_data)
        return base64.b64encode(quantized.tobytes(order='C')).decode('ascii')

    bounds_cfg = volume_cfg.get('bounds') or {}
    x_bounds = _coerce_range(bounds_cfg.get('x'), [-0.5 * data_shape_zyx[2], 0.5 * data_shape_zyx[2]])
    y_bounds = _coerce_range(bounds_cfg.get('y'), [-0.5 * data_shape_zyx[1], 0.5 * data_shape_zyx[1]])
    z_bounds = _coerce_range(bounds_cfg.get('z'), [-0.5 * data_shape_zyx[0], 0.5 * data_shape_zyx[0]])
    apply_center_offset = bool(volume_cfg.get('apply_center_offset', False))
    bound_offset = _coerce_threejs_volume_bound_offset(volume_cfg.get('bound_offset'))

    sample_setting = volume_cfg.get('samples', volume_cfg.get('steps', 192))
    default_steps = int(np.clip(_coerce_float(sample_setting, 192), 24, 768))
    default_opacity = float(np.clip(_coerce_float(volume_cfg.get('opacity'), 0.15), 0.0, 1.0))
    default_alpha_coef = float(np.clip(_coerce_float(volume_cfg.get('alpha_coef'), 50.0), 1.0, 200.0))
    default_gradient_step = float(np.clip(_coerce_float(volume_cfg.get('gradient_step'), 0.005), 1e-4, 0.05))
    default_stretch = _normalize_threejs_volume_stretch(volume_cfg.get('stretch', 'linear'))

    name = str(volume_cfg.get('name') or f'Volume {index + 1}')
    key = str(volume_cfg.get('key') or f'volume-{index}')
    state_key = _threejs_volume_state_key(volume_cfg, key)
    state_name = str(volume_cfg.get('state_name') or volume_cfg.get('legend_name') or name)
    time_myr = _coerce_threejs_volume_time_myr(volume_cfg.get('time_myr'))
    opacity_function = _normalize_threejs_volume_opacity_function(volume_cfg.get('opacity_function'))
    colormap_options = _build_threejs_volume_colormap_options(
        volume_cfg.get('colormap', 'inferno'),
        opacity_function=opacity_function,
    )

    return {
        'key': key,
        'state_key': state_key,
        'state_name': state_name,
        'time_myr': time_myr,
        'name': name,
        'path': '<inline>',
        'hdu': 'inline',
        'data_b64': _encode_sampled_uint8(sampled),
        'data_encoding': 'uint8',
        'shape': {
            'x': int(sampled.shape[2]),
            'y': int(sampled.shape[1]),
            'z': int(sampled.shape[0]),
        },
        'source_shape': {
            'x': int(data_shape_zyx[2]),
            'y': int(data_shape_zyx[1]),
            'z': int(data_shape_zyx[0]),
        },
        'downsample_step': {
            'x': float(data_shape_zyx[2]) / float(sampled.shape[2]),
            'y': float(data_shape_zyx[1]) / float(sampled.shape[1]),
            'z': float(data_shape_zyx[0]) / float(sampled.shape[0]),
        },
        'downsample_method': 'scipy_zoom' if sampled.shape != data.shape else 'inline',
        'bounds': {
            'x': [
                float((x_bounds[0] + bound_offset['x']) - center_offset['x'])
                if apply_center_offset else float(x_bounds[0] + bound_offset['x']),
                float((x_bounds[1] + bound_offset['x']) - center_offset['x'])
                if apply_center_offset else float(x_bounds[1] + bound_offset['x']),
            ],
            'y': [
                float((y_bounds[0] + bound_offset['y']) - center_offset['y'])
                if apply_center_offset else float(y_bounds[0] + bound_offset['y']),
                float((y_bounds[1] + bound_offset['y']) - center_offset['y'])
                if apply_center_offset else float(y_bounds[1] + bound_offset['y']),
            ],
            'z': [
                float((z_bounds[0] + bound_offset['z']) - center_offset['z'])
                if apply_center_offset else float(z_bounds[0] + bound_offset['z']),
                float((z_bounds[1] + bound_offset['z']) - center_offset['z'])
                if apply_center_offset else float(z_bounds[1] + bound_offset['z']),
            ],
        },
        'data_range': [float(data_min), float(data_max)],
        'value_unit': str(volume_cfg.get('unit_label') or volume_cfg.get('value_unit') or '').strip(),
        'legend_color': colormap_options[0].get('legend_color'),
        'visible': bool(volume_cfg.get('visible', True)),
        'only_at_t0': bool(volume_cfg.get('only_at_t0', time_myr is None)),
        'interpolation': bool(volume_cfg.get('interpolation', True)),
        'sky_overlay_data_b64': None,
        'sky_overlay_data_encoding': None,
        'sky_overlay_shape': None,
        'sky_overlay_atlas_tiles': None,
        'sky_overlay_downsample_step': None,
        'default_controls': {
            'vmin': float(default_vmin),
            'vmax': float(default_vmax),
            'opacity': float(default_opacity),
            'steps': int(default_steps),
            'samples': int(default_steps),
            'alpha_coef': float(default_alpha_coef),
            'gradient_step': float(default_gradient_step),
            'stretch': str(default_stretch),
            'colormap': colormap_options[0]['name'],
        },
        'colormap_options': colormap_options,
    }


def _build_threejs_volume_layer_spec(volume_cfg, center_offset=None, index=0, include_sky_overlay=False):
    from astropy.io import fits

    center_offset = center_offset or {'x': 0.0, 'y': 0.0, 'z': 0.0}
    path = volume_cfg.get('path') or volume_cfg.get('fits_file')
    if not path:
        return _build_threejs_inline_volume_layer_spec(volume_cfg, center_offset=center_offset, index=index)

    hdu_selector = volume_cfg.get('hdu')
    max_resolution = int(np.clip(_coerce_float(volume_cfg.get('max_resolution'), 96), 24, 512))
    sky_overlay_max_resolution = None
    if include_sky_overlay:
        sky_overlay_max_resolution = int(
            np.clip(
                _coerce_float(volume_cfg.get('sky_overlay_max_resolution'), max(max_resolution, 256)),
                64,
                512,
            )
        )
    apply_center_offset = bool(volume_cfg.get('apply_center_offset', True))
    bound_offset = _coerce_threejs_volume_bound_offset(volume_cfg.get('bound_offset'))
    sample_setting = volume_cfg.get('samples', volume_cfg.get('steps', 192))
    default_steps = int(np.clip(_coerce_float(sample_setting, 192), 24, 768))
    default_opacity = float(np.clip(_coerce_float(volume_cfg.get('opacity'), 0.15), 0.0, 1.0))
    default_alpha_coef = float(np.clip(_coerce_float(volume_cfg.get('alpha_coef'), 50.0), 1.0, 200.0))
    default_gradient_step = float(np.clip(_coerce_float(volume_cfg.get('gradient_step'), 0.005), 1e-4, 0.05))
    default_stretch = _normalize_threejs_volume_stretch(volume_cfg.get('stretch', 'linear'))

    with fits.open(path, memmap=True) as hdul:
        hdu_info = _resolve_threejs_volume_hdu(hdul, hdu_selector)
        hdu = hdu_info['hdu']
        data = hdu_info['cube_data_zyx']
        data_shape_zyx = tuple(int(v) for v in data.shape)
        target_shape_zyx = tuple(min(max_resolution, dim) for dim in data_shape_zyx)
        sampled = _downsample_threejs_volume_with_zoom(data, target_shape_zyx)
        sky_overlay_shape_zyx = None
        sky_overlay_sampled = None
        if sky_overlay_max_resolution is not None:
            sky_overlay_shape_zyx = tuple(min(sky_overlay_max_resolution, dim) for dim in data_shape_zyx)
            if sky_overlay_shape_zyx == target_shape_zyx:
                sky_overlay_sampled = sampled
            else:
                sky_overlay_sampled = _downsample_threejs_volume_with_zoom(data, sky_overlay_shape_zyx)
        header = hdu.header.copy()
        axis_numbers_zyx = hdu_info['axis_numbers_zyx']
        resolved_hdu_label = hdu_info['resolved_hdu']

    stats_source = sky_overlay_sampled if sky_overlay_sampled is not None else sampled
    finite_mask = np.isfinite(stats_source)
    positive_mask = finite_mask & (stats_source > 0)
    stats_values = stats_source[positive_mask] if np.any(positive_mask) else stats_source[finite_mask]
    if stats_values.size == 0:
        raise ValueError(f"Volume FITS cube contains no finite samples: {path}")

    data_min = float(np.nanmin(stats_values))
    data_max = float(np.nanmax(stats_values))
    if not np.isfinite(data_min) or not np.isfinite(data_max):
        raise ValueError(f"Volume FITS cube statistics are not finite: {path}")

    lower_quantile = float(np.clip(_coerce_float(volume_cfg.get('default_vmin_quantile'), 0.70), 0.0, 1.0))
    upper_quantile = float(np.clip(_coerce_float(volume_cfg.get('default_vmax_quantile'), 0.995), 0.0, 1.0))
    if upper_quantile < lower_quantile:
        upper_quantile = lower_quantile

    default_vmin = volume_cfg.get('vmin')
    default_vmax = volume_cfg.get('vmax')
    if default_vmin is None:
        default_vmin = float(np.nanquantile(stats_values, lower_quantile))
    else:
        default_vmin = float(default_vmin)
    if default_vmax is None:
        default_vmax = float(np.nanquantile(stats_values, upper_quantile))
    else:
        default_vmax = float(default_vmax)

    default_vmin = float(np.clip(default_vmin, data_min, data_max))
    default_vmax = float(np.clip(default_vmax, data_min, data_max))
    if not default_vmax > default_vmin:
        if data_max > data_min:
            default_vmin = data_min
            default_vmax = data_max
        else:
            default_vmax = default_vmin + 1.0

    def _quantize_sampled_uint8(sampled_data):
        sampled_finite_mask = np.isfinite(sampled_data)
        if data_max > data_min:
            normalized = np.zeros_like(sampled_data, dtype=np.float32)
            normalized[sampled_finite_mask] = (sampled_data[sampled_finite_mask] - data_min) / (data_max - data_min)
            normalized = np.clip(normalized, 0.0, 1.0)
        else:
            normalized = np.zeros_like(sampled_data, dtype=np.float32)
            normalized[sampled_finite_mask] = 1.0
        quantized = np.zeros_like(sampled_data, dtype=np.uint8)
        quantized[sampled_finite_mask] = np.rint(normalized[sampled_finite_mask] * 255.0).astype(np.uint8)
        return np.ascontiguousarray(quantized)

    def _encode_sampled_uint8(sampled_data):
        quantized = _quantize_sampled_uint8(sampled_data)
        return base64.b64encode(quantized.tobytes(order='C')).decode('ascii')

    def _encode_sampled_uint8_png_atlas(sampled_data):
        try:
            from PIL import Image
        except ImportError:
            return _encode_sampled_uint8(sampled_data), 'uint8', None

        quantized = _quantize_sampled_uint8(sampled_data)
        nz, ny, nx = (int(quantized.shape[0]), int(quantized.shape[1]), int(quantized.shape[2]))
        tile_cols = int(np.ceil(np.sqrt(max(nz, 1))))
        tile_rows = int(np.ceil(nz / tile_cols))
        atlas = np.zeros((tile_rows * ny, tile_cols * nx), dtype=np.uint8)
        for z_index in range(nz):
            row = z_index // tile_cols
            col = z_index % tile_cols
            atlas[row * ny:(row + 1) * ny, col * nx:(col + 1) * nx] = quantized[z_index]

        buffer = io.BytesIO()
        Image.fromarray(atlas, mode='L').save(
            buffer,
            format='PNG',
            optimize=True,
            compress_level=9,
        )
        return (
            base64.b64encode(buffer.getvalue()).decode('ascii'),
            'png_atlas_uint8',
            {'x': tile_cols, 'y': tile_rows},
        )

    data_b64 = _encode_sampled_uint8(sampled)
    sky_overlay_b64 = None
    sky_overlay_encoding = None
    sky_overlay_tiles = None
    if sky_overlay_sampled is not None:
        sky_overlay_b64, sky_overlay_encoding, sky_overlay_tiles = _encode_sampled_uint8_png_atlas(sky_overlay_sampled)

    x_bounds = _threejs_volume_axis_bounds(header, axis_number=int(axis_numbers_zyx[2]), axis_size=int(data_shape_zyx[2]))
    y_bounds = _threejs_volume_axis_bounds(header, axis_number=int(axis_numbers_zyx[1]), axis_size=int(data_shape_zyx[1]))
    z_bounds = _threejs_volume_axis_bounds(header, axis_number=int(axis_numbers_zyx[0]), axis_size=int(data_shape_zyx[0]))

    opacity_function = _normalize_threejs_volume_opacity_function(volume_cfg.get('opacity_function'))
    colormap_options = _build_threejs_volume_colormap_options(
        volume_cfg.get('colormap', 'inferno'),
        opacity_function=opacity_function,
    )
    name = str(volume_cfg.get('name') or str(path).rsplit('/', 1)[-1].rsplit('.', 1)[0] or f'Volume {index + 1}')
    value_unit = str(volume_cfg.get('unit_label') or header.get('BUNIT') or '').strip()

    key = str(volume_cfg.get('key') or f'volume-{index}')
    state_key = _threejs_volume_state_key(volume_cfg, key)
    time_myr = _coerce_threejs_volume_time_myr(volume_cfg.get('time_myr'))
    state_name = str(volume_cfg.get('state_name') or volume_cfg.get('legend_name') or name)

    return {
        'key': key,
        'state_key': state_key,
        'state_name': state_name,
        'time_myr': time_myr,
        'name': name,
        'path': str(path),
        'hdu': str(resolved_hdu_label),
        'data_b64': data_b64,
        'data_encoding': 'uint8',
        'shape': {
            'x': int(sampled.shape[2]),
            'y': int(sampled.shape[1]),
            'z': int(sampled.shape[0]),
        },
        'source_shape': {
            'x': int(data_shape_zyx[2]),
            'y': int(data_shape_zyx[1]),
            'z': int(data_shape_zyx[0]),
        },
        'downsample_step': {
            'x': float(data_shape_zyx[2]) / float(sampled.shape[2]),
            'y': float(data_shape_zyx[1]) / float(sampled.shape[1]),
            'z': float(data_shape_zyx[0]) / float(sampled.shape[0]),
        },
        'downsample_method': 'scipy_zoom',
        'bounds': {
            'x': [
                float((x_bounds[0] + bound_offset['x']) - center_offset['x'])
                if apply_center_offset else float(x_bounds[0] + bound_offset['x']),
                float((x_bounds[1] + bound_offset['x']) - center_offset['x'])
                if apply_center_offset else float(x_bounds[1] + bound_offset['x']),
            ],
            'y': [
                float((y_bounds[0] + bound_offset['y']) - center_offset['y'])
                if apply_center_offset else float(y_bounds[0] + bound_offset['y']),
                float((y_bounds[1] + bound_offset['y']) - center_offset['y'])
                if apply_center_offset else float(y_bounds[1] + bound_offset['y']),
            ],
            'z': [
                float((z_bounds[0] + bound_offset['z']) - center_offset['z'])
                if apply_center_offset else float(z_bounds[0] + bound_offset['z']),
                float((z_bounds[1] + bound_offset['z']) - center_offset['z'])
                if apply_center_offset else float(z_bounds[1] + bound_offset['z']),
            ],
        },
        'data_range': [float(data_min), float(data_max)],
        'value_unit': value_unit,
        'legend_color': colormap_options[0].get('legend_color'),
        'visible': bool(volume_cfg.get('visible', True)),
        'only_at_t0': bool(volume_cfg.get('only_at_t0', time_myr is None)),
        'interpolation': bool(volume_cfg.get('interpolation', True)),
        'sky_overlay_data_b64': sky_overlay_b64,
        'sky_overlay_data_encoding': sky_overlay_encoding,
        'sky_overlay_shape': (
            {
                'x': int(sky_overlay_sampled.shape[2]),
                'y': int(sky_overlay_sampled.shape[1]),
                'z': int(sky_overlay_sampled.shape[0]),
            }
            if sky_overlay_sampled is not None
            else None
        ),
        'sky_overlay_atlas_tiles': sky_overlay_tiles,
        'sky_overlay_downsample_step': (
            {
                'x': float(data_shape_zyx[2]) / float(sky_overlay_sampled.shape[2]),
                'y': float(data_shape_zyx[1]) / float(sky_overlay_sampled.shape[1]),
                'z': float(data_shape_zyx[0]) / float(sky_overlay_sampled.shape[0]),
            }
            if sky_overlay_sampled is not None
            else None
        ),
        'default_controls': {
            'vmin': float(default_vmin),
            'vmax': float(default_vmax),
            'opacity': float(default_opacity),
            'steps': int(default_steps),
            'samples': int(default_steps),
            'alpha_coef': float(default_alpha_coef),
            'gradient_step': float(default_gradient_step),
            'stretch': str(default_stretch),
            'colormap': colormap_options[0]['name'],
        },
        'colormap_options': colormap_options,
    }


def _resolve_threejs_volume_hdu(hdul, hdu_selector):
    if hdu_selector not in (None, '', 'auto', 'AUTO'):
        resolved_hdu = _resolve_threejs_volume_hdu_explicit(hdul, hdu_selector)
        cube_data_zyx, axis_numbers_zyx = _coerce_threejs_volume_cube(resolved_hdu.data)
        return {
            'hdu': resolved_hdu,
            'resolved_hdu': _threejs_volume_hdu_label(hdul, resolved_hdu),
            'cube_data_zyx': cube_data_zyx,
            'axis_numbers_zyx': axis_numbers_zyx,
            'auto_selected': False,
        }

    candidates = _threejs_volume_hdu_candidates(hdul)
    if not candidates:
        available = ', '.join(
            f"[{idx}] {type(hdu).__name__} {getattr(hdu, 'name', '')!s}".strip()
            for idx, hdu in enumerate(hdul)
        )
        raise ValueError(
            "Could not auto-detect a 3D FITS image cube. "
            f"Available HDUs: {available or 'none'}."
        )

    selected = max(candidates, key=_threejs_volume_hdu_candidate_sort_key)
    return {
        'hdu': selected['hdu'],
        'resolved_hdu': selected['label'],
        'cube_data_zyx': selected['cube_data_zyx'],
        'axis_numbers_zyx': selected['axis_numbers_zyx'],
        'auto_selected': True,
    }


def _resolve_threejs_volume_hdu_explicit(hdul, hdu_selector):
    if isinstance(hdu_selector, str):
        target = hdu_selector.strip().lower()
        for hdu in hdul:
            if str(getattr(hdu, 'name', '')).strip().lower() == target:
                return hdu
        if target.isdigit():
            return hdul[int(target)]
        available = ', '.join(_threejs_volume_hdu_label(hdul, hdu) for hdu in hdul)
        raise KeyError(
            f"Could not find FITS HDU named '{hdu_selector}'. "
            f"Available HDUs: {available or 'none'}. "
            "Omit 'hdu' to auto-select a 3D cube."
        )
    return hdul[int(hdu_selector)]


def _threejs_volume_hdu_candidates(hdul):
    candidates = []
    for idx, hdu in enumerate(hdul):
        data = getattr(hdu, 'data', None)
        if data is None:
            continue
        try:
            cube_data_zyx, axis_numbers_zyx = _coerce_threejs_volume_cube(data)
        except ValueError:
            continue
        shape_zyx = tuple(int(v) for v in cube_data_zyx.shape)
        name = str(getattr(hdu, 'name', '') or '').strip()
        label = _threejs_volume_hdu_label(hdul, hdu)
        name_lower = name.lower()
        positive_hints = ('mean', 'dust', 'density', 'map', 'cube', 'data', 'signal', 'primary')
        negative_hints = ('std', 'sigma', 'var', 'variance', 'error', 'err', 'uncert', 'mask', 'weight')
        candidates.append({
            'index': int(idx),
            'hdu': hdu,
            'label': label,
            'cube_data_zyx': cube_data_zyx,
            'axis_numbers_zyx': axis_numbers_zyx,
            'shape_zyx': shape_zyx,
            'voxel_count': int(np.prod(shape_zyx, dtype=np.int64)),
            'positive_name_hint': any(hint in name_lower for hint in positive_hints),
            'negative_name_hint': any(hint in name_lower for hint in negative_hints),
            'ndim': int(np.ndim(data)),
            'image_hdu': hasattr(hdu, 'is_image') and bool(hdu.is_image),
        })
    return candidates


def _threejs_volume_hdu_candidate_sort_key(candidate):
    return (
        int(bool(candidate.get('positive_name_hint'))),
        -int(bool(candidate.get('negative_name_hint'))),
        int(bool(candidate.get('image_hdu'))),
        int(candidate.get('voxel_count', 0)),
        -abs(int(candidate.get('ndim', 3)) - 3),
        -int(candidate.get('index', 0)),
    )


def _threejs_volume_hdu_label(hdul, hdu):
    try:
        index = int(hdul.index_of(hdu))
    except Exception:
        index = next((idx for idx, item in enumerate(hdul) if item is hdu), -1)
    name = str(getattr(hdu, 'name', '') or '').strip()
    if name:
        return name
    return str(index if index >= 0 else 0)


def _coerce_threejs_volume_cube(data):
    arr = np.asarray(data)
    if arr.size == 0:
        raise ValueError("FITS volume HDU is empty.")
    if not np.issubdtype(arr.dtype, np.number):
        raise ValueError("FITS volume HDU data is not numeric.")

    fits_axis_numbers = [arr.ndim - idx for idx in range(arr.ndim)]
    kept_axis_numbers = [
        fits_axis_numbers[idx]
        for idx, size in enumerate(arr.shape)
        if int(size) != 1
    ]
    squeezed = np.squeeze(arr)
    if squeezed.ndim != 3:
        raise ValueError("FITS volume HDU is not 3D after removing singleton axes.")
    if len(kept_axis_numbers) != 3:
        raise ValueError("Could not map squeezed FITS cube axes.")
    return np.asarray(squeezed), tuple(int(v) for v in kept_axis_numbers)


def _downsample_threejs_volume_with_zoom(data, target_shape_zyx):
    from scipy.ndimage import zoom

    source_shape = tuple(int(v) for v in np.shape(data))
    if source_shape != tuple(int(v) for v in target_shape_zyx):
        zoom_factors = tuple(
            float(target_dim) / float(source_dim)
            for source_dim, target_dim in zip(source_shape, target_shape_zyx)
        )
        sampled = zoom(
            data,
            zoom=zoom_factors,
            order=1,
            mode='nearest',
            prefilter=False,
            grid_mode=True,
            output=np.float32,
        )
        return np.asarray(sampled, dtype=np.float32)
    return np.asarray(data, dtype=np.float32)


def _threejs_volume_axis_bounds(header, axis_number, axis_size):
    cdelt_key = f'CDELT{axis_number}'
    crval_key = f'CRVAL{axis_number}'
    crpix_key = f'CRPIX{axis_number}'
    if (cdelt_key not in header) or (crval_key not in header) or (crpix_key not in header):
        half_size = 0.5 * float(axis_size)
        return float(-half_size), float(half_size)

    delta = float(header.get(cdelt_key, 1.0))
    crval = float(header.get(crval_key, 0.0))
    crpix = float(header.get(crpix_key, 1.0))
    if (not np.isfinite(delta)) or (abs(delta) <= 0.0) or (not np.isfinite(crval)) or (not np.isfinite(crpix)):
        half_size = 0.5 * float(axis_size)
        return float(-half_size), float(half_size)
    first_center = crval + ((1.0 - crpix) * delta)
    last_center = crval + ((float(axis_size) - crpix) * delta)
    edge_pad = 0.5 * delta
    lower = min(first_center, last_center) - abs(edge_pad)
    upper = max(first_center, last_center) + abs(edge_pad)
    return float(lower), float(upper)


def _normalize_threejs_volume_opacity_function(opacity_function):
    if opacity_function is None or opacity_function is False:
        return np.asarray([[0.0, 0.0], [1.0, 1.0]], dtype=float)
    if isinstance(opacity_function, str) and not opacity_function.strip():
        return np.asarray([[0.0, 0.0], [1.0, 1.0]], dtype=float)

    values = np.asarray(opacity_function, dtype=float)
    if values.ndim == 1:
        if values.size < 4 or values.size % 2 != 0:
            raise ValueError(
                "Three.js volume opacity_function must contain an even number of position/alpha values."
            )
        values = values.reshape(-1, 2)
    elif values.ndim != 2 or values.shape[1] != 2:
        raise ValueError(
            "Three.js volume opacity_function must be a flat [x0, a0, x1, a1, ...] sequence or an Nx2 array."
        )

    values = values[np.all(np.isfinite(values), axis=1)]
    if values.shape[0] < 2:
        return np.asarray([[0.0, 0.0], [1.0, 1.0]], dtype=float)

    positions = np.clip(values[:, 0], 0.0, 1.0)
    alphas = np.clip(values[:, 1], 0.0, 1.0)
    order = np.argsort(positions, kind='mergesort')
    merged = np.column_stack((positions[order], alphas[order]))

    dedup_positions = []
    dedup_alphas = []
    for position, alpha in merged:
        if dedup_positions and abs(position - dedup_positions[-1]) <= 1e-12:
            dedup_alphas[-1] = float(alpha)
            continue
        dedup_positions.append(float(position))
        dedup_alphas.append(float(alpha))

    merged = np.column_stack((dedup_positions, dedup_alphas))
    if merged[0, 0] > 0.0:
        merged = np.vstack(([0.0, merged[0, 1]], merged))
    if merged[-1, 0] < 1.0:
        merged = np.vstack((merged, [1.0, merged[-1, 1]]))
    return np.asarray(merged, dtype=float)


def _build_threejs_volume_colormap_options(selected_colormap, opacity_function=None):
    options = []
    seen = set()
    requested = [selected_colormap, *DEFAULT_THREEJS_VOLUME_COLORMAPS]

    for candidate in requested:
        if candidate in (None, ''):
            continue
        sampled = _sample_threejs_volume_colormap(candidate, opacity_function=opacity_function)
        if sampled is None:
            continue
        name, rgba = sampled
        key = str(name).strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        legend_idx = int(np.clip(round(0.72 * (rgba.shape[0] - 1)), 0, rgba.shape[0] - 1))
        legend_rgb = rgba[legend_idx, :3].astype(int).tolist()
        options.append({
            'name': str(name),
            'label': str(name),
            'lut_b64': base64.b64encode(np.ascontiguousarray(rgba).tobytes(order='C')).decode('ascii'),
            'legend_color': f'rgb({legend_rgb[0]}, {legend_rgb[1]}, {legend_rgb[2]})',
        })

    if not options:
        raise ValueError("Could not resolve any colormap options for the Three.js volume renderer.")
    return options


def _sample_threejs_volume_colormap(colormap_value, n_samples=1024, opacity_function=None):
    values = np.linspace(0.0, 1.0, int(n_samples), dtype=float)

    cmap_obj = _resolve_threejs_volume_colormap_object(colormap_value)
    if cmap_obj is None:
        return None

    label, cmap_callable = cmap_obj
    sampled = np.asarray(cmap_callable(values), dtype=float)
    if sampled.ndim != 2 or sampled.shape[0] != values.size:
        return None
    if sampled.shape[1] == 3:
        alpha = np.ones((sampled.shape[0], 1), dtype=float)
        sampled = np.concatenate((sampled, alpha), axis=1)
    elif sampled.shape[1] != 4:
        return None

    opacity_curve = _normalize_threejs_volume_opacity_function(opacity_function)
    sampled[:, 3] *= np.interp(values, opacity_curve[:, 0], opacity_curve[:, 1])
    rgba = np.clip(np.rint(sampled * 255.0), 0.0, 255.0).astype(np.uint8)
    return str(label), np.ascontiguousarray(rgba)


def _resolve_threejs_volume_colormap_object(colormap_value):
    if colormap_value is None:
        colormap_value = 'inferno'

    if hasattr(colormap_value, '__call__') and not isinstance(colormap_value, str):
        label = getattr(colormap_value, 'name', None) or 'custom'
        return str(label), colormap_value

    if not isinstance(colormap_value, str):
        return None

    requested = str(colormap_value).strip()
    if not requested:
        return None

    candidates = []
    for name in (requested, requested.lower(), requested.capitalize()):
        if name not in candidates:
            candidates.append(name)

    try:
        from matplotlib import colormaps as mpl_colormaps

        for name in candidates:
            try:
                cmap = mpl_colormaps[name]
                return getattr(cmap, 'name', name), cmap
            except Exception:
                continue
    except Exception:
        pass

    try:
        import colormaps as installed_colormaps

        for name in candidates:
            cmap = getattr(installed_colormaps, name, None)
            if callable(cmap):
                return getattr(cmap, 'name', name), cmap
    except Exception:
        pass

    return None


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
