"""Build time-dependent Oviz scenes from traces, volumes, and sky layers."""

import base64
import contextlib
import copy
import functools
import gc
import importlib.resources
import io
import math
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import webcolors
import yaml
from astropy.coordinates import SkyCoord

from . import orbit_maker
from .spiral_models import CASTRO_GINARD2021_ARM_TABLE, spiral_arm_coordinates
from .threejs_figure import ThreeJSFigure
from .threejs_profiles import build_threejs_profile, merge_threejs_profile
from .threejs_scene import build_threejs_scene_spec


# Constants -------------------------------------------------------------------

DEFAULT_VIEWER = "oviz"
SUN_TRACE_NAME = 'Sun'

GALACTIC_GUIDE_TRACE_NAMES = {
    'Galactic Quadrants',
    'Galactic l Labels',
}

GALACTIC_RADIUS_CIRCLE_DEFS = (
    (4.0, 'R = 4 kpc'),
    (8.122, 'R = 8.12 kpc'),
    (12.0, 'R = 12 kpc'),
)

GALACTIC_RADIUS_TRACE_NAMES = (
    {name for _, name in GALACTIC_RADIUS_CIRCLE_DEFS}
    | {f'{name} Label' for _, name in GALACTIC_RADIUS_CIRCLE_DEFS}
    | {'GC Ring'}
)

GALACTIC_SIMPLE_ALLOWED_TRACE_NAMES = {
    'Sun',
    'Clusters (< 60 Myr)',
    'R = 8.12 kpc',
}

# The Castro-Ginard et al. (2021) arms drawn by ``make_plot(include_spiral_arms=True)``.
SPIRAL_ARMS = {
    arm: {
        "theta_ref_deg": theta_ref,
        "theta_range_deg": theta_range,
        "Rref_kpc": r_ref,
        "psi_deg": pitch,
        "Omega_p": pattern_speed,
    }
    for arm, theta_ref, theta_range, r_ref, pitch, pattern_speed in CASTRO_GINARD2021_ARM_TABLE
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
CUSTOMDATA_IDX_N_STARS = 10
CUSTOMDATA_IDX_CLUSTER_ALIASES = 11
MAX_SELECTED_MEMBER_POINTS = 1200

DEFAULT_THREEJS_VOLUME_COLORMAPS = (
    'inferno',
    'magma',
    'plasma',
    'viridis',
    'cividis',
    'turbo',
    'gist_heat',
    'Greys',
)

DEFAULT_THREEJS_TRACE_COLORMAP = 'turbo'
DEFAULT_THREEJS_VOLUME_SAMPLE_STEPS = 100
DEFAULT_THREEJS_VOLUME_MAX_RESOLUTION_CAP = 512
DEFAULT_THREEJS_AR_VOLUME_MAX_RESOLUTION = 64


# Trace, frame and layout records ---------------------------------------------


class _OvizAttrDict(dict):
    """A dict with attribute access: the trace, frame and layout records ``make_plot`` builds."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __setattr__(self, name, value):
        self[name] = value

    def to_scene_json(self):
        return copy.deepcopy(dict(self))


def _attrdict_factory(kind):
    """Constructor of ``kind`` records: keyword fields, or a dict that is deep-copied."""
    def _factory(*args, **kwargs):
        data = {}
        if args:
            if len(args) != 1 or not isinstance(args[0], dict):
                raise TypeError(f"{kind} accepts at most one dict positional argument")
            data.update(copy.deepcopy(args[0]))
        data.update(kwargs)
        return _OvizAttrDict(data)
    return _factory


_scatter3d = _attrdict_factory("Scatter3d")
_scatter = _attrdict_factory("Scatter")
_layout = _attrdict_factory("Layout")
_frame = _attrdict_factory("Frame")


# The default Sun -------------------------------------------------------------


def _default_sun_trace():
    """The Sun as a one-object trace at the origin of heliocentric coordinates.

    Its orbit is integrated like any cluster's, so in the LSR-centred frame it
    moves with the Sun's peculiar motion. It matches the Sun traces the
    figure scripts build (yellow, 4–8 px, always born).
    """
    from .traces import Trace

    sun = pd.DataFrame({
        'name': [SUN_TRACE_NAME], 'age_myr': [4600.0], 'n_stars': [1],
        'x': [0.0], 'y': [0.0], 'z': [0.0], 'U': [0.0], 'V': [0.0], 'W': [0.0],
    })
    return Trace(
        sun, data_name=SUN_TRACE_NAME, min_size=4.0, max_size=8.0, color='yellow',
        opacity=1.0, marker_style='circle', show_tracks=False, size_by_n_stars=False,
    )


class _CollectionWithSun:
    """The user's trace collection plus Oviz's default Sun (listed first), for
    building one figure. The user's collection itself is never changed: its
    traces, their order and indices stay as they were."""

    def __init__(self, base, sun):
        self.base = base
        self.sun = sun

    def get_all_clusters(self):
        return [self.sun, *self.base.get_all_clusters()]

    def get_cluster(self, identifier):
        if isinstance(identifier, str) and identifier == self.sun.data_name:
            return self.sun
        return self.base.get_cluster(identifier)

    def __getattr__(self, name):
        return getattr(self.base, name)


# Scene builder ---------------------------------------------------------------


class Animate3D:
    """Build an interactive, time-dependent figure from a trace collection.

    :meth:`make_plot` integrates the traces' orbits with galpy, builds one
    frame per time step (with Galactic guides, static traces and volumes) and
    returns the figure. :class:`oviz.Scene3D` is the same class under an
    astronomy-facing name.

    Parameters
    ----------
    data_collection : TraceCollection
        The traces, in legend order.
    xyz_widths : tuple of float
        Half-widths (pc) of the initial x, y and z ranges.
    xyz_ranges : tuple of (low, high) pairs, optional
        Explicit x, y and z ranges (pc), used instead of ``xyz_widths``.
    figure_title : str, optional
        Title shown above the figure.
    figure_theme : str
        A theme in ``oviz/themes``: ``"dark"``, ``"light"``, ``"gray"`` or
        ``"solarized_light"``.
    trace_grouping_dict : dict, optional
        Legend groups, each a list of trace names. An ``"All"`` group with
        every trace is always added.
    potential : galpy potential, optional
        Defaults to ``MWPotential2014``.
    vo, ro, zo : float
        Circular velocity (km/s), solar radius (kpc) and solar height (kpc).

    ``light_template`` sets the light theme's layout template;
    ``figure_layout`` and ``figure_layout_dict`` are kept for compatibility
    and have no effect (the theme defines the layout).
    """

    def __init__(
        self,
        data_collection,
        xyz_widths=(1000, 1000, 300),
        xyz_ranges=None,
        figure_title=None,
        figure_theme=None,
        light_template=None,
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
        self.light_template = light_template
        self.trace_grouping_dict = trace_grouping_dict or {}
        self.figure_title = figure_title

        # Read in the layout from the theme
        self.figure_layout_dict = read_theme(self)

        # Configure the main figure title if provided
        if self.figure_title:
            font_color = 'white' if self.figure_theme == 'dark' else 'black'
            self.figure_layout_dict['title'] = dict(
                text=self.figure_title,
                x=0.5,
                font=dict(family='Helvetica', size=20, color=font_color)
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
        spiral_arm_models=None,
        show_age_kde_inset=False,
        age_kde_bandwidth_myr=3.0,
        camera_zoom_factor=1.0,
        galactic_reference_opacity=0.5,
        renderer='threejs',
        show_milky_way_model=False,
        enable_sky_panel=False,
        sky_radius_deg=1.0,
        sky_frame='galactic',
        sky_survey='P/DSS2/color',
        cluster_members_file=None,
        show_cluster_members_in_sky=False,
        volumes=None,
        threejs_initial_state=None,
        preset=None,
        actions=None,
        compress_scene_spec="auto",
        scene_spec_compression_threshold_bytes=None,
        viewer=None,
        viewer_mode=None,
        camera_anchor=None,
        show_sun=True,
    ):
        """Integrate the data, build timeline frames, and return a figure.

        Parameters
        ----------
        time : array-like
            Timeline values in Myr. The array must include zero.
        renderer : str, default ``"threejs"``
            Output renderer. Three.js is the maintained interactive path.
        galactic_mode : bool
            Use the Galactic-scale layout and reference geometry.
        enable_sky_panel : bool
            Enable the registered Aladin Lite Sky view.
        cluster_members_file : path-like, optional
            CSV containing cluster labels and Galactic or ICRS member positions.
        show_cluster_members_in_sky : bool
            Replace eligible bulk cluster markers with member stars in Sky mode.
        volumes : sequence of mappings, optional
            Scalar-volume definitions, including data source, bounds, and style.
        threejs_initial_state : mapping, optional
            Initial camera, control, layer, panel, and presentation settings.
        actions : sequence of mappings, optional
            Declarative camera, legend, State, or timeline Actions.
        compress_scene_spec : bool or ``"auto"``
            Classic viewer only: store the scene payload as gzip-compressed
            base64. The Oviz viewer always stores compact binary data.
        viewer : {"oviz", "classic"}, optional
            ``"oviz"`` (the default) writes the WebGL2 Oviz viewer: compact
            binary payloads, GPU time interpolation, Views & story, video
            capture. ``"classic"`` writes the previous Three.js runtime
            byte-for-byte (Slides and Paper remain classic-only).
        viewer_mode : {"focus", "detailed"}, optional
            Oviz viewer only: the mode the figure opens in. ``"focus"`` (the
            default) shows the figure, its key and one quiet bar;
            ``"detailed"`` keeps the layers panel, details and every control
            in view. Readers switch modes at any time with **U**.
        camera_anchor : {"lsr", "sun", "free"}, optional
            Oviz viewer only: what the camera orbits, zooms toward and moves
            with through time. By default it is the Local Standard of Rest,
            the origin of the LSR-centred frame (when the home view orbits
            it), so the Sun and the clusters move around a steady camera.
            ``"sun"`` rides along with the Sun trace; ``"free"`` zooms toward
            the pointer. Readers change it under Display settings.
        show_sun : bool, default True
            Oviz viewer only: add the Sun (a one-object "Sun" trace at the
            heliocentric origin, integrated like the clusters) when the data
            have no trace named "Sun". It shows in every layer group and is
            left out in Sky view, which looks out from the Sun.
        spiral_arm_models : sequence, optional
            Published spiral arms to draw (from ``oviz.spiral_models``, e.g.
            ``KHALIL2025_ARMS`` and ``CASTRO_GINARD2021_ARMS``). Each model is
            one line trace holding all of its arms, turning at its pattern
            speeds through the timeline. The traces start hidden (listed in
            every layer group, switched off) and are left out in Sky view.
        show : bool
            Display the figure after construction.
        save_name : path-like, optional
            Write the generated figure to this HTML path.

        Returns
        -------
        OvizFigure or ThreeJSFigure
            The interactive figure for the selected viewer.
        """
        renderer_name = _normalize_renderer_name(renderer)
        self.viewer_name = _normalize_viewer_name(viewer)
        self.viewer_mode = viewer_mode

        time = orbit_maker.normalize_time_grid(time)

        # Determine the reference frame center
        if reference_frame_center is None:
            reference_frame_center = self.set_focus(focus_group)
        # The default frame is centred on the LSR's orbit (orbit_maker), so the
        # viewer can anchor its camera there; a focus group's orbit is not it.
        self.frame_is_lsr = reference_frame_center is None
        self.camera_anchor = camera_anchor
        self.spiral_arm_models = tuple(spiral_arm_models or ())
        # Oviz figures always have the Sun (the classic viewer stays as it was).
        if isinstance(self.data_collection, _CollectionWithSun):
            self.data_collection = self.data_collection.base
        default_sun = self._default_sun_for(bool(show_sun) and self.viewer_name == 'oviz')

        # Re-integrate whenever the requested grid differs from the cached one.
        # Orbit frames and marker sizes are both functions of this timeline.
        cached_time = getattr(self.data_collection, 'time', None)
        needs_integration = cached_time is None
        if not needs_integration:
            try:
                needs_integration = not np.array_equal(
                    np.asarray(cached_time, dtype=float),
                    time,
                )
            except (TypeError, ValueError):
                needs_integration = True
        if needs_integration:
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
        if default_sun is not None:
            self._prepare_default_sun(default_sun, reference_frame_center, fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time)
            self.data_collection = _CollectionWithSun(self.data_collection, default_sun)

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
        self.show_cluster_members_in_sky = bool(show_cluster_members_in_sky)
        self.sky_members_by_cluster = None
        self.volume_configs = _normalize_threejs_volume_configs(volumes)
        profile_initial_state = build_threejs_profile(preset) if preset else {}
        caller_initial_state = copy.deepcopy(threejs_initial_state) if threejs_initial_state else {}
        self.threejs_initial_state = merge_threejs_profile(profile_initial_state, caller_initial_state)
        self.threejs_actions = copy.deepcopy(actions) if actions else []
        self.threejs_compress_scene_spec = compress_scene_spec
        self.threejs_scene_spec_compression_threshold_bytes = scene_spec_compression_threshold_bytes
        lite_mode_enabled = bool(
            self.threejs_initial_state.get("lite_mode_enabled")
            or self.threejs_initial_state.get("minimal_mode_enabled")
        )
        if age_kde_bandwidth_myr <= 0:
            raise ValueError("age_kde_bandwidth_myr must be > 0.")
        self.age_kde_bandwidth_myr = float(age_kde_bandwidth_myr)
        if renderer_name != 'threejs' and self.volume_configs:
            raise NotImplementedError(
                "volumes are currently only supported with renderer='threejs'."
            )
        if renderer_name != 'threejs' and self.threejs_actions:
            raise ValueError("actions are currently only supported with renderer='threejs'.")
        if renderer_name == 'threejs' and self.threejs_actions and not lite_mode_enabled:
            raise ValueError("actions are currently only supported for lite threejs exports.")
        if self.enable_sky_panel:
            if self.sky_radius_deg <= 0:
                raise ValueError("sky_radius_deg must be > 0.")
            if self.sky_frame not in ('galactic', 'icrs'):
                raise ValueError("sky_frame must be either 'galactic' or 'icrs'.")
            if self.sky_frame != 'galactic':
                raise ValueError("sky_frame='icrs' is not yet supported; use 'galactic'.")
            if self.show_cluster_members_in_sky and not cluster_members_file:
                raise ValueError(
                    "show_cluster_members_in_sky=True requires a readable "
                    "cluster_members_file with cluster names and sky coordinates."
                )
        elif self.show_cluster_members_in_sky:
            raise ValueError(
                "show_cluster_members_in_sky=True requires enable_sky_panel=True."
            )

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

        self.figure_layout = _layout(layout_dict)
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
            self._trace_name(st) for st in static_traces if self._trace_name(st)
        ]

        # Build frames for each time step
        cluster_groups = self.data_collection.get_all_clusters()
        frames = []
        present_day_cache = {}

        # Frames and the scene spec are millions of small, acyclic objects.
        with _cyclic_gc_paused():
            for i, t_i in enumerate(self.time):
                traces = self._generate_scatter_list(
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
                    coord_system=self.coord_system,
                    present_day_cache=present_day_cache,
                )
                self._add_static_traces(
                    traces, static_traces, static_traces_times,
                    reference_frame_center, t_i
                )
                # Frames never override the legend group's visibility.
                for trace in traces:
                    trace.pop('visible', None)
                # The traces are new for every frame, so the frame need not copy them.
                frames.append(_frame(data=traces, name=str(t_i)))

            # Initialize the figure at t=0 (base data)
            self._initialize_figure(frames)

            self._build_threejs_figure(frames)

        # Show or save
        if show:
            self.figure.show()
        if save_name:
            self.figure.write_html(save_name)

        if isinstance(self.data_collection, _CollectionWithSun):
            self.data_collection = self.data_collection.base
        return self.figure

    # Frames: per-time traces and legend visibility ---------------------------

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
        coord_system='centered',
        present_day_cache=None,
    ):
        """All traces of the frame at time ``t``: one per cluster group, then the guides.

        ``x_rf, y_rf, z_rf`` is the reference frame's position at ``t``. Outside
        Galactic mode ``show_gc_line`` adds the R = 8.12 kpc circle; in Galactic
        mode the two ``show_galactic_*`` flags control the circles and guides.
        ``present_day_cache`` (a dict) keeps each table's t = 0 sky quantities
        between the frames of one build.
        """
        scatter_list = []
        x_col, y_col, z_col = _xyz_columns(coord_system)

        sun_x = 0.0
        sun_y = 0.0

        for cluster_group in cluster_groups:
            assert cluster_group.integrated
            df_int = cluster_group.df_int
            if df_int.empty:
                continue

            # Rows per time step, indexed once per integrated table instead of
            # rescanning every row for every frame (isclose semantics kept).
            df_t = _rows_at_time(df_int, float(t))
            if df_t.empty:
                continue

            # Present-day (t = 0) sky quantities, for click -> Sky selections.
            # They are the same in every frame, so one build computes them once.
            cached = present_day_cache.get(id(df_int)) if present_day_cache is not None else None
            if cached is None:
                df_t0 = _rows_at_time(df_int, 0.0)
                cached = (len(df_t0), _present_day_sky(df_t0, x_col, y_col, z_col))
                if present_day_cache is not None:
                    present_day_cache[id(df_int)] = cached
            if cached[0] == len(df_t):
                x0, y0, z0, l0, b0, dist0 = cached[1]
            else:
                # Fallback for any unexpected ordering/shape mismatch.
                x0, y0, z0, l0, b0, dist0 = _present_day_sky(df_t, x_col, y_col, z_col)

            if cluster_group.data_name.strip().lower() == 'sun':
                sun_x = float(np.nanmedian(df_t[x_col].to_numpy(dtype=float)))
                sun_y = float(np.nanmedian(df_t[y_col].to_numpy(dtype=float)))

            age_at_t = df_t['age_myr'] + t
            age_present = df_t['age_myr']
            # Literal pieces are joined first: each Series + str is a pass over the column.
            hovertext = (
                '<b style="font-size:16px;">' + df_t['name'].str.replace('_', ' ').astype(str)
                + ('</b><br>' + cluster_group.data_name + '<br>Age (now) = ')
                + age_present.round(1).astype(str)
                + ' Myr<br>Age (t) = '
                + age_at_t.round(1).astype(str)
                + ' Myr<br>'
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
            if 'n_stars' in df_t.columns:
                n_star_values = pd.to_numeric(df_t['n_stars'], errors='coerce').to_numpy(dtype=float)
            else:
                n_star_values = np.full(len(df_t), np.nan, dtype=float)
            if 'name_all' in df_t.columns:
                cluster_alias_values = (
                    df_t['name_all'].fillna('').astype(str).to_numpy(dtype=object)
                )
            else:
                cluster_alias_values = np.repeat('', len(df_t)).astype(object)

            trace_meta = {
                'trace_kind': 'cluster',
                'size_by_n_stars': bool(getattr(cluster_group, 'size_by_n_stars', False)),
            }
            if cluster_group.colormap:
                trace_meta.update({
                    'color_by': 'age',
                    'color_label': 'Age (Myr)',
                    'colormap': str(cluster_group.colormap),
                })

            scatter_list.append(
                _scatter3d(
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
                        trace_colors,
                        n_star_values,
                        cluster_alias_values,
                    )),
                    meta=trace_meta,
                    hovertext=hovertext,
                    hoverinfo='text',  # This removes default x, y, z
                    hovertemplate='%{hovertext}<extra></extra>',  # This ensures only custom hovertext is shown
                    name=cluster_group.data_name
                )
            )

        show_reference_lines = show_gc_line if not galactic_mode else show_galactic_center_circles

        if show_reference_lines:
            if galactic_mode:
                scatter_list.append(
                    self._galactic_center_ring_trace(
                        t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
                    )
                )
                scatter_list.append(
                    self._galactic_center_label_trace(
                        t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
                    )
                )
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
            scatter_list.extend(
                self._build_galactic_guide_traces(
                    sun_x=sun_x,
                    sun_y=sun_y,
                    plane_z_model=self._galactic_plane_z_model(
                        t=t,
                        x_rf=x_rf,
                        y_rf=y_rf,
                        z_rf=z_rf,
                        coord_system=coord_system,
                    ),
                )
            )

        if galactic_mode and include_spiral_arms:
            scatter_list.extend(
                self._build_spiral_arm_traces(
                    t=t, x_rf=x_rf, y_rf=y_rf, z_rf=z_rf, coord_system=coord_system
                )
            )

        for model in getattr(self, 'spiral_arm_models', ()):
            scatter_list.append(
                self._spiral_arm_model_trace(model, t, x_rf, y_rf, z_rf, coord_system=coord_system)
            )

        if self.show_age_kde_inset:
            scatter_list.append(self._build_kde_time_marker_trace(t))
            scatter_list.extend(self._build_kde_inset_traces())

        return scatter_list

    def _add_static_traces(self, traces, static_traces, static_traces_times, reference_frame_center, t):
        """Append each static trace (``meta.static``) to ``traces``, or an empty one when hidden at ``t``.

        ``static_traces_times[i]`` lists the times at which trace ``i`` shows.
        With a focus group, traces other than tracks are recentred on it.
        """
        for i, st in enumerate(static_traces):
            trace_name = self._trace_name(st)
            if t not in static_traces_times[i]:
                traces.append(_scatter3d(
                    x=[], y=[], z=[],
                    name=trace_name,
                    visible=False,
                    meta={'static': True}
                ))
                continue

            st_copy = copy.deepcopy(st)
            existing_meta = st_copy.get('meta') if isinstance(st_copy, dict) else getattr(st_copy, 'meta', None)
            existing_meta = existing_meta if isinstance(existing_meta, dict) else {}
            st_copy['meta'] = {**existing_meta, 'static': True}

            # Re-center if focusing on a group (except for tracks)
            if (self.focus_group is not None) and trace_name and not trace_name.endswith('Track'):
                for axis_idx, axis_key in enumerate(('x', 'y', 'z')):
                    axis_values = st_copy.get(axis_key) if isinstance(st_copy, dict) else getattr(st_copy, axis_key, None)
                    if axis_values is not None:
                        st_copy[axis_key] = np.array(axis_values) - reference_frame_center[axis_idx]
            traces.append(st_copy)

    def _initialize_figure(self, frames):
        """Keep the t = 0 traces, with the first group's visibility, as ``initial_data``."""
        default_group_key = list(self.trace_grouping_dict.keys())[0]  # e.g. "All"
        grouping_0 = self.trace_grouping_dict[default_group_key]

        # A copy of the t = 0 frame's traces, which keep their own visibility.
        idx_zero = np.where(self.time == 0)[0][0]
        self.initial_data = copy.deepcopy(frames[idx_zero]['data'])
        for trace in self.initial_data:
            self._set_trace_visible(trace, self.get_visibility(self._trace_name(trace), grouping_0))

    def _ordered_slider_times(self):
        """Return times in the same order used by the time slider."""
        time_neg = self.time[self.time < 0]
        time_pos = self.time[self.time >= 0]

        if (len(time_neg) > 0) and (len(time_pos) > 1):
            return np.append(time_neg, time_pos)
        if (len(time_neg) > 0) and (len(time_pos) == 1):
            return np.flip(self.time)
        return self.time

    def get_visibility(self, trace_name: str, grouping: list):
        """Initial visibility of a trace in a legend group: True, False or ``"legendonly"``.

        A data trace shows when its name is in ``grouping``, and its orbit
        track (``"<name> Track"``) with it. Other static traces, the default
        Sun and the Galactic guides show in every group (static traces and
        tracks as ``"legendonly"`` with ``static_traces_legendonly``).
        Published spiral-arm models are listed in every group but start hidden.
        """
        if trace_name is None:
            return False

        if trace_name == KDE_TIME_MARKER_TRACE_NAME:
            return True

        if any(trace_name == model.name for model in getattr(self, 'spiral_arm_models', ())):
            return "legendonly"

        kde_source_trace = self._kde_source_trace_name(trace_name)
        if kde_source_trace is not None:
            return kde_source_trace in grouping

        if trace_name.endswith(" Track"):
            if trace_name.replace(" Track", "") not in grouping:
                return False
            return "legendonly" if self.static_traces_legendonly else True
        if trace_name in self.base_static_trace_names:
            return "legendonly" if self.static_traces_legendonly else True

        # The default Sun is in every group, as the Galactic guides are.
        if trace_name == SUN_TRACE_NAME and getattr(self, '_default_sun', None) is not None:
            return True
        if (
            trace_name == 'GC'
            or trace_name in GALACTIC_RADIUS_TRACE_NAMES
            or trace_name in GALACTIC_GUIDE_TRACE_NAMES
            or trace_name in SPIRAL_ARM_TRACE_NAMES
        ):
            return True
        return trace_name in grouping

    def set_focus(self, focus_group):
        """Median ``x, y, z, U, V, W`` of the named trace, or None without a focus group."""
        if not focus_group:
            return None

        focus_group_data = self.data_collection.get_cluster(focus_group).df
        coords = focus_group_data[['x', 'y', 'z', 'U', 'V', 'W']].median().values
        return coords

    def _trace_name(self, trace):
        """The ``name`` of a trace given as a dict or as an object with attributes."""
        if isinstance(trace, dict):
            return trace.get('name')
        return getattr(trace, 'name', None)

    def _set_trace_visible(self, trace, visible_flag):
        """Set ``visible`` on a trace given as a dict or as an object with attributes."""
        if isinstance(trace, dict):
            trace['visible'] = visible_flag
        else:
            trace.visible = visible_flag

    def _default_sun_for(self, wanted):
        """Oviz's default Sun for this plot, or None: only when wanted, for a
        real trace collection (not a test double), and only when the data
        bring no trace named "Sun" of their own."""
        from .traces import TraceCollection

        base = self.data_collection
        if not wanted or not isinstance(base, TraceCollection):
            self._default_sun = None
            return None
        if any(str(getattr(c, 'data_name', '')).strip().lower() == 'sun' for c in base.get_all_clusters()):
            self._default_sun = None
            return None
        if getattr(self, '_default_sun', None) is None:
            self._default_sun = _default_sun_trace()
        return self._default_sun

    def _prepare_default_sun(self, sun, reference_frame_center, fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time):
        """Integrate the default Sun on the collection's time grid and frame
        (again only when either changed) and size it like the other traces."""
        time = np.asarray(self.data_collection.time, dtype=float)
        key = (time.tobytes(), repr(reference_frame_center), self.vo, self.ro, self.zo, id(self.potential))
        if getattr(self, '_default_sun_key', None) != key:
            sun.integrate_orbits(
                time, reference_frame_center=reference_frame_center,
                potential=self.potential, vo=self.vo, ro=self.ro, zo=self.zo,
            )
            self._default_sun_key = key
        if not sun.sizes_set:
            sun.set_age_based_sizes(fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time)

    # Galactic reference geometry ---------------------------------------------

    def _coordFIX_to_coordROT(self, x_gc_pc, y_gc_pc, z_gc_pc, time_myr):
        """Galactocentric pc to the co-rotating frame of ``coord_system='rot'``."""
        return orbit_maker._rotating_frame_xyz(
            x_gc_pc, y_gc_pc, z_gc_pc, time_myr, r_sun=self.ro, v_sun=self.vo
        )

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
        if coord_system == 'rot':
            x_gc_pc, y_gc_pc, z_gc_pc = _radius_circle_xyz_pc(float(radius_kpc), 'galactocentric')
            x_vals, y_vals, z_vals = self._coordFIX_to_coordROT(
                x_gc_pc, y_gc_pc, z_gc_pc, float(t_myr)
            )
        else:
            x_gal, y_gal, z_gal = _radius_circle_xyz_pc(float(radius_kpc), 'galactic')
            x_vals = x_gal - float(x_sub)
            y_vals = y_gal - float(y_sub)
            z_vals = z_gal - float(z_sub)

        line_color = self._reference_line_color()
        return _scatter3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            mode='lines',
            line=dict(
                color=line_color,
                width=self._reference_line_width(),
                dash='solid',
            ),
            opacity=self._reference_opacity(),
            visible=True,
            name=trace_name,
            showlegend=False,
            hovertext=trace_name
        )

    def _reference_line_color(self):
        """Common reference-line color used by GC circle and galactic guides."""
        return '#94a3b8' if self.figure_theme == 'dark' else '#475569'

    def _reference_label_color(self):
        """High-contrast companion color for the Galactic Center annotation."""
        return '#e2e8f0' if self.figure_theme == 'dark' else '#1e293b'

    def _reference_line_width(self):
        """One restrained screen-space stroke for every Galactic guide line."""
        return 1.0

    def _reference_opacity(self):
        """One restrained opacity for circles and quadrant lines."""
        base_opacity = float(getattr(self, 'galactic_reference_opacity', 0.5))
        return min(max(base_opacity * 0.68, 0.0), 1.0)

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

    def _galactic_center_label_trace(
        self,
        t,
        x_rf,
        y_rf,
        z_rf,
        coord_system='centered',
        z_offset_pc=300.0,
    ):
        """Center anchor and label for the Galactic Center."""
        x_gc, y_gc, z_gc = self._galactic_center_position(
            t=t,
            x_rf=x_rf,
            y_rf=y_rf,
            z_rf=z_rf,
            coord_system=coord_system,
        )
        return _scatter3d(
            x=[x_gc],
            y=[y_gc],
            z=[float(z_gc) + float(z_offset_pc)],
            mode='text',
            text=['GALACTIC CENTER'],
            textposition='middle center',
            textfont=dict(
                color=self._reference_label_color(),
                size=18,
                family='Inter, Helvetica Neue, Arial, sans-serif',
            ),
            opacity=0.82,
            name='GC',
            showlegend=False,
            hoverinfo='skip',
            meta={
                'screen_stable_text': True,
                'screen_px': 18.0,
            },
        )

    def _galactic_center_ring_trace(
        self,
        t,
        x_rf,
        y_rf,
        z_rf,
        coord_system='centered',
        radius_pc=250.0,
        npts=96,
    ):
        """Small open-circle ring around the Galactic Center."""
        x_gc, y_gc, z_gc = self._galactic_center_position(
            t=t,
            x_rf=x_rf,
            y_rf=y_rf,
            z_rf=z_rf,
            coord_system=coord_system,
        )
        theta = np.linspace(0.0, 2.0 * np.pi, int(npts), endpoint=True)
        x_vals = float(x_gc) + float(radius_pc) * np.cos(theta)
        y_vals = float(y_gc) + float(radius_pc) * np.sin(theta)
        z_vals = np.full_like(x_vals, float(z_gc))
        line_color = self._reference_line_color()
        return _scatter3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            mode='lines',
            line=dict(color=line_color, width=self._reference_line_width(), dash='solid'),
            opacity=self._reference_opacity(),
            visible=True,
            name='GC Ring',
            showlegend=False,
            hoverinfo='skip'
        )

    def _radius_label_trace(
        self,
        radius_kpc,
        label_text,
        x_center,
        y_center,
        z_center,
        angle_deg=0.0
    ):
        """Place a radius label at a defined angle where 0 deg points +y away from GC."""
        angle_rad = np.deg2rad(float(angle_deg))
        radius_pc = float(radius_kpc) * 1000.0
        x_label = float(x_center) + radius_pc * np.sin(angle_rad)
        y_label = float(y_center) + radius_pc * np.cos(angle_rad)
        z_label = float(z_center)

        return _scatter3d(
            x=[x_label],
            y=[y_label],
            z=[z_label],
            mode='text',
            text=[label_text],
            textposition='middle left',
            textfont=dict(
                color=self._reference_label_color(),
                size=28,
                family='Inter, Helvetica Neue, Arial, sans-serif',
            ),
            opacity=max(self._reference_opacity(), 0.5),
            name=f'{label_text} Label',
            showlegend=False,
            hoverinfo='skip',
            meta={
                'screen_stable_text': True,
                'screen_px': 28.0,
            },
        )

    def _build_galactic_circles_with_labels(self, t, x_rf, y_rf, z_rf, coord_system='centered'):
        """Build a restrained, labelled set of galactocentric radius circles."""
        traces = []
        x_gc, y_gc, z_gc = self._galactic_center_position(
            t=t,
            x_rf=x_rf,
            y_rf=y_rf,
            z_rf=z_rf,
            coord_system=coord_system,
        )
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
                )
            )
        return traces

    def _galactic_plane_z_model(
        self,
        t,
        x_rf,
        y_rf,
        z_rf,
        coord_system='centered',
    ):
        """Return z = ax + by + c for the transformed Galactocentric midplane."""
        if coord_system == 'rot':
            x_plane, y_plane, z_plane = self._coordFIX_to_coordROT(
                np.array([0.0, 1000.0, 0.0], dtype=float),
                np.array([0.0, 0.0, 1000.0], dtype=float),
                np.zeros(3, dtype=float),
                float(t),
            )
        else:
            plane_x, plane_y, plane_z = _galactic_plane_points_pc()
            x_plane = plane_x - float(x_rf)
            y_plane = plane_y - float(y_rf)
            z_plane = plane_z - float(z_rf)
        matrix = np.column_stack((
            np.asarray(x_plane, dtype=float),
            np.asarray(y_plane, dtype=float),
            np.ones(3, dtype=float),
        ))
        try:
            a, b, c = np.linalg.solve(matrix, np.asarray(z_plane, dtype=float))
        except np.linalg.LinAlgError:
            a, b, c = 0.0, 0.0, float(np.nanmedian(z_plane))
        return float(a), float(b), float(c)

    def _build_galactic_guide_traces(
        self,
        sun_x=0.0,
        sun_y=0.0,
        plane_z_model=(0.0, 0.0, 0.0),
    ):
        """Build four simple present-day Galactic quadrant boundaries."""
        guide_color = self._reference_line_color()
        label_color = self._reference_label_color()
        plane_a, plane_b, plane_c = [
            float(value)
            for value in plane_z_model
        ]

        def plane_z(x_value, y_value):
            return (
                plane_a * float(x_value)
                + plane_b * float(y_value)
                + plane_c
            )

        try:
            x_min, x_max = [float(v) for v in self.figure_layout['scene']['xaxis']['range']]
        except Exception:
            x_min, x_max = -10000.0, 10000.0
        try:
            y_min, y_max = [float(v) for v in self.figure_layout['scene']['yaxis']['range']]
        except Exception:
            y_min, y_max = -10000.0, 10000.0
        xy_span = min(x_max - x_min, y_max - y_min)
        ray_radius = 0.68 * np.hypot(x_max - x_min, y_max - y_min)
        major_angles = np.deg2rad([0.0, 90.0, 180.0, 270.0])

        major_x, major_y, major_z = [], [], []
        for angle in major_angles:
            end_x = float(sun_x) + ray_radius * np.cos(angle)
            end_y = float(sun_y) + ray_radius * np.sin(angle)
            major_x.extend([
                float(sun_x),
                end_x,
                None,
            ])
            major_y.extend([
                float(sun_y),
                end_y,
                None,
            ])
            major_z.extend([
                plane_z(sun_x, sun_y),
                plane_z(end_x, end_y),
                None,
            ])

        # The four cardinal longitudes form the quadrant boundaries.
        quadrants = _scatter3d(
            x=major_x,
            y=major_y,
            z=major_z,
            mode='lines',
            line=dict(
                color=guide_color,
                width=self._reference_line_width(),
                dash='solid',
            ),
            opacity=self._reference_opacity(),
            name='Galactic Quadrants',
            showlegend=False,
            hoverinfo='skip'
        )

        # Put the four cardinal longitude labels far enough from the Sun to
        # remain legible over the cluster field. They are intentionally large
        # and screen-stable.
        label_radius = 0.30 * xy_span
        labels = ['ℓ = 0°', 'ℓ = 90°', 'ℓ = 180°', 'ℓ = 270°']
        x_labels = [
            float(sun_x) + label_radius * np.cos(angle)
            for angle in major_angles
        ]
        y_labels = [
            float(sun_y) + label_radius * np.sin(angle)
            for angle in major_angles
        ]

        l_labels = _scatter3d(
            x=x_labels,
            y=y_labels,
            z=[
                plane_z(x_value, y_value)
                for x_value, y_value in zip(x_labels, y_labels)
            ],
            mode='text',
            text=labels,
            textposition='middle center',
            textfont=dict(
                color=label_color,
                size=28,
                family='Inter, Helvetica Neue, Arial, sans-serif',
            ),
            opacity=0.92,
            name='Galactic l Labels',
            showlegend=False,
            hoverinfo='skip',
            meta={
                'screen_stable_text': True,
                'screen_px': 28.0,
            },
        )
        return [quadrants, l_labels]

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
                _scatter3d(
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

    def _spiral_arm_model_trace(self, model, t, x_rf, y_rf, z_rf, coord_system='centered'):
        """One published arm model (``oviz.spiral_models``) at time t: all its arms in one line trace."""
        helio, galcen = spiral_arm_coordinates(model, float(t), ro=self.ro, zo=self.zo)
        if coord_system == 'rot':
            x_vals, y_vals, z_vals = self._coordFIX_to_coordROT(*galcen, float(t))
        else:
            x_vals = helio[0] - float(x_rf)
            y_vals = helio[1] - float(y_rf)
            z_vals = helio[2] - float(z_rf)
        return _scatter3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            mode='lines',
            line=dict(color=model.color, width=float(model.width), dash=model.dash),
            opacity=float(model.opacity),
            name=model.name,
            showlegend=True,
            hoverinfo='skip',
            # Arms in the Galactic plane are 3D context; Sky view leaves them out.
            meta={'oviz_hide_in_sky': True},
        )

    # Age KDE inset -----------------------------------------------------------

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
        """Precompute one age KDE per trace and add the inset's axes to ``layout_dict``."""
        self.kde_density_by_trace = {}
        self.kde_trace_name_by_trace = {}

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

        # The KDE x-range matches the timeline's span.
        self.kde_x_grid = np.linspace(x_min, x_max, 300)
        self.kde_x_range = [float(x_min), float(x_max)]

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

        axis_color = 'gray' if self.figure_theme == 'dark' else 'black'
        self.kde_axis_color = axis_color

        # Centre the inset horizontally.
        x2_domain = [0.30, 0.70]
        y2_domain = [0.03, 0.17]
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
            title=dict(
                text='',
                font=dict(color=axis_color, size=10, family='helvetica'),
            ),
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
            title=dict(
                text='Relative SFH',
                font=dict(color=axis_color, size=10, family='helvetica'),
            ),
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
                _scatter(
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
        return _scatter(
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

    # Scene spec for the viewers ----------------------------------------------

    def _build_threejs_figure(self, frames):
        """Build the standalone figure wrapper from the current frame data."""
        scene_spec = self._build_threejs_scene_spec(frames)
        self.fig = scene_spec
        self.fig_dict = scene_spec
        if getattr(self, "viewer_name", DEFAULT_VIEWER) == "oviz":
            from .viewer.figure import OvizFigure

            self.figure = OvizFigure(
                scene_spec,
                mode=getattr(self, "viewer_mode", None) or "focus",
                camera_anchor=getattr(self, "camera_anchor", None),
                lsr_origin=getattr(self, "frame_is_lsr", True),
            )
            return
        self.figure = ThreeJSFigure(
            scene_spec,
            compress_scene_spec=getattr(self, "threejs_compress_scene_spec", "auto"),
            scene_spec_compression_threshold_bytes=getattr(
                self,
                "threejs_scene_spec_compression_threshold_bytes",
                None,
            ),
        )

    def _build_threejs_scene_spec(self, frames):
        """Serialize the current animation into a renderer-agnostic scene spec."""
        # Every trace of every frame offers the same colormaps: sample them once per build.
        self._colormap_options_cache = {}
        try:
            return build_threejs_scene_spec(
                self,
                frames,
                trace_to_scene_json=_trace_to_scene_json,
                coerce_range=_coerce_range,
                format_time_label=_format_time_label,
                coerce_float=_coerce_float,
                file_to_data_url=_threejs_file_to_data_url,
                catalog_from_frame_spec=_threejs_catalog_from_frame_spec,
                annotate_point_motion_ranges=_annotate_threejs_point_motion_ranges,
            )
        finally:
            self._colormap_options_cache = None

    def _build_threejs_sky_panel_spec(self, default_catalog=None):
        """Sky view settings and the member-star catalog per cluster."""
        if not getattr(self, 'enable_sky_panel', False):
            return {'enabled': False}
        if (
            getattr(self, 'show_cluster_members_in_sky', False)
            and self.sky_members_by_cluster is None
        ):
            # Full Hunt member catalogs can contain more than a million rows.
            # Read them lazily and retain only clusters represented by this
            # figure (including their aliases) before embedding the HTML.
            requested_cluster_names = set((default_catalog or {}).keys())
            self.sky_members_by_cluster = _load_threejs_cluster_catalog(
                self.cluster_members_file,
                cluster_names=requested_cluster_names,
            )
            if not self.sky_members_by_cluster:
                raise ValueError(
                    "show_cluster_members_in_sky=True requires a readable "
                    "cluster_members_file with members matching figure clusters."
                )
        merged_catalog = _merge_threejs_member_catalogs(default_catalog, self.sky_members_by_cluster)

        return {
            'enabled': True,
            'radius_deg': float(self.sky_radius_deg),
            'frame': self.sky_frame,
            'survey': self.sky_survey,
            'show_cluster_members_in_sky': bool(
                getattr(self, 'show_cluster_members_in_sky', False)
            ),
            'member_point_size_denominator': 30,
            # A 1/20-scale point can project to substantially less than one
            # physical pixel in an all-sky view. Keep the requested data-space
            # scale while enforcing a small render-only visibility floor.
            'member_min_screen_size_px': 0.0,
            'member_distance_policy': 'stellar_parallax_preferred',
            'members_by_cluster': merged_catalog,
        }

    def _build_threejs_age_kde_spec(self, trace_key_by_name=None):
        """The age KDE widget: one curve per trace and every cluster's present age."""
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
            data_df = _cluster_table(cluster_group)
            if data_df is None or ('age_myr' not in data_df.columns):
                continue

            df_present = _present_day_rows(data_df)
            age_values = pd.to_numeric(df_present['age_myr'], errors='coerce').to_numpy(dtype=float)
            cluster_names = _row_names(df_present, trace_name)

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
        """The cluster filter widget: one entry per cluster with its age and star count."""
        trace_key_by_name = trace_key_by_name or {}
        entries = []
        seen_keys = set()

        for cluster_group in self.data_collection.get_all_clusters():
            trace_name = str(getattr(cluster_group, 'data_name', '') or '')
            trace_key = trace_key_by_name.get(trace_name)
            data_df = _cluster_table(cluster_group)
            if data_df is None or data_df.empty or ('age_myr' not in data_df.columns):
                continue

            df_present = _present_day_rows(data_df)
            if df_present.empty:
                continue

            age_values = pd.to_numeric(df_present['age_myr'], errors='coerce').to_numpy(dtype=float)
            if 'n_stars' in df_present.columns:
                n_stars_values = pd.to_numeric(df_present['n_stars'], errors='coerce').to_numpy(dtype=float)
            else:
                n_stars_values = np.full(len(df_present), np.nan, dtype=float)
            cluster_names = _row_names(df_present, trace_name)

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
        """The Birth tree widget: each cluster's birth position and sampled track."""
        trace_key_by_name = trace_key_by_name or {}
        entries = []
        trace_options = []
        x_col, y_col, z_col = _xyz_columns(getattr(self, 'coord_system', 'centered'))

        for cluster_group in self.data_collection.get_all_clusters():
            trace_name = str(getattr(cluster_group, 'data_name', '') or '')
            trace_key = trace_key_by_name.get(trace_name)
            if not trace_name or not trace_key:
                continue

            data_df = _cluster_table(cluster_group)
            if data_df is None or data_df.empty or ('age_myr' not in data_df.columns):
                continue

            if not {x_col, y_col, z_col}.issubset(data_df.columns):
                continue

            cluster_color = getattr(cluster_group, 'color', None)
            if not isinstance(cluster_color, str) or not cluster_color:
                cluster_color = '#ffffff'

            # Rows per cluster, grouped as ``groupby(name.astype(str), sort=False)``
            # does (missing names dropped), read from whole-column arrays.
            if 'name' in data_df.columns:
                keys = data_df[[]].copy()
                keys['__cluster_name'] = data_df['name'].astype(str)
                cluster_rows = keys.groupby('__cluster_name', sort=False).indices
            else:
                cluster_rows = {trace_name: np.arange(len(data_df))}
            if 'time' in data_df.columns:
                time_values = pd.to_numeric(data_df['time'], errors='coerce').to_numpy(dtype=float)
            else:
                time_values = np.zeros(len(data_df))
            age_values_all = pd.to_numeric(data_df['age_myr'], errors='coerce').to_numpy(dtype=float)
            x_values = pd.to_numeric(data_df[x_col], errors='coerce').to_numpy(dtype=float)
            y_values = pd.to_numeric(data_df[y_col], errors='coerce').to_numpy(dtype=float)
            z_values = pd.to_numeric(data_df[z_col], errors='coerce').to_numpy(dtype=float)

            trace_entry_count = 0
            trace_age_values = age_values_all[np.isfinite(age_values_all)]
            trace_max_age = float(np.nanmax(trace_age_values)) if trace_age_values.size else 0.0

            for cluster_name, rows in cluster_rows.items():
                age_values = age_values_all[rows]
                age_values = age_values[np.isfinite(age_values)]
                if age_values.size == 0:
                    continue
                age_now = float(np.nanmax(age_values))
                birth_time_myr = -age_now

                time_samples = time_values[rows]
                x_samples = x_values[rows]
                y_samples = y_values[rows]
                z_samples = z_values[rows]
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
                    'time_samples': time_samples.tolist(),
                    'x_samples': x_samples.tolist(),
                    'y_samples': y_samples.tolist(),
                    'z_samples': z_samples.tolist(),
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
            'default_threshold_mode': 'distance_pc',
            'default_connection_mode': 'birth_to_older_track',
            'default_threshold_pc': 100.0,
            'threshold_min_pc': 0.0,
            'threshold_max_pc': max(5000.0, max_age * 50.0 if np.isfinite(max_age) else 5000.0),
            'default_threshold_age_myr': 5.0,
            'threshold_min_age_myr': 0.0,
            'threshold_max_age_myr': max(10.0, max_age + 5.0 if np.isfinite(max_age) else 10.0),
            'max_age_myr': max_age,
            'traces': trace_options,
            'entries': entries,
        }

    def _build_threejs_volume_layers(self):
        """Volume layer specs, centred on the reference frame at t = 0."""
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
        return None

    def _threejs_frame_decorations(
        self,
        frame_json,
        time_value,
        x_range,
        y_range,
        z_range,
        fallback_center,
        volume_layers=None,
        galactic_simple_config=None,
        galaxy_image_config=None,
    ):
        """Frame-local decorations: the volume layers shown at ``time_value``,
        galaxy image planes, the Galactic-lite guides and the Milky Way model."""
        at_present = np.isclose(float(time_value), 0.0, atol=1e-9)
        decorations = []
        for layer in volume_layers or []:
            layer_time = _coerce_threejs_volume_time_myr(layer.get('time_myr'))
            if layer_time is not None:
                if not np.isclose(float(time_value), float(layer_time), atol=1e-9):
                    continue
            elif (
                layer.get('supports_show_all_times', False)
                and layer.get('time_myr') in (None, '', False)
            ):
                pass
            elif layer.get('only_at_t0', True) and not at_present:
                continue
            decorations.append({
                'kind': 'volume_layer',
                'key': layer.get('key'),
                'state_key': layer.get('state_key') or layer.get('key'),
            })

        if galaxy_image_config and galaxy_image_config.get('enabled'):
            image_size_pc = float(_coerce_float(galaxy_image_config.get('size_pc'), 40000.0))
            plane_center = self._threejs_milky_way_center(frame_json, fallback_center)
            if galaxy_image_config.get('image_data_url') and image_size_pc > 0.0:
                if bool(galaxy_image_config.get('only_at_t0', True)):
                    fade_alpha = 1.0 if at_present else 0.0
                else:
                    fade_alpha = 1.0
                decorations.append(
                    _image_plane_decoration(galaxy_image_config, 'galaxy-image-overlay', plane_center, fade_alpha)
                )

        if galactic_simple_config and galactic_simple_config.get('enabled'):
            image_size_pc = float(_coerce_float(galactic_simple_config.get('size_pc'), 40000.0))
            plane_center = self._threejs_milky_way_center(frame_json, fallback_center)
            if galactic_simple_config.get('image_data_url') and image_size_pc > 0.0:
                decorations.append(_image_plane_decoration(
                    galactic_simple_config, 'galactic-plane-overlay', plane_center, 1.0 if at_present else 0.0
                ))

            xy_span = max(min(float(x_range[1]) - float(x_range[0]), float(y_range[1]) - float(y_range[0])), 1.0)
            z_span = max(float(z_range[1]) - float(z_range[0]), 1.0)
            circle_radius_pc = min(8122.0, 0.48 * xy_span)
            decorations.append({
                'kind': 'galactic_center_axes',
                'key': 'galactic-center-axes',
                'center': {
                    'x': float(plane_center.get('x', 0.0)),
                    'y': float(plane_center.get('y', 0.0)),
                    'z': 0.0,
                },
                'half_length_xy': circle_radius_pc,
                'half_length_z': 0.42 * z_span,
                'color': '#7a8089',
                'opacity': 0.62,
                'width_px': 1.0,
                'render_order': -8,
            })

            sun_center = self._threejs_sun_position(frame_json, fallback_center)
            if at_present:
                decorations.append({
                    'kind': 'solar_system_marker',
                    'key': 'solar-system-marker',
                    'base': {
                        'x': float(sun_center.get('x', 0.0)),
                        'y': float(sun_center.get('y', 0.0)),
                        'z': 0.0,
                    },
                    'bottom_z': -140.0,
                    'top_z': 760.0,
                    'label': 'Solar System',
                    'label_z': 900.0,
                    'color': '#ffe45c',
                    'hide_below_scale_bar_pc': float(_coerce_float(galactic_simple_config.get('hide_below_scale_bar_pc'), 400.0)),
                    'fade_start_scale_bar_pc': float(_coerce_float(galactic_simple_config.get('fade_start_scale_bar_pc'), 1000.0)),
                    'render_order': 8,
                })

        if not getattr(self, 'show_milky_way_model', False) or not at_present:
            return decorations

        x_span = max(float(x_range[1]) - float(x_range[0]), 1.0)
        y_span = max(float(y_range[1]) - float(y_range[0]), 1.0)
        z_span = max(float(z_range[1]) - float(z_range[0]), 1.0)
        disc_radius_pc = float(np.clip(0.49 * min(x_span, y_span), 4500.0, 14000.0))
        center = self._threejs_milky_way_center(frame_json, fallback_center)

        decorations.append({
            'kind': 'milky_way_model',
            'opacity_scale': 1.0,
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
        """The plotted Galactic-centre marker of the frame, else ``fallback_center``."""
        point = _first_trace_point(frame_json, 'GC')
        if point is not None:
            return point
        return {
            'x': float(fallback_center['x']),
            'y': float(fallback_center['y']),
            'z': float(fallback_center['z']),
        }

    def _threejs_sun_position(self, frame_json, fallback_center):
        """The plotted Sun of the frame, else the origin at the fallback height."""
        point = _first_trace_point(frame_json, 'Sun')
        if point is not None:
            return point
        return {
            'x': 0.0,
            'y': 0.0,
            'z': float(fallback_center.get('z', 0.0)),
        }

    def _threejs_theme(self, layout_json):
        """Translate the layout theme into a smaller scene theme."""
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

    def _threejs_trace_spec(self, trace_key, trace_json, minimal_mode=False, galactic_simple_mode=False):
        """Convert a trace dictionary into a simpler primitive spec for three.js."""
        if trace_json.get('type', 'scatter3d') != 'scatter3d':
            return None

        trace_name = str(trace_json.get('name') or '')
        if galactic_simple_mode and trace_name not in _galactic_simple_allowed_trace_names(self):
            return None
        if galactic_simple_mode and trace_name in (GALACTIC_RADIUS_TRACE_NAMES - {'R = 8.12 kpc'}):
            return None

        mode = str(trace_json.get('mode', 'markers'))
        trace_meta = trace_json.get('meta') if isinstance(trace_json.get('meta'), dict) else {}
        spec = {
            'key': trace_key,
            'name': trace_json.get('name') or trace_key,
            'showlegend': bool(trace_json.get('showlegend', True)),
            'size_by_n_stars_default': bool(trace_meta.get('size_by_n_stars')),
        }
        # meta={"oviz_role": "guide"}: an annotation (labels, a length bar, a
        # model curve) that the Oviz viewer lists under Guides, not as data.
        if str(trace_meta.get('oviz_role') or '').strip().lower() == 'guide':
            spec['role'] = 'guide'
        # meta={"oviz_pickable": False}: points that sample a surface or a
        # model rather than objects; drawn, but never clicked or lassoed.
        if trace_meta.get('oviz_pickable') is False:
            spec['pickable'] = False
        # meta={"oviz_hide_in_sky": True}: 3D-only context (a shell around the
        # Sun, say) that Sky view, looking out from the Sun, leaves out.
        if trace_meta.get('oviz_hide_in_sky') is True:
            spec['hide_in_sky'] = True
        legend_color = None

        if 'lines' in mode:
            segments = _line_segments_from_trace(trace_json)
            if segments:
                line_color, line_opacity = _color_to_css_and_opacity(
                    trace_json.get('line', {}).get('color'),
                    base_opacity=float(trace_json.get('opacity', 1.0)),
                )
                spec['segments'] = segments
                line_width = float(trace_json.get('line', {}).get('width', 1.0) or 1.0)
                if galactic_simple_mode and trace_name == 'R = 8.12 kpc':
                    line_width = 1.15
                spec['line'] = {
                    'color': line_color,
                    'width': line_width,
                    'dash': trace_json.get('line', {}).get('dash', 'solid'),
                }
                spec['opacity'] = line_opacity
                legend_color = line_color

        if 'markers' in mode:
            points = _points_from_trace(
                trace_json,
                default_opacity=float(trace_json.get('opacity', 1.0)),
                include_selection=not minimal_mode,
                include_hovertext=not minimal_mode,
                include_motion=(not minimal_mode) or galactic_simple_mode,
                include_n_stars=not minimal_mode,
            )
            if points:
                if galactic_simple_mode:
                    point_size_boost = 1.0
                    if trace_name == 'Sun':
                        point_size_boost = 1.35
                    elif trace_name == 'Clusters (< 60 Myr)':
                        point_size_boost = 1.18
                    if not np.isclose(point_size_boost, 1.0):
                        for point in points:
                            point_size = _coerce_float(point.get('size'), math.nan)
                            if math.isfinite(point_size) and point_size > 0.0:
                                point['size'] = float(point_size * point_size_boost)
                spec['points'] = points
                point_sizes = [
                    float(point.get('size'))
                    for point in points
                    if math.isfinite(point.get('size')) and float(point.get('size')) > 0.0
                ]
                if point_sizes:
                    spec['default_point_size'] = float(np.median(point_sizes))
                spec['has_n_stars'] = any(
                    math.isfinite(_coerce_float(point.get('n_stars'), math.nan)) for point in points
                ) if not minimal_mode else False
                point_opacities = [
                    float(point.get('opacity'))
                    for point in points
                    if math.isfinite(point.get('opacity'))
                ]
                if point_opacities:
                    spec['default_opacity'] = float(np.median(point_opacities))
                if legend_color is None:
                    legend_color = points[0].get('color')
                color_by = _threejs_trace_color_by_spec(
                    trace_json, points, colormap_cache=getattr(self, '_colormap_options_cache', None)
                )
                if color_by is not None:
                    spec['color_by'] = color_by
                    if color_by.get('default_color_mode') == 'by_value':
                        legend_color = color_by.get('legend_color') or legend_color

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


# Theme and orbit tracks ------------------------------------------------------


def read_theme(plot):
    """
    Reads the theme configuration from a YAML file and sets up figure layout based on
    the provided figure_theme (e.g., 'light', 'dark', etc.).
    """
    theme_resource = importlib.resources.files("oviz.themes").joinpath(f"{plot.figure_theme}.yaml")
    with theme_resource.open("r", encoding="utf-8") as file:
        layout = yaml.safe_load(file)

    if plot.figure_theme == 'light':
        layout['template'] = (
            plot.light_template if plot.light_template else 'default'
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
    x_col, y_col, z_col = _xyz_columns(coord_system)

    tracks = _scatter3d(
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


# Small helpers ---------------------------------------------------------------


@contextlib.contextmanager
def _cyclic_gc_paused():
    """Pause the cyclic garbage collector, restoring its previous state afterwards.

    Building a scene creates millions of small dicts and lists without
    reference cycles; full collections would only traverse them again.
    Reference counting still frees everything as usual.
    """
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if was_enabled:
            gc.enable()


def _normalize_viewer_name(viewer):
    """``"oviz"`` (the default) or ``"classic"`` from the names ``make_plot`` accepts."""
    if viewer is None:
        return DEFAULT_VIEWER
    name = str(viewer).strip().lower()
    if name in ("oviz", "new", "v2", "webgl2", "default"):
        return "oviz"
    if name in ("classic", "legacy", "threejs", "three.js"):
        return "classic"
    raise ValueError("viewer must be 'oviz' or 'classic'.")


def _normalize_renderer_name(renderer):
    """``"threejs"``, the only renderer; anything else raises ``ValueError``."""
    renderer_name = str(renderer).strip().lower()
    if renderer_name in ('three', 'threejs', 'three.js'):
        return 'threejs'
    raise ValueError("renderer must be 'threejs'.")


def _trace_to_scene_json(trace):
    """A deep copy of a trace dict (or of a record with ``to_scene_json``)."""
    if isinstance(trace, dict):
        return copy.deepcopy(trace)
    if hasattr(trace, 'to_scene_json'):
        return trace.to_scene_json()
    raise TypeError(f'Unsupported trace type for serialization: {type(trace)!r}')


def _coerce_range(values, default):
    """``[low, high]`` floats from a two-item sequence, else ``default``."""
    if not isinstance(values, (list, tuple, np.ndarray)) or len(values) != 2:
        return [float(default[0]), float(default[1])]
    try:
        return [float(values[0]), float(values[1])]
    except Exception:
        return [float(default[0]), float(default[1])]


def _format_time_label(time_value):
    """Frame label for a time in Myr: an integer when whole, else one decimal."""
    rounded = round(float(time_value), 10)
    if np.isclose(rounded, round(rounded), atol=1e-9):
        return str(int(round(rounded)))
    return f'{rounded:.1f}'.rstrip('0').rstrip('.')


def _is_sequence_value(value):
    return isinstance(value, (list, tuple, np.ndarray, pd.Series))


def _expand_value(value, length, default=None):
    """A list of ``length`` values from a scalar or a sequence (padded with its last item or cut)."""
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
    """``float(value)`` when finite, else ``float(default)``."""
    try:
        out = float(value)
        if math.isfinite(out):
            return out
    except Exception:
        pass
    return float(default)


def _as_object_list(value):
    """``value`` as a flat Python list (None gives an empty list, a scalar one item)."""
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


def _threejs_file_to_data_url(path_value):
    """An image file as a base64 ``data:`` URL, or None when there is no such file."""
    if not path_value:
        return None
    path = Path(path_value).expanduser()
    if not path.exists() or not path.is_file():
        return None
    suffix = path.suffix.lower()
    mime_type = {
        '.jpg': 'image/jpeg',
        '.jpeg': 'image/jpeg',
        '.png': 'image/png',
        '.webp': 'image/webp',
        '.gif': 'image/gif',
    }.get(suffix, 'application/octet-stream')
    encoded = base64.b64encode(path.read_bytes()).decode('ascii')
    return f'data:{mime_type};base64,{encoded}'


# Integrated tables -----------------------------------------------------------


def _present_day_sky(rows, x_col, y_col, z_col):
    """Positions ``x, y, z`` and heliocentric ``l``, ``b`` (deg) and distance (pc) of table rows."""
    x0 = pd.to_numeric(rows[x_col], errors='coerce').to_numpy(dtype=float)
    y0 = pd.to_numeric(rows[y_col], errors='coerce').to_numpy(dtype=float)
    z0 = pd.to_numeric(rows[z_col], errors='coerce').to_numpy(dtype=float)

    x_helio0 = pd.to_numeric(rows['x_helio'], errors='coerce').to_numpy(dtype=float)
    y_helio0 = pd.to_numeric(rows['y_helio'], errors='coerce').to_numpy(dtype=float)
    z_helio0 = pd.to_numeric(rows['z_helio'], errors='coerce').to_numpy(dtype=float)
    dist0 = np.sqrt(x_helio0 ** 2 + y_helio0 ** 2 + z_helio0 ** 2)

    with np.errstate(invalid='ignore', divide='ignore'):
        l0 = np.rad2deg(np.arctan2(y_helio0, x_helio0))
        l0 = np.mod(l0, 360.0)
        b0 = np.rad2deg(np.arcsin(np.clip(z_helio0 / np.where(dist0 > 0, dist0, np.nan), -1.0, 1.0)))
    return x0, y0, z0, l0, b0, dist0


def _xyz_columns(coord_system):
    """Position columns of an integrated table for ``coord_system`` (``'rot'`` or centred)."""
    if coord_system == 'rot':
        return 'x_rot', 'y_rot', 'z_rot'
    return 'x', 'y', 'z'


def _cluster_table(cluster_group):
    """A trace's integrated table, or its input table before integration (None without either)."""
    if getattr(cluster_group, 'df_int', None) is not None:
        return cluster_group.df_int
    return getattr(cluster_group, 'df', None)


def _present_day_rows(data_df):
    """The rows at t = 0 when the table has a ``time`` column with any; otherwise all rows."""
    if 'time' in data_df.columns:
        time_values = pd.to_numeric(data_df['time'], errors='coerce').to_numpy(dtype=float)
        zero_mask = np.isclose(time_values, 0.0, rtol=0.0, atol=1e-9)
        if np.any(zero_mask):
            return data_df.loc[zero_mask]
    return data_df


def _row_names(data_df, trace_name):
    """Each row's ``name`` as a string, or ``trace_name`` for tables without names."""
    if 'name' in data_df.columns:
        return data_df['name'].astype(str).to_numpy(dtype=object)
    return np.array([trace_name] * len(data_df), dtype=object)


_TIME_INDEX_CACHE: "dict[int, tuple]" = {}


def _rows_at_time(df_int, t):
    """``df_int[np.isclose(df_int['time'], t, atol=1e-9)]`` in O(rows at t).

    The index is rebuilt whenever the table object or its time column
    changes, so re-integrated traces never see stale rows.
    """

    times = df_int['time'].to_numpy(dtype=float)
    key = id(df_int)
    cached = _TIME_INDEX_CACHE.get(key)
    if cached is None or cached[0] is not df_int or cached[1] != len(times) or not np.array_equal(cached[2], times):
        order = np.argsort(times, kind='stable')
        sorted_times = times[order]
        cached = (df_int, len(times), times.copy(), order, sorted_times)
        if len(_TIME_INDEX_CACHE) >= 64:
            _TIME_INDEX_CACHE.clear()
        _TIME_INDEX_CACHE[key] = cached
    _, _, _, order, sorted_times = cached
    lo = np.searchsorted(sorted_times, t - 1e-9, side='left')
    hi = np.searchsorted(sorted_times, t + 1e-9, side='right')
    rows = np.sort(order[lo:hi])
    return df_int.iloc[rows]


# Galactic geometry and sky coordinates ---------------------------------------


@functools.lru_cache(maxsize=32)
def _radius_circle_xyz_pc(radius_kpc, frame):
    """Cartesian pc of a 1000-point galactocentric circle, in ``frame``.

    The geometry is identical for every timeline frame (only the per-frame
    recentring differs), so the astropy transform runs once per radius.
    """

    n_marks = 1000
    R = -float(radius_kpc) * np.ones(n_marks)
    phi = np.linspace(-180, 180, n_marks)
    circle = SkyCoord(
        rho=R * u.kpc,
        phi=phi * u.deg,
        z=[0.] * len(R) * u.pc,
        frame='galactocentric',
        representation_type='cylindrical'
    )
    target = circle.galactocentric if frame == 'galactocentric' else circle.galactic
    xyz = [np.asarray(c.value, dtype=float) * 1000.0 for c in (target.cartesian.x, target.cartesian.y, target.cartesian.z)]
    for arr in xyz:
        arr.setflags(write=False)
    return tuple(xyz)


@functools.lru_cache(maxsize=1)
def _galactic_plane_points_pc():
    """Galactic Cartesian pc of three points spanning the Galactocentric midplane.

    The points (the Galactic centre and 1 kpc along x and y) never change,
    so the astropy transform runs once instead of once per frame.
    """
    plane = SkyCoord(
        x=np.array([0.0, 1000.0, 0.0]) * u.pc,
        y=np.array([0.0, 0.0, 1000.0]) * u.pc,
        z=np.zeros(3) * u.pc,
        frame='galactocentric',
        representation_type='cartesian',
    ).galactic.cartesian
    xyz = tuple(np.array(axis.to_value(u.pc), dtype=float) for axis in (plane.x, plane.y, plane.z))
    for arr in xyz:
        arr.setflags(write=False)
    return xyz


@functools.lru_cache(maxsize=1)
def _galactic_to_icrs_matrix():
    """Astropy's Galactic → ICRS rotation, evaluated once.

    The chain (Galactic → FK5 J2000 → ICRS) is a fixed rotation, so
    transforming the three basis vectors reproduces astropy to ~1e-15
    without building a SkyCoord per point (which dominated make_plot).
    """

    # Images of the Galactic x, y, z unit vectors are the matrix columns.
    basis = SkyCoord(l=[0.0, 90.0, 0.0] * u.deg, b=[0.0, 0.0, 90.0] * u.deg, frame='galactic')
    return np.asarray(basis.icrs.cartesian.xyz.value, dtype=float)


@functools.lru_cache(maxsize=1 << 16)
def _galactic_to_icrs_deg(l_deg, b_deg):
    """(l, b) → (ra, dec) in degrees; NaN pair when the input is invalid."""

    if not (np.isfinite(l_deg) and np.isfinite(b_deg)):
        return (np.nan, np.nan)
    try:
        matrix = _galactic_to_icrs_matrix()
    except Exception:
        return (np.nan, np.nan)
    l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)
    cb = np.cos(b_rad)
    v = matrix @ np.array([cb * np.cos(l_rad), cb * np.sin(l_rad), np.sin(b_rad)])
    ra = float(np.degrees(np.arctan2(v[1], v[0])) % 360.0)
    dec = float(np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0))))
    return (ra, dec)


# Colours ---------------------------------------------------------------------


_COLOR_SCALE_STOPS = {
    "viridis": ("#440154", "#31688e", "#35b779", "#fde725"),
    "plasma": ("#0d0887", "#9c179e", "#ed7953", "#f0f921"),
    "magma": ("#000004", "#51127c", "#b73779", "#fcfdbf"),
    "inferno": ("#000004", "#57106e", "#bc3754", "#fcffa4"),
    "cividis": ("#00224e", "#575d6d", "#a59c74", "#fee838"),
    "turbo": ("#30123b", "#28a5f5", "#7ef658", "#fca636", "#7a0403"),
    "greys": ("#000000", "#777777", "#ffffff"),
    "gist_heat": ("#000000", "#b00000", "#ffff00", "#ffffff"),
}


def _hex_to_rgb_tuple(value):
    """``(r, g, b)`` of a ``#rgb`` or ``#rrggbb`` colour; mid-grey when unreadable."""
    value = str(value or "").strip().lstrip("#")
    if len(value) == 3:
        value = "".join(ch * 2 for ch in value)
    if len(value) != 6:
        return (128, 128, 128)
    try:
        return tuple(int(value[idx : idx + 2], 16) for idx in (0, 2, 4))
    except ValueError:
        return (128, 128, 128)


def _sample_colorscale(colorscale, positions):
    """``rgb(...)`` strings of a named or listed colour scale at ``positions`` in [0, 1]."""
    if isinstance(colorscale, str):
        stops = _COLOR_SCALE_STOPS.get(colorscale.strip().lower(), _COLOR_SCALE_STOPS["viridis"])
    elif isinstance(colorscale, (list, tuple)) and colorscale:
        raw_stops = []
        for item in colorscale:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                raw_stops.append(item[1])
            else:
                raw_stops.append(item)
        stops = tuple(str(item) for item in raw_stops) or _COLOR_SCALE_STOPS["viridis"]
    else:
        stops = _COLOR_SCALE_STOPS["viridis"]

    rgbs = [_hex_to_rgb_tuple(stop) for stop in stops]
    if len(rgbs) == 1:
        rgbs = [rgbs[0], rgbs[0]]
    out = []
    for position in positions:
        t = min(max(float(position), 0.0), 1.0)
        scaled = t * (len(rgbs) - 1)
        lo = int(np.floor(scaled))
        hi = min(lo + 1, len(rgbs) - 1)
        frac = scaled - lo
        rgb = tuple(
            int(round(rgbs[lo][channel] + (rgbs[hi][channel] - rgbs[lo][channel]) * frac))
            for channel in range(3)
        )
        out.append(f"rgb({rgb[0]}, {rgb[1]}, {rgb[2]})")
    return out


def _color_to_css_and_opacity(color, base_opacity=1.0):
    """A CSS colour and an opacity (``base_opacity`` times any alpha) from a colour value."""
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
    """Per-point CSS colours and alphas of a marker: fixed, listed, or numbers on its colour scale."""
    color_value = marker.get('color')
    if color_value is None:
        css, alpha = _color_to_css_and_opacity('#808080', default_opacity)
        return [css] * length, [alpha] * length

    numeric_values = _marker_numeric_color_values(marker, length)
    if numeric_values is not None:
        colorscale = marker.get('colorscale', 'Viridis')
        cmin = _coerce_float(marker.get('cmin'), np.nan)
        cmax = _coerce_float(marker.get('cmax'), np.nan)
        finite_mask = np.isfinite(numeric_values)
        if not np.isfinite(cmin):
            cmin = float(np.nanmin(numeric_values[finite_mask])) if np.any(finite_mask) else 0.0
        if not np.isfinite(cmax):
            cmax = float(np.nanmax(numeric_values[finite_mask])) if np.any(finite_mask) else 1.0
        scaled = np.zeros_like(numeric_values, dtype=float)
        if not np.isclose(cmax, cmin, atol=1e-12):
            scaled[finite_mask] = np.clip((numeric_values[finite_mask] - cmin) / (cmax - cmin), 0.0, 1.0)
        css_values = ['#808080'] * length
        alpha_values = [float(default_opacity)] * length
        if np.any(finite_mask):
            finite_indices = np.flatnonzero(finite_mask)
            sampled = _sample_colorscale(colorscale, scaled[finite_mask].tolist())
            for point_idx, item in zip(finite_indices, sampled):
                css, alpha = _color_to_css_and_opacity(item, default_opacity)
                css_values[int(point_idx)] = css
                alpha_values[int(point_idx)] = alpha
        return css_values, alpha_values

    if _is_sequence_value(color_value):
        color_values = list(np.asarray(color_value, dtype=object).tolist())
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


def _marker_numeric_color_values(marker, length):
    """The marker's colours as numbers (NaN for blanks), or None unless all are numeric and some finite."""
    color_value = marker.get('color') if isinstance(marker, dict) else None
    if color_value is None or length <= 0:
        return None

    if _is_sequence_value(color_value):
        values = _as_object_list(color_value)
        if len(values) not in (1, length):
            return None
        expanded = _expand_value(values, length, np.nan)
    else:
        expanded = [color_value] * length

    numeric_values = []
    has_finite = False
    for value in expanded:
        if value in (None, ''):
            numeric_values.append(np.nan)
            continue
        try:
            numeric_value = float(value)
        except Exception:
            return None
        if np.isfinite(numeric_value):
            has_finite = True
            numeric_values.append(float(numeric_value))
        else:
            numeric_values.append(np.nan)

    if not has_finite:
        return None
    return np.asarray(numeric_values, dtype=float)


def _threejs_trace_color_by_spec(trace_json, points, colormap_cache=None):
    """The trace's colour-by-value settings (scale, range, colormaps), or None without scalars.

    ``colormap_cache`` (a dict) reuses the sampled colormaps between calls.
    """
    marker = trace_json.get('marker', {}) if isinstance(trace_json.get('marker'), dict) else {}
    trace_meta = trace_json.get('meta') if isinstance(trace_json.get('meta'), dict) else {}
    marker_color_scalars = _marker_numeric_color_values(marker, len(points))
    scalar_values = np.asarray([
        _coerce_float(point.get('color_scalar'), np.nan)
        for point in points
    ], dtype=float)
    finite_mask = np.isfinite(scalar_values)
    if not np.any(finite_mask):
        return None

    mode = str(trace_meta.get('color_by') or '').strip().lower()
    if mode in ('none', 'false', 'off'):
        return None
    if not mode:
        mode = 'age' if any(point.get('color_scalar_kind') == 'age' for point in points) else 'value'

    label = str(trace_meta.get('color_label') or '').strip()
    if not label:
        label = 'Age (Myr)' if mode == 'age' else 'Value'

    cmin = _coerce_float(marker.get('cmin'), np.nan)
    cmax = _coerce_float(marker.get('cmax'), np.nan)
    if not np.isfinite(cmin):
        cmin = float(np.nanmin(scalar_values[finite_mask]))
    if not np.isfinite(cmax):
        cmax = float(np.nanmax(scalar_values[finite_mask]))
    if not cmax > cmin:
        cmax = float(cmin + 1.0)

    selected_colormap = trace_meta.get('colormap') or marker.get('colorscale') or DEFAULT_THREEJS_TRACE_COLORMAP
    if not isinstance(selected_colormap, str):
        selected_colormap = DEFAULT_THREEJS_TRACE_COLORMAP
    try:
        colormap_options = _trace_colormap_options(selected_colormap, colormap_cache)
    except ValueError:
        colormap_options = _trace_colormap_options(DEFAULT_THREEJS_TRACE_COLORMAP, colormap_cache)

    default_color_mode = str(trace_meta.get('default_color_mode') or '').strip().lower()
    if default_color_mode not in ('fixed', 'by_value'):
        default_color_mode = 'by_value' if marker_color_scalars is not None else 'fixed'

    selected_option = colormap_options[0]
    return {
        'mode': mode,
        'label': label,
        'cmin': float(cmin),
        'cmax': float(cmax),
        'colormap': selected_option['name'],
        'colormap_options': colormap_options,
        'legend_color': selected_option.get('legend_color'),
        'default_color_mode': default_color_mode,
    }


def _trace_colormap_options(name, cache=None):
    """Colormap options for colour-by-value traces, sampled once per ``cache``; each call gets copies."""
    if cache is None:
        return _build_threejs_volume_colormap_options(name)
    if name not in cache:
        cache[name] = _build_threejs_volume_colormap_options(name)
    return [dict(option) for option in cache[name]]


# Viewer primitives: points, line segments, labels, decorations ---------------


def _galactic_simple_allowed_trace_names(plot):
    """Trace names Galactic-lite exports keep: the defaults plus every grouped trace."""
    allowed = set(GALACTIC_SIMPLE_ALLOWED_TRACE_NAMES)
    trace_grouping_dict = getattr(plot, 'trace_grouping_dict', {}) or {}
    for grouped_trace_names in trace_grouping_dict.values():
        if isinstance(grouped_trace_names, str):
            grouped_trace_names = [grouped_trace_names]
        if not isinstance(grouped_trace_names, (list, tuple, set, np.ndarray, pd.Series)):
            continue
        for trace_name in grouped_trace_names:
            if trace_name not in (None, ''):
                allowed.add(str(trace_name))
    return allowed


def _line_segments_from_trace(trace_json):
    """``[x0, y0, z0, x1, y1, z1]`` segments joining consecutive finite points (None breaks a line)."""
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
        if not (math.isfinite(point[0]) and math.isfinite(point[1]) and math.isfinite(point[2])):
            prev = None
            continue
        if prev is not None:
            segments.append(prev + point)
        prev = point
    return segments


def _points_from_trace(
    trace_json,
    default_opacity=1.0,
    include_selection=False,
    include_hovertext=True,
    include_motion=True,
    include_n_stars=True,
):
    """Point records of a markers trace: position, size, symbol, colour and opacity,
    plus hover text, selection, motion and star count as requested."""
    x_vals = _as_object_list(trace_json.get('x'))
    y_vals = _as_object_list(trace_json.get('y'))
    z_vals = _as_object_list(trace_json.get('z'))
    n_points = min(len(x_vals), len(y_vals), len(z_vals))
    if n_points == 0:
        return []

    marker = trace_json.get('marker', {})
    sizes = [_coerce_float(val, 4.0) for val in _expand_value(marker.get('size', 4.0), n_points, 4.0)]
    symbols = [str(val) for val in _expand_value(marker.get('symbol', 'circle'), n_points, 'circle')]
    hovertext = _expand_value(trace_json.get('hovertext'), n_points, '') if include_hovertext else None

    base_marker_opacity = marker.get('opacity')
    if base_marker_opacity is None:
        marker_opacity = [float(default_opacity)] * n_points
    else:
        marker_opacity = [
            _coerce_float(val, default_opacity) * float(default_opacity)
            for val in _expand_value(base_marker_opacity, n_points, default_opacity)
        ]

    colors, color_opacity = _resolve_marker_color_values(marker, n_points, default_opacity=1.0)
    color_scalars = _marker_numeric_color_values(marker, n_points)
    customdata = trace_json.get('customdata')
    custom_rows = None
    if customdata is not None and (include_selection or include_motion or include_n_stars):
        custom_arr = np.asarray(customdata, dtype=object)
        if custom_arr.ndim == 1:
            custom_arr = custom_arr.reshape(-1, 1)
        if custom_arr.ndim == 2 and custom_arr.shape[0] >= n_points:
            custom_rows = custom_arr[:n_points].tolist()

    # One vectorized clip of the products gives the same values as one clip per point.
    opacities = np.clip(
        np.asarray(marker_opacity, dtype=float) * np.asarray(color_opacity, dtype=float), 0.0, 1.0
    ).tolist()
    trace_name = trace_json.get('name')
    isfinite = math.isfinite
    points = []
    for idx in range(n_points):
        try:
            x_val = float(x_vals[idx])
            y_val = float(y_vals[idx])
            z_val = float(z_vals[idx])
        except Exception:
            continue
        if not (isfinite(x_val) and isfinite(y_val) and isfinite(z_val)):
            continue

        point = {
            'x': x_val,
            'y': y_val,
            'z': z_val,
            'size': float(max(sizes[idx], 0.0)),
            'symbol': symbols[idx],
            'color': colors[idx],
            'opacity': opacities[idx],
        }
        points.append(point)
        if color_scalars is not None and isfinite(color_scalars[idx]):
            point['color_scalar'] = float(color_scalars[idx])
        if include_hovertext:
            point['hovertext'] = str(hovertext[idx]) if hovertext[idx] is not None else ''
        if custom_rows is not None:
            row = custom_rows[idx]
            if color_scalars is None:
                age_at_t = _coerce_float(row[CUSTOMDATA_IDX_AGE_AT_T], math.nan)
                if isfinite(age_at_t):
                    point['color_scalar'] = age_at_t
                    point['color_scalar_kind'] = 'age'
            selection = _selection_from_customdata_row(row)
            if selection is not None:
                selection['trace_name'] = trace_name
                if include_motion:
                    motion = _motion_from_selection(selection)
                    if motion is not None:
                        point['motion'] = motion
                if include_n_stars:
                    n_stars = _coerce_float(selection.get('n_stars'), math.nan)
                    if isfinite(n_stars):
                        point['n_stars'] = n_stars
            if include_selection:
                point['selection'] = selection

    return points


def _annotate_threejs_point_motion_ranges(frame_specs):
    """Add each moving object's size and opacity range over all frames to its ``motion`` records."""
    isfinite = math.isfinite
    ranges = {}  # motion key -> [size_min, size_max, opacity_min, opacity_max]
    moving = []
    for frame_spec in frame_specs:
        for trace in frame_spec.get('traces', []):
            for point in trace.get('points', []):
                motion = point.get('motion')
                if not isinstance(motion, dict):
                    continue
                motion_key = str(motion.get('key') or '').strip()
                if not motion_key:
                    continue
                moving.append((motion, motion_key))
                entry = ranges.get(motion_key)
                if entry is None:
                    entry = ranges[motion_key] = [math.inf, -math.inf, math.inf, -math.inf]
                point_size = _coerce_float(point.get('size'), math.nan)
                point_opacity = _coerce_float(point.get('opacity'), math.nan)
                if isfinite(point_size):
                    entry[0] = min(entry[0], point_size)
                    entry[1] = max(entry[1], point_size)
                if isfinite(point_opacity):
                    entry[2] = min(entry[2], point_opacity)
                    entry[3] = max(entry[3], point_opacity)

    # A key without one finite size (or opacity) falls back to 0 (or 1): no
    # point of that object has a usable value of its own.
    bounds = {
        motion_key: {
            'size_min': size_min if isfinite(size_min) else 0.0,
            'size_max': size_max if isfinite(size_max) else 0.0,
            'opacity_min': opacity_min if isfinite(opacity_min) else 1.0,
            'opacity_max': opacity_max if isfinite(opacity_max) else 1.0,
        }
        for motion_key, (size_min, size_max, opacity_min, opacity_max) in ranges.items()
    }
    for motion, motion_key in moving:
        motion.update(bounds[motion_key])


def _labels_from_trace(trace_json):
    """Text labels of a ``text`` trace, with font, colour and screen-size settings."""
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
    trace_meta = trace_json.get('meta') if isinstance(trace_json.get('meta'), dict) else {}
    is_gc_label = trace_name == 'GC'
    is_galactic_radius_label = trace_name.startswith('R=') and trace_name.endswith(' Label')
    is_screen_stable_label = bool(
        is_gc_label
        or is_galactic_radius_label
        or trace_meta.get('screen_stable_text')
    )
    screen_px = None
    if is_gc_label:
        try:
            screen_px = float(trace_meta.get('screen_px', 18.0))
        except Exception:
            screen_px = 18.0
    elif is_galactic_radius_label:
        screen_px = 9.0
    elif trace_meta.get('screen_px') is not None:
        try:
            screen_px = float(trace_meta.get('screen_px'))
        except Exception:
            screen_px = None

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
            'screen_stable': is_screen_stable_label,
            'screen_px': screen_px,
        })

    return labels


def _first_trace_point(frame_json, trace_name):
    """``{'x', 'y', 'z'}`` of the first point of the frame's scatter3d trace ``trace_name``, or None."""
    for trace_json in frame_json.get('data', []):
        if trace_json.get('type', 'scatter3d') != 'scatter3d':
            continue
        if trace_json.get('name') != trace_name:
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
    return None


def _image_plane_decoration(config, default_key, plane_center, opacity_scale):
    """An ``image_plane`` frame decoration centred on ``plane_center`` in the Galactic plane."""
    return {
        'kind': 'image_plane',
        'key': str(config.get('key') or default_key),
        'center': {
            'x': float(plane_center.get('x', 0.0)),
            'y': float(plane_center.get('y', 0.0)),
            'z': 0.0,
        },
        'opacity': float(np.clip(config.get('opacity', 0.6), 0.0, 1.0)),
        'opacity_scale': opacity_scale,
        'render_order': -20,
    }


# Selections and Sky member catalogs ------------------------------------------


def _selection_from_customdata_row(row):
    """A cluster's selection record (sky position, ages, name, colour) from its customdata row."""
    if row is None:
        return None

    values = _as_object_list(row)
    count = len(values)
    if count <= CUSTOMDATA_IDX_Z0:
        return None

    l_deg = _coerce_float(values[CUSTOMDATA_IDX_L0_DEG], math.nan)
    b_deg = _coerce_float(values[CUSTOMDATA_IDX_B0_DEG], math.nan)
    if not (math.isfinite(l_deg) and math.isfinite(b_deg)):
        return None

    # _coerce_float gives finite numbers or NaN, so no further checks are needed.
    age_now = _coerce_float(values[CUSTOMDATA_IDX_AGE_NOW], math.nan)
    age_at_t = _coerce_float(values[CUSTOMDATA_IDX_AGE_AT_T], math.nan)
    click_time_myr = age_at_t - age_now
    ra_deg, dec_deg = _galactic_to_icrs_deg(l_deg, b_deg)

    cluster_name = None
    if count > CUSTOMDATA_IDX_CLUSTER_NAME and values[CUSTOMDATA_IDX_CLUSTER_NAME] not in (None, ''):
        cluster_name = str(values[CUSTOMDATA_IDX_CLUSTER_NAME])
    cluster_color = None
    if count > CUSTOMDATA_IDX_CLUSTER_COLOR and values[CUSTOMDATA_IDX_CLUSTER_COLOR] not in (None, ''):
        cluster_color = str(values[CUSTOMDATA_IDX_CLUSTER_COLOR])
    n_stars = math.nan
    if count > CUSTOMDATA_IDX_N_STARS:
        n_stars = _coerce_float(values[CUSTOMDATA_IDX_N_STARS], math.nan)
    cluster_aliases = ''
    if count > CUSTOMDATA_IDX_CLUSTER_ALIASES and values[CUSTOMDATA_IDX_CLUSTER_ALIASES] not in (None, ''):
        cluster_aliases = str(values[CUSTOMDATA_IDX_CLUSTER_ALIASES])

    return {
        'l_deg': l_deg,
        'b_deg': b_deg,
        'dist_pc': _coerce_float(values[CUSTOMDATA_IDX_DIST0_PC], math.nan),
        'x0': _coerce_float(values[CUSTOMDATA_IDX_X0], math.nan),
        'y0': _coerce_float(values[CUSTOMDATA_IDX_Y0], math.nan),
        'z0': _coerce_float(values[CUSTOMDATA_IDX_Z0], math.nan),
        'age_now_myr': age_now,
        'age_at_t_myr': age_at_t,
        'click_time_myr': click_time_myr if math.isfinite(click_time_myr) else math.nan,
        'ra_deg': ra_deg,
        'dec_deg': dec_deg,
        'cluster_name': cluster_name,
        'cluster_color': cluster_color,
        'n_stars': n_stars,
        'name_all': cluster_aliases,
    }


def _selection_identity_key(selection):
    """A stable key for a selection: cluster name, else position, else sky position, else trace."""
    if not isinstance(selection, dict):
        return ''

    cluster_name = str(selection.get('cluster_name') or '').strip()
    if cluster_name:
        return cluster_name

    x0 = _coerce_float(selection.get('x0'), math.nan)
    y0 = _coerce_float(selection.get('y0'), math.nan)
    z0 = _coerce_float(selection.get('z0'), math.nan)
    if math.isfinite(x0) and math.isfinite(y0) and math.isfinite(z0):
        return f'{x0:.6f}|{y0:.6f}|{z0:.6f}'

    ra_deg = _coerce_float(selection.get('ra_deg'), math.nan)
    dec_deg = _coerce_float(selection.get('dec_deg'), math.nan)
    if math.isfinite(ra_deg) and math.isfinite(dec_deg):
        return f'{ra_deg:.6f}|{dec_deg:.6f}'

    return str(selection.get('trace_name') or '').strip()


def _motion_from_selection(selection):
    """The record the viewer uses to fade a cluster in time, or None without both ages."""
    if not isinstance(selection, dict):
        return None

    age_now_myr = _coerce_float(selection.get('age_now_myr'), math.nan)
    age_at_t_myr = _coerce_float(selection.get('age_at_t_myr'), math.nan)
    if not (math.isfinite(age_now_myr) and math.isfinite(age_at_t_myr)):
        return None

    return {
        'key': _selection_identity_key(selection),
        'cluster_name': str(selection.get('cluster_name') or '').strip(),
        'trace_name': str(selection.get('trace_name') or '').strip(),
        'age_now_myr': float(age_now_myr),
        'age_at_t_myr': float(age_at_t_myr),
        'time_myr': float(age_at_t_myr - age_now_myr),
    }


def _normalize_threejs_cluster_catalog_key(value):
    """A cluster name reduced to lower-case letters and digits, for matching member catalogs."""
    return ''.join(character for character in str(value or '').strip().lower() if character.isalnum())


def _load_threejs_cluster_catalog(cluster_members_file, cluster_names=None):
    """Member stars per cluster from a CSV (name plus l/b or ra/dec), or None.

    With ``cluster_names`` only those clusters are read (matched on their
    normalized names). Each keeps at most ``MAX_SELECTED_MEMBER_POINTS`` stars.
    """
    if not cluster_members_file:
        return None

    try:
        header = pd.read_csv(cluster_members_file, nrows=0)
    except Exception:
        return None

    columns = list(header.columns)
    columns_by_lower = {str(column).strip().lower(): column for column in columns}

    def first_column(*candidates):
        for candidate in candidates:
            matched = columns_by_lower.get(str(candidate).strip().lower())
            if matched is not None:
                return matched
        return None

    cluster_column = first_column('name', 'cluster_name', 'cluster', 'group_name')
    l_column = first_column('l', 'l_deg', 'glon', 'galactic_l')
    b_column = first_column('b', 'b_deg', 'glat', 'galactic_b')
    ra_column = first_column('ra', 'ra_deg', 'ra_icrs', '_ra_icrs', 'raj2000')
    dec_column = first_column('dec', 'dec_deg', 'de_icrs', '_de_icrs', 'dej2000')
    source_id_column = first_column('source_id', 'gaia_source_id', 'gaiadr3', 'id')
    pmra_column = first_column('pmra', 'pm_ra_cosdec', 'pm_ra')
    pmdec_column = first_column('pmdec', 'pm_dec')
    parallax_column = first_column('parallax', 'parallax_mas')
    distance_column = first_column(
        'r_med_geo',
        'distance_pc',
        'dist_pc',
        'distance',
    )
    if cluster_column is None:
        return None
    has_galactic = l_column is not None and b_column is not None
    has_icrs = ra_column is not None and dec_column is not None
    if not has_galactic and not has_icrs:
        return None

    usecols = [cluster_column]
    for column in (
        l_column,
        b_column,
        ra_column,
        dec_column,
        source_id_column,
        pmra_column,
        pmdec_column,
        parallax_column,
        distance_column,
    ):
        if column is not None and column not in usecols:
            usecols.append(column)

    dtype = {cluster_column: str}
    if source_id_column is not None:
        dtype[source_id_column] = str
    requested_keys = {
        _normalize_threejs_cluster_catalog_key(value)
        for value in (cluster_names or [])
        if _normalize_threejs_cluster_catalog_key(value)
    }
    try:
        if requested_keys:
            matched_chunks = []
            for chunk in pd.read_csv(
                cluster_members_file,
                usecols=usecols,
                dtype=dtype,
                chunksize=250_000,
            ):
                normalized_names = chunk[cluster_column].map(
                    _normalize_threejs_cluster_catalog_key
                )
                matched = chunk.loc[normalized_names.isin(requested_keys)]
                if not matched.empty:
                    matched_chunks.append(matched)
            if not matched_chunks:
                return None
            df = pd.concat(matched_chunks, ignore_index=True)
        else:
            df = pd.read_csv(cluster_members_file, usecols=usecols, dtype=dtype)
    except Exception:
        return None

    if df.empty:
        return None
    compact_member_payload = len(df) > 100_000

    rename_columns = {cluster_column: 'cluster_name'}
    if l_column is not None:
        rename_columns[l_column] = 'l_deg'
    if b_column is not None:
        rename_columns[b_column] = 'b_deg'
    if ra_column is not None:
        rename_columns[ra_column] = 'ra_deg'
    if dec_column is not None:
        rename_columns[dec_column] = 'dec_deg'
    if source_id_column is not None:
        rename_columns[source_id_column] = 'source_id'
    if pmra_column is not None:
        rename_columns[pmra_column] = 'pmra_masyr'
    if pmdec_column is not None:
        rename_columns[pmdec_column] = 'pmdec_masyr'
    if parallax_column is not None:
        rename_columns[parallax_column] = 'parallax_mas'
    if distance_column is not None:
        rename_columns[distance_column] = 'distance_pc'
    df = df.rename(columns=rename_columns)
    df['cluster_name'] = df['cluster_name'].astype(str)
    for coordinate_column in (
        'l_deg',
        'b_deg',
        'ra_deg',
        'dec_deg',
        'pmra_masyr',
        'pmdec_masyr',
        'parallax_mas',
        'distance_pc',
    ):
        if coordinate_column in df.columns:
            df[coordinate_column] = pd.to_numeric(df[coordinate_column], errors='coerce')

    if has_galactic:
        valid_coordinates = df['l_deg'].notnull() & df['b_deg'].notnull()
    else:
        valid_coordinates = df['ra_deg'].notnull() & df['dec_deg'].notnull()
    df = df.loc[valid_coordinates].copy()
    if df.empty:
        return None

    grouped = {}
    for cluster_name, grp in df.groupby('cluster_name', sort=False):
        if has_galactic:
            l_vals = grp['l_deg'].to_numpy(dtype=np.float64)
            b_vals = grp['b_deg'].to_numpy(dtype=np.float64)
            icrs_coords = SkyCoord(
                l=l_vals * u.deg,
                b=b_vals * u.deg,
                frame='galactic',
            ).icrs
            ra_vals = np.asarray(icrs_coords.ra.deg, dtype=np.float64)
            dec_vals = np.asarray(icrs_coords.dec.deg, dtype=np.float64)
        else:
            ra_vals = grp['ra_deg'].to_numpy(dtype=np.float64)
            dec_vals = grp['dec_deg'].to_numpy(dtype=np.float64)
            galactic_coords = SkyCoord(
                ra=ra_vals * u.deg,
                dec=dec_vals * u.deg,
                frame='icrs',
            ).galactic
            l_vals = np.asarray(galactic_coords.l.deg, dtype=np.float64)
            b_vals = np.asarray(galactic_coords.b.deg, dtype=np.float64)

        idx = np.arange(l_vals.size, dtype=int)
        if idx.size > MAX_SELECTED_MEMBER_POINTS:
            # Gaia source IDs encode sky position, so catalog rows commonly
            # arrive in a spatially structured order. Evenly spaced row
            # sampling therefore imprints artificial stripes or grids on the
            # rendered sky. Select a stable pseudo-random subset from the
            # member identity and exact coordinates instead.
            sample_columns = [
                column
                for column in (
                    'source_id',
                    'l_deg',
                    'b_deg',
                    'ra_deg',
                    'dec_deg',
                )
                if column in grp.columns
            ]
            sample_hashes = pd.util.hash_pandas_object(
                grp[sample_columns],
                index=False,
                categorize=True,
            ).to_numpy(dtype=np.uint64)
            chosen = np.argpartition(
                sample_hashes,
                int(MAX_SELECTED_MEMBER_POINTS) - 1,
            )[: int(MAX_SELECTED_MEMBER_POINTS)]
            # Preserve source ordering only after the unbiased subset has been
            # chosen. This keeps exported payloads deterministic.
            idx = np.sort(chosen.astype(int))

        source_ids = (
            grp['source_id'].fillna('').astype(str).to_numpy()
            if 'source_id' in grp.columns
            else None
        )
        pmra_values = (
            grp['pmra_masyr'].to_numpy(dtype=np.float64)
            if 'pmra_masyr' in grp.columns
            else None
        )
        pmdec_values = (
            grp['pmdec_masyr'].to_numpy(dtype=np.float64)
            if 'pmdec_masyr' in grp.columns
            else None
        )
        parallax_values = (
            grp['parallax_mas'].to_numpy(dtype=np.float64)
            if 'parallax_mas' in grp.columns
            else None
        )
        distance_values = (
            grp['distance_pc'].to_numpy(dtype=np.float64)
            if 'distance_pc' in grp.columns
            else None
        )
        members = []
        for i in idx:
            member = {
                'ra': float(np.round(ra_vals[i], 6)),
                'dec': float(np.round(dec_vals[i], 6)),
                'is_cluster_member': True,
            }
            if not compact_member_payload:
                member['l'] = float(np.round(l_vals[i], 6))
                member['b'] = float(np.round(b_vals[i], 6))
            if source_ids is not None and source_ids[i]:
                member['source_id'] = str(source_ids[i])
            if pmra_values is not None and np.isfinite(pmra_values[i]):
                member['pmra_masyr'] = float(np.round(pmra_values[i], 6))
            if pmdec_values is not None and np.isfinite(pmdec_values[i]):
                member['pmdec_masyr'] = float(np.round(pmdec_values[i], 6))
            # Proper motions are angular rates. Prefer the individual star's
            # parallax distance when converting them to transverse velocities;
            # catalog distance estimates are retained only as a fallback for a
            # missing or non-positive parallax.
            distance_pc = np.nan
            if parallax_values is not None:
                parallax_mas = float(parallax_values[i])
                if np.isfinite(parallax_mas) and parallax_mas > 0.0:
                    distance_pc = 1000.0 / parallax_mas
            if not np.isfinite(distance_pc):
                distance_pc = (
                    float(distance_values[i])
                    if distance_values is not None
                    and np.isfinite(distance_values[i])
                    and distance_values[i] > 0.0
                    else np.nan
                )
            if np.isfinite(distance_pc) and distance_pc > 0.0:
                member['distance_pc'] = float(np.round(distance_pc, 3))
            members.append(member)
        grouped[str(cluster_name).strip()] = members

    return grouped or None


def _catalog_point_from_selection(selection, default_label=''):
    """A Sky catalog point (l, b, ra, dec, label) from a selection, or None."""
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
    """Copies of at most ``MAX_SELECTED_MEMBER_POINTS`` evenly spaced points."""
    if not points:
        return []
    if len(points) <= MAX_SELECTED_MEMBER_POINTS:
        return copy.deepcopy(points)

    idx = np.linspace(0, len(points) - 1, int(MAX_SELECTED_MEMBER_POINTS), dtype=int)
    return [copy.deepcopy(points[int(i)]) for i in idx]


def _threejs_catalog_from_frame_spec(frame_spec):
    """Sky catalog points per trace, cluster name and alias from a frame's selections."""
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
            cluster_aliases = selection.get('name_all') if isinstance(selection, dict) else None
            if cluster_aliases:
                for alias in str(cluster_aliases).replace(';', ',').replace('|', ',').split(','):
                    alias = alias.strip()
                    if alias:
                        grouped_points.setdefault(alias, []).append(catalog_point)

        if trace_name and trace_points:
            catalogs[str(trace_name)] = _limit_catalog_points(trace_points)
        for cluster_name, points in grouped_points.items():
            catalogs.setdefault(cluster_name, _limit_catalog_points(points))

    return catalogs


def _merge_threejs_member_catalogs(*catalogs):
    """Merge catalogs (later ones win per key), copying the point lists."""
    merged = {}
    for catalog in catalogs:
        if not isinstance(catalog, dict):
            continue
        for key, points in catalog.items():
            if not key or not isinstance(points, list) or not points:
                continue
            merged[str(key)] = copy.deepcopy(points)
    return merged


# Volumes: configuration ------------------------------------------------------


def _normalize_threejs_volume_configs(volumes):
    """``make_plot(volumes=...)`` as a list of config dicts (a path becomes ``{'path': ...}``)."""
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


def _normalize_threejs_volume_lighting_mode(mode):
    value = str(mode or 'standard').strip().lower().replace('-', '_')
    if value in {'galactic', 'galactic_center', 'galactic_centre'}:
        return 'galactic'
    return 'standard'


def _coerce_threejs_volume_galactic_center(value):
    if isinstance(value, dict):
        source = [value.get('x'), value.get('y'), value.get('z')]
    elif isinstance(value, (list, tuple, np.ndarray)) and len(value) >= 3:
        source = value[:3]
    else:
        source = [8122.0, 0.0, 0.0]
    return [
        float(_coerce_float(source[0], 8122.0)),
        float(_coerce_float(source[1], 0.0)),
        float(_coerce_float(source[2], 0.0)),
    ]


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


def _threejs_volume_variant_metadata(volume_cfg):
    metadata = {}
    for source_key in ('variant_group', 'variant_label', 'base_state_name'):
        value = str(volume_cfg.get(source_key) or '').strip()
        if value:
            metadata[source_key] = value
    if volume_cfg.get('variant_order') is not None:
        try:
            variant_order = float(volume_cfg.get('variant_order'))
        except (TypeError, ValueError):
            variant_order = None
        if variant_order is not None and np.isfinite(variant_order):
            metadata['variant_order'] = int(variant_order) if variant_order.is_integer() else float(variant_order)
    return metadata


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


def _normalize_threejs_volume_clip_bounds(volume_cfg):
    clip_cfg = (
        volume_cfg.get('clip_bounds')
        or volume_cfg.get('crop_bounds')
        or volume_cfg.get('data_bounds')
    )
    if not isinstance(clip_cfg, dict):
        return {}

    normalized = {}
    for axis_name in ('x', 'y', 'z'):
        raw_bounds = clip_cfg.get(axis_name)
        if not isinstance(raw_bounds, (list, tuple, np.ndarray)) or len(raw_bounds) != 2:
            continue
        try:
            lower = float(raw_bounds[0])
            upper = float(raw_bounds[1])
        except Exception:
            continue
        if not np.isfinite(lower) or not np.isfinite(upper) or np.isclose(lower, upper):
            continue
        normalized[axis_name] = [float(min(lower, upper)), float(max(lower, upper))]
    return normalized


def _threejs_volume_max_resolution_cap(volume_cfg, data_shape_zyx, *, minimum):
    """The largest resampled size allowed per axis: the config's cap within [minimum, data size]."""
    data_limit = max(int(v) for v in data_shape_zyx)
    requested_cap = volume_cfg.get(
        'max_resolution_cap',
        volume_cfg.get('max_resolution_limit', DEFAULT_THREEJS_VOLUME_MAX_RESOLUTION_CAP),
    )
    return int(
        np.clip(
            _coerce_float(requested_cap, DEFAULT_THREEJS_VOLUME_MAX_RESOLUTION_CAP),
            int(minimum),
            data_limit,
        )
    )


# Volumes: layer specs --------------------------------------------------------


def _xyz_shape(shape_zyx):
    """``{'x', 'y', 'z'}`` sizes of a (z, y, x) array shape."""
    return {'x': int(shape_zyx[2]), 'y': int(shape_zyx[1]), 'z': int(shape_zyx[0])}


def _xyz_downsample_step(source_shape_zyx, sampled_shape_zyx):
    """Source voxels per sampled voxel along x, y and z."""
    return {
        axis: float(source_shape_zyx[index]) / float(sampled_shape_zyx[index])
        for axis, index in (('x', 2), ('y', 1), ('z', 0))
    }


def _quantize_volume_uint8(data, data_min, data_max):
    """Map ``data`` linearly from [data_min, data_max] to 0-255; non-finite voxels become 0."""
    finite_mask = np.isfinite(data)
    normalized = np.zeros_like(data, dtype=np.float32)
    if data_max > data_min:
        normalized[finite_mask] = (data[finite_mask] - data_min) / (data_max - data_min)
        normalized = np.clip(normalized, 0.0, 1.0)
    else:
        normalized[finite_mask] = 1.0
    quantized = np.zeros_like(data, dtype=np.uint8)
    quantized[finite_mask] = np.rint(normalized[finite_mask] * 255.0).astype(np.uint8)
    return np.ascontiguousarray(quantized)


def _uint8_b64(quantized):
    """Raw bytes of a uint8 array as base64 text."""
    return base64.b64encode(quantized.tobytes(order='C')).decode('ascii')


# iOS Safari fails (silently, via canvas limits) on very large decode surfaces
# and spikes hundreds of MB of RGBA on getImageData, so atlases stay small.
PNG_ATLAS_MAX_SIDE_PX = 2048


def _encode_volume_png_atlas(quantized):
    """Encode a uint8 cube as PNG slice atlases.

    Returns ``(data_b64, encoding, tiles, slabs)``: one atlas when it fits
    within ``PNG_ATLAS_MAX_SIDE_PX``, otherwise a list of slab atlases.
    Without Pillow it falls back to raw base64 bytes (``'uint8'``).
    """
    try:
        from PIL import Image
    except ImportError:
        return _uint8_b64(quantized), 'uint8', None, None

    def _png_b64(array_2d):
        buffer = io.BytesIO()
        Image.fromarray(array_2d, mode='L').save(
            buffer,
            format='PNG',
            optimize=True,
            compress_level=9,
        )
        return base64.b64encode(buffer.getvalue()).decode('ascii')

    nz, ny, nx = (int(quantized.shape[0]), int(quantized.shape[1]), int(quantized.shape[2]))
    tile_cols = int(np.ceil(np.sqrt(max(nz, 1))))
    tile_rows = int(np.ceil(nz / tile_cols))
    max_atlas_side_px = PNG_ATLAS_MAX_SIDE_PX
    if tile_cols * nx <= max_atlas_side_px and tile_rows * ny <= max_atlas_side_px:
        atlas = np.zeros((tile_rows * ny, tile_cols * nx), dtype=np.uint8)
        for z_index in range(nz):
            row = z_index // tile_cols
            col = z_index % tile_cols
            atlas[row * ny:(row + 1) * ny, col * nx:(col + 1) * nx] = quantized[z_index]
        return (
            _png_b64(atlas),
            'png_atlas_uint8',
            {'x': tile_cols, 'y': tile_rows},
            None,
        )

    slab_cols = max(1, min(nz, max_atlas_side_px // max(nx, 1)))
    slab_rows = max(1, min(int(np.ceil(nz / slab_cols)), max_atlas_side_px // max(ny, 1)))
    slices_per_slab = max(1, slab_cols * slab_rows)
    slabs = []
    for start in range(0, nz, slices_per_slab):
        chunk = quantized[start:start + slices_per_slab]
        rows_needed = int(np.ceil(chunk.shape[0] / slab_cols))
        atlas = np.zeros((rows_needed * ny, slab_cols * nx), dtype=np.uint8)
        for offset in range(chunk.shape[0]):
            row = offset // slab_cols
            col = offset % slab_cols
            atlas[row * ny:(row + 1) * ny, col * nx:(col + 1) * nx] = chunk[offset]
        slabs.append(_png_b64(atlas))
    return (
        '',
        'png_atlas_uint8',
        {'x': slab_cols, 'y': slab_rows, 'slices_per_slab': slices_per_slab},
        slabs,
    )


def _build_threejs_volume_ar_proxy(sampled_data, data_min, data_max, volume_cfg):
    """A block-max downsampled uint8 copy of the volume for AR exports, or None."""
    if volume_cfg.get('ar_proxy_enabled') is False:
        return None

    source = np.asarray(sampled_data, dtype=np.float32)
    if source.ndim != 3 or source.size == 0:
        return None

    max_resolution = int(
        np.clip(
            _coerce_float(
                volume_cfg.get('ar_proxy_max_resolution'),
                DEFAULT_THREEJS_AR_VOLUME_MAX_RESOLUTION,
            ),
            8,
            96,
        )
    )
    source_shape = tuple(int(value) for value in source.shape)
    steps = tuple(max(1, int(np.ceil(value / max_resolution))) for value in source_shape)
    pooled = source
    for axis, step in enumerate(steps):
        if step <= 1:
            continue
        indices = np.arange(0, pooled.shape[axis], step, dtype=np.intp)
        pooled = np.fmax.reduceat(pooled, indices, axis=axis)

    quantized = _quantize_volume_uint8(pooled, data_min, data_max)
    return {
        'data_b64': _uint8_b64(quantized),
        'data_encoding': 'uint8',
        'shape': _xyz_shape(quantized.shape),
        'downsample_step': _xyz_downsample_step(source_shape, quantized.shape),
        'method': 'block_max',
    }


def _threejs_volume_stats_values(sampled):
    """Finite positive samples (all finite ones if none is positive), for value ranges."""
    finite_mask = np.isfinite(sampled)
    positive_mask = finite_mask & (sampled > 0)
    return sampled[positive_mask] if np.any(positive_mask) else sampled[finite_mask]


def _threejs_volume_default_range(volume_cfg, stats_values, data_min, data_max):
    """The (vmin, vmax) a volume opens with: ``vmin``/``vmax`` or quantiles of ``stats_values``."""
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
    return default_vmin, default_vmax


def _threejs_volume_default_controls(volume_cfg, vmin, vmax, colormap_name):
    """The rendering controls a volume opens with, clipped to the viewer's ranges."""
    def clipped(key, default, low, high):
        return float(np.clip(_coerce_float(volume_cfg.get(key), default), low, high))

    sample_setting = volume_cfg.get('samples', volume_cfg.get('steps', DEFAULT_THREEJS_VOLUME_SAMPLE_STEPS))
    steps = int(np.clip(_coerce_float(sample_setting, DEFAULT_THREEJS_VOLUME_SAMPLE_STEPS), 24, 768))
    return {
        'vmin': float(vmin),
        'vmax': float(vmax),
        'opacity': clipped('opacity', 0.15, 0.0, 1.0),
        'steps': steps,
        'samples': steps,
        'alpha_coef': clipped('alpha_coef', 50.0, 1.0, 200.0),
        'gradient_step': clipped('gradient_step', 0.005, 1e-4, 0.05),
        'stretch': str(_normalize_threejs_volume_stretch(volume_cfg.get('stretch', 'linear'))),
        'colormap': colormap_name,
        'show_all_times': bool(volume_cfg.get('show_all_times', False)),
        'lighting_mode': str(_normalize_threejs_volume_lighting_mode(volume_cfg.get('lighting_mode'))),
        'galactic_center': _coerce_threejs_volume_galactic_center(volume_cfg.get('galactic_center')),
        'galactic_light_intensity': clipped('galactic_light_intensity', 1.35, 0.0, 4.0),
        'galactic_ambient': clipped('galactic_ambient', 0.22, 0.0, 1.0),
        'galactic_extinction': clipped('galactic_extinction', 2.4, 0.0, 8.0),
        'galactic_scattering': clipped('galactic_scattering', 0.55, 0.0, 2.0),
        'galactic_anisotropy': clipped('galactic_anisotropy', 0.45, 0.0, 0.9),
        'galactic_warmth': clipped('galactic_warmth', 0.72, 0.0, 1.0),
    }


def _threejs_volume_bounds(axis_bounds, bound_offset, center_offset, apply_center_offset):
    """Display bounds per axis: the data bounds plus ``bound_offset``, recentred when asked."""
    bounds = {}
    for axis, (lower, upper) in zip(('x', 'y', 'z'), axis_bounds):
        bounds[axis] = [
            float((value + bound_offset[axis]) - center_offset[axis])
            if apply_center_offset else float(value + bound_offset[axis])
            for value in (lower, upper)
        ]
    return bounds


def _threejs_volume_identity(volume_cfg, key, name, time_myr):
    """The leading keys of a volume layer spec: key, State key/name, variants, time, name."""
    return {
        'key': key,
        'state_key': _threejs_volume_state_key(volume_cfg, key),
        'state_name': str(volume_cfg.get('state_name') or volume_cfg.get('legend_name') or name),
        **_threejs_volume_variant_metadata(volume_cfg),
        'time_myr': time_myr,
        'name': name,
    }


def _threejs_volume_display(volume_cfg, time_myr, data_min, data_max, value_unit, colormap_options):
    """Value range, legend colour, visibility and timing keys of a volume layer spec."""
    return {
        'data_range': [float(data_min), float(data_max)],
        'value_unit': value_unit,
        'legend_color': colormap_options[0].get('legend_color'),
        'visible': bool(volume_cfg.get('visible', True)),
        'only_at_t0': bool(volume_cfg.get('only_at_t0', time_myr is None)),
        'supports_show_all_times': bool(volume_cfg.get('supports_show_all_times', time_myr is None)),
        'co_rotate_with_frame': bool(volume_cfg.get('co_rotate_with_frame', False)),
        'reference_time_myr': float(_coerce_float(volume_cfg.get('reference_time_myr'), 0.0)),
        'interpolation': bool(volume_cfg.get('interpolation', True)),
    }


def _build_threejs_inline_volume_layer_spec(volume_cfg, center_offset=None, index=0):
    """Volume layer spec for an in-memory (z, y, x) array given as ``volume_cfg['data']``."""
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
        max_resolution_cap = _threejs_volume_max_resolution_cap(volume_cfg, data_shape_zyx, minimum=8)
        max_resolution = int(np.clip(_coerce_float(max_resolution, max(data_shape_zyx)), 8, max_resolution_cap))
        target_shape_zyx = tuple(min(max_resolution, dim) for dim in data_shape_zyx)
    sampled = (
        _downsample_threejs_volume_with_zoom(data, target_shape_zyx)
        if target_shape_zyx != data_shape_zyx
        else np.ascontiguousarray(data)
    )

    stats_values = _threejs_volume_stats_values(sampled)
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
    default_vmin, default_vmax = _threejs_volume_default_range(volume_cfg, stats_values, data_min, data_max)

    bounds_cfg = volume_cfg.get('bounds') or {}
    axis_bounds = (
        _coerce_range(bounds_cfg.get('x'), [-0.5 * data_shape_zyx[2], 0.5 * data_shape_zyx[2]]),
        _coerce_range(bounds_cfg.get('y'), [-0.5 * data_shape_zyx[1], 0.5 * data_shape_zyx[1]]),
        _coerce_range(bounds_cfg.get('z'), [-0.5 * data_shape_zyx[0], 0.5 * data_shape_zyx[0]]),
    )
    apply_center_offset = bool(volume_cfg.get('apply_center_offset', False))
    bound_offset = _coerce_threejs_volume_bound_offset(volume_cfg.get('bound_offset'))

    name = str(volume_cfg.get('name') or f'Volume {index + 1}')
    key = str(volume_cfg.get('key') or f'volume-{index}')
    time_myr = _coerce_threejs_volume_time_myr(volume_cfg.get('time_myr'))
    opacity_function = _normalize_threejs_volume_opacity_function(volume_cfg.get('opacity_function'))
    colormap_options = _build_threejs_volume_colormap_options(
        volume_cfg.get('colormap', 'inferno'),
        opacity_function=opacity_function,
    )
    value_unit = str(volume_cfg.get('unit_label') or volume_cfg.get('value_unit') or '').strip()

    return {
        **_threejs_volume_identity(volume_cfg, key, name, time_myr),
        'path': '<inline>',
        'hdu': 'inline',
        'data_b64': _uint8_b64(_quantize_volume_uint8(sampled, data_min, data_max)),
        'data_encoding': 'uint8',
        'ar_proxy': _build_threejs_volume_ar_proxy(sampled, data_min, data_max, volume_cfg),
        'shape': _xyz_shape(sampled.shape),
        'source_shape': _xyz_shape(data_shape_zyx),
        'downsample_step': _xyz_downsample_step(data_shape_zyx, sampled.shape),
        'downsample_method': 'scipy_zoom' if sampled.shape != data.shape else 'inline',
        'bounds': _threejs_volume_bounds(axis_bounds, bound_offset, center_offset, apply_center_offset),
        **_threejs_volume_display(volume_cfg, time_myr, data_min, data_max, value_unit, colormap_options),
        'sky_overlay_data_b64': None,
        'sky_overlay_data_encoding': None,
        'sky_overlay_shape': None,
        'sky_overlay_atlas_tiles': None,
        'sky_overlay_downsample_step': None,
        'default_controls': _threejs_volume_default_controls(
            volume_cfg, default_vmin, default_vmax, colormap_options[0]['name']
        ),
        'colormap_options': colormap_options,
    }


def _build_threejs_volume_layer_spec(volume_cfg, center_offset=None, index=0, include_sky_overlay=False):
    """Volume layer spec for a FITS cube (``path``), or an inline array without one.

    The cube is clipped to ``clip_bounds``, resampled to at most
    ``max_resolution`` voxels per axis and quantized to uint8; with
    ``include_sky_overlay`` a finer copy is added for the Sky view.
    """
    from astropy.io import fits

    center_offset = center_offset or {'x': 0.0, 'y': 0.0, 'z': 0.0}
    path = volume_cfg.get('path') or volume_cfg.get('fits_file')
    if not path:
        return _build_threejs_inline_volume_layer_spec(volume_cfg, center_offset=center_offset, index=index)

    hdu_selector = volume_cfg.get('hdu')
    requested_max_resolution = volume_cfg.get('max_resolution')
    requested_sky_overlay_max_resolution = volume_cfg.get('sky_overlay_max_resolution')
    apply_center_offset = bool(volume_cfg.get('apply_center_offset', True))
    bound_offset = _coerce_threejs_volume_bound_offset(volume_cfg.get('bound_offset'))

    with fits.open(path, memmap=True) as hdul:
        hdu_info = _resolve_threejs_volume_hdu(hdul, hdu_selector)
        hdu = hdu_info['hdu']
        data = hdu_info['cube_data_zyx']
        data_shape_zyx = tuple(int(v) for v in data.shape)
        original_data_shape_zyx = data_shape_zyx
        header = hdu.header.copy()
        axis_numbers_zyx = hdu_info['axis_numbers_zyx']
        resolved_hdu_label = hdu_info['resolved_hdu']

        # Axis index in (z, y, x) order for each display axis.
        axis_indices = {'x': 2, 'y': 1, 'z': 0}
        axis_bounds = tuple(
            _threejs_volume_axis_bounds(
                header,
                axis_number=int(axis_numbers_zyx[axis_indices[axis]]),
                axis_size=int(data_shape_zyx[axis_indices[axis]]),
            )
            for axis in ('x', 'y', 'z')
        )

        clip_bounds = _normalize_threejs_volume_clip_bounds(volume_cfg)
        if clip_bounds:
            slices = {}
            clipped_bounds = []
            for axis in ('x', 'y', 'z'):
                slices[axis], axis_bound = _threejs_volume_clip_axis_slice(
                    header,
                    axis_number=int(axis_numbers_zyx[axis_indices[axis]]),
                    axis_size=int(data_shape_zyx[axis_indices[axis]]),
                    clip_bounds=clip_bounds.get(axis),
                    center_offset=float(center_offset.get(axis, 0.0)),
                    bound_offset=float(bound_offset[axis]),
                    apply_center_offset=apply_center_offset,
                    axis_name=axis,
                    path=path,
                )
                clipped_bounds.append(axis_bound)
            axis_bounds = tuple(clipped_bounds)
            data = data[slices['z'], slices['y'], slices['x']]
            data_shape_zyx = tuple(int(v) for v in data.shape)

        max_resolution_cap = _threejs_volume_max_resolution_cap(volume_cfg, data_shape_zyx, minimum=24)
        max_resolution = int(np.clip(_coerce_float(requested_max_resolution, 96), 24, max_resolution_cap))
        target_shape_zyx = tuple(min(max_resolution, dim) for dim in data_shape_zyx)
        sampled = _downsample_threejs_volume_with_zoom(data, target_shape_zyx)
        sky_overlay_sampled = None
        if include_sky_overlay:
            sky_overlay_max_resolution = int(
                np.clip(
                    _coerce_float(requested_sky_overlay_max_resolution, max(max_resolution, 256)),
                    64,
                    max_resolution_cap,
                )
            )
            sky_overlay_shape_zyx = tuple(min(sky_overlay_max_resolution, dim) for dim in data_shape_zyx)
            if sky_overlay_shape_zyx == target_shape_zyx:
                sky_overlay_sampled = sampled
            else:
                sky_overlay_sampled = _downsample_threejs_volume_with_zoom(data, sky_overlay_shape_zyx)

    stats_values = _threejs_volume_stats_values(
        sky_overlay_sampled if sky_overlay_sampled is not None else sampled
    )
    if stats_values.size == 0:
        raise ValueError(f"Volume FITS cube contains no finite samples: {path}")

    data_min = float(np.nanmin(stats_values))
    data_max = float(np.nanmax(stats_values))
    if not np.isfinite(data_min) or not np.isfinite(data_max):
        raise ValueError(f"Volume FITS cube statistics are not finite: {path}")
    default_vmin, default_vmax = _threejs_volume_default_range(volume_cfg, stats_values, data_min, data_max)

    data_encoding_request = str(volume_cfg.get('data_encoding') or '').strip().lower()
    use_png_atlas = bool(
        data_encoding_request == 'png_atlas_uint8'
        or volume_cfg.get('encode_as_png_atlas')
        or volume_cfg.get('use_png_atlas')
    )
    quantized = _quantize_volume_uint8(sampled, data_min, data_max)
    data_atlas_tiles = None
    data_b64_slabs = None
    if use_png_atlas:
        data_b64, data_encoding, data_atlas_tiles, data_b64_slabs = _encode_volume_png_atlas(quantized)
    else:
        data_b64 = _uint8_b64(quantized)
        data_encoding = 'uint8'
    sky_overlay_b64 = None
    sky_overlay_encoding = None
    sky_overlay_tiles = None
    if sky_overlay_sampled is not None:
        sky_overlay_b64, sky_overlay_encoding, sky_overlay_tiles, _sky_overlay_slabs = _encode_volume_png_atlas(
            _quantize_volume_uint8(sky_overlay_sampled, data_min, data_max)
        )

    opacity_function = _normalize_threejs_volume_opacity_function(volume_cfg.get('opacity_function'))
    colormap_options = _build_threejs_volume_colormap_options(
        volume_cfg.get('colormap', 'inferno'),
        opacity_function=opacity_function,
    )
    name = str(volume_cfg.get('name') or str(path).rsplit('/', 1)[-1].rsplit('.', 1)[0] or f'Volume {index + 1}')
    value_unit = str(volume_cfg.get('unit_label') or header.get('BUNIT') or '').strip()
    key = str(volume_cfg.get('key') or f'volume-{index}')
    time_myr = _coerce_threejs_volume_time_myr(volume_cfg.get('time_myr'))

    return {
        **_threejs_volume_identity(volume_cfg, key, name, time_myr),
        'path': str(path),
        'hdu': str(resolved_hdu_label),
        'data_b64': data_b64,
        'data_b64_slabs': data_b64_slabs,
        'data_encoding': data_encoding,
        'data_atlas_tiles': data_atlas_tiles,
        'ar_proxy': _build_threejs_volume_ar_proxy(sampled, data_min, data_max, volume_cfg),
        'shape': _xyz_shape(sampled.shape),
        'source_shape': _xyz_shape(data_shape_zyx),
        'original_source_shape': _xyz_shape(original_data_shape_zyx),
        'downsample_step': _xyz_downsample_step(data_shape_zyx, sampled.shape),
        'downsample_method': 'scipy_zoom',
        'bounds': _threejs_volume_bounds(axis_bounds, bound_offset, center_offset, apply_center_offset),
        **_threejs_volume_display(volume_cfg, time_myr, data_min, data_max, value_unit, colormap_options),
        'sky_overlay_data_b64': sky_overlay_b64,
        'sky_overlay_data_encoding': sky_overlay_encoding,
        'sky_overlay_shape': (
            _xyz_shape(sky_overlay_sampled.shape) if sky_overlay_sampled is not None else None
        ),
        'sky_overlay_atlas_tiles': sky_overlay_tiles,
        'sky_overlay_downsample_step': (
            _xyz_downsample_step(data_shape_zyx, sky_overlay_sampled.shape)
            if sky_overlay_sampled is not None
            else None
        ),
        'default_controls': _threejs_volume_default_controls(
            volume_cfg, default_vmin, default_vmax, colormap_options[0]['name']
        ),
        'colormap_options': colormap_options,
    }


# Volumes: FITS cubes ---------------------------------------------------------


def _resolve_threejs_volume_hdu(hdul, hdu_selector):
    """The FITS HDU holding the cube: the requested one, or the best 3D candidate."""
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
    """A numeric 3D (z, y, x) cube without singleton axes, and its FITS axis numbers."""
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
    """``data`` resampled (linear) to ``target_shape_zyx`` as float32."""
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
    """Outer edges of a FITS axis from its WCS keywords, or +/- half its size without them."""
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


def _threejs_volume_axis_centers(header, axis_number, axis_size):
    """Voxel-centre coordinates along a FITS axis (indices without WCS keywords)."""
    cdelt_key = f'CDELT{axis_number}'
    crval_key = f'CRVAL{axis_number}'
    crpix_key = f'CRPIX{axis_number}'
    if (cdelt_key not in header) or (crval_key not in header) or (crpix_key not in header):
        return np.arange(int(axis_size), dtype=float)

    delta = float(header.get(cdelt_key, 1.0))
    crval = float(header.get(crval_key, 0.0))
    crpix = float(header.get(crpix_key, 1.0))
    if (not np.isfinite(delta)) or (abs(delta) <= 0.0) or (not np.isfinite(crval)) or (not np.isfinite(crpix)):
        return np.arange(int(axis_size), dtype=float)
    pixel_numbers = np.arange(1, int(axis_size) + 1, dtype=float)
    return crval + ((pixel_numbers - crpix) * delta)


def _threejs_volume_clip_axis_slice(
    header,
    *,
    axis_number,
    axis_size,
    clip_bounds,
    center_offset,
    bound_offset,
    apply_center_offset,
    axis_name,
    path,
):
    full_bounds = _threejs_volume_axis_bounds(header, axis_number=axis_number, axis_size=axis_size)
    if not clip_bounds:
        return slice(None), full_bounds

    display_lower, display_upper = [float(v) for v in clip_bounds]
    raw_lower = display_lower - float(bound_offset)
    raw_upper = display_upper - float(bound_offset)
    if apply_center_offset:
        raw_lower += float(center_offset)
        raw_upper += float(center_offset)
    raw_lower, raw_upper = min(raw_lower, raw_upper), max(raw_lower, raw_upper)

    centers = _threejs_volume_axis_centers(header, axis_number=axis_number, axis_size=axis_size)
    mask = (centers >= raw_lower) & (centers <= raw_upper)
    if not np.any(mask):
        raise ValueError(
            f"Volume clip_bounds for axis {axis_name!r} do not overlap FITS cube {path}: "
            f"{clip_bounds!r}"
        )

    indices = np.flatnonzero(mask)
    start = int(indices[0])
    stop = int(indices[-1]) + 1
    return slice(start, stop), (float(raw_lower), float(raw_upper))


# Volumes: colormaps ----------------------------------------------------------


def _normalize_threejs_volume_opacity_function(opacity_function):
    """An Nx2 array of sorted (position, alpha) points spanning [0, 1]; linear by default."""
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
    """The selected colormap followed by the defaults (and their reverses) as 1024-step RGBA LUTs."""
    options = []
    seen = set()
    requested = []

    def _append_requested(value):
        if value in (None, ''):
            return
        requested.append(value)

    _append_requested(selected_colormap)
    for base_name in DEFAULT_THREEJS_VOLUME_COLORMAPS:
        _append_requested(base_name)
        _append_requested(_threejs_reversed_colormap_name(base_name))

    for candidate in requested:
        sampled = _sample_threejs_volume_colormap(candidate, opacity_function=opacity_function)
        if sampled is None:
            continue
        name, rgba = sampled
        key = _normalize_threejs_volume_colormap_name(name)
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
    """``(name, RGBA uint8 LUT)`` of a colormap with the opacity curve applied, or None."""
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
    """``(label, callable)`` for a colormap name (matplotlib or ``colormaps``) or callable, or None."""
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

    normalized_requested = _normalize_threejs_volume_colormap_name(requested)
    if not normalized_requested:
        return None

    is_reversed = normalized_requested.endswith('_r')
    base_normalized = normalized_requested[:-2] if is_reversed else normalized_requested
    candidates = _threejs_volume_colormap_name_candidates(requested)
    if is_reversed:
        candidates.extend(_threejs_volume_colormap_name_candidates(base_normalized))
    ordered_candidates = []
    for name in candidates:
        if name not in ordered_candidates:
            ordered_candidates.append(name)

    try:
        from matplotlib import colormaps as mpl_colormaps

        for name in ordered_candidates:
            try:
                cmap = mpl_colormaps[name]
                return _finalize_threejs_volume_colormap(name, cmap, reversed_requested=is_reversed)
            except Exception:
                continue
    except Exception:
        pass

    try:
        import colormaps as installed_colormaps

        for name in ordered_candidates:
            cmap = getattr(installed_colormaps, name, None)
            if callable(cmap):
                return _finalize_threejs_volume_colormap(name, cmap, reversed_requested=is_reversed)
    except Exception:
        pass

    return None


def _normalize_threejs_volume_colormap_name(name):
    if name in (None, ''):
        return ''
    normalized = str(name).strip()
    if not normalized:
        return ''
    return normalized.replace('-', '_').replace(' ', '_').lower()


def _threejs_reversed_colormap_name(name):
    normalized = _normalize_threejs_volume_colormap_name(name)
    if not normalized:
        return ''
    if normalized.endswith('_r'):
        return normalized
    return f'{normalized}_r'


def _threejs_volume_colormap_name_candidates(name):
    normalized = _normalize_threejs_volume_colormap_name(name)
    if not normalized:
        return []

    candidates = [normalized]
    plain_name = normalized[:-2] if normalized.endswith('_r') else normalized
    if plain_name:
        title_name = plain_name.title()
        candidates.extend((plain_name, title_name))
        if normalized.endswith('_r'):
            candidates.append(f'{title_name}_r')

    ordered = []
    for candidate in candidates:
        if candidate and candidate not in ordered:
            ordered.append(candidate)
    return ordered


def _finalize_threejs_volume_colormap(name, cmap, reversed_requested=False):
    normalized_name = _normalize_threejs_volume_colormap_name(name)
    label = getattr(cmap, 'name', None) or normalized_name or 'custom'

    if not reversed_requested:
        return str(label), cmap

    if normalized_name.endswith('_r'):
        return str(label), cmap

    if hasattr(cmap, 'reversed'):
        try:
            reversed_cmap = cmap.reversed()
            reversed_label = getattr(reversed_cmap, 'name', None) or _threejs_reversed_colormap_name(label)
            return str(reversed_label), reversed_cmap
        except Exception:
            pass

    def reversed_callable(values, _base_cmap=cmap):
        return _base_cmap(np.clip(1.0 - np.asarray(values, dtype=float), 0.0, 1.0))

    reversed_callable.name = _threejs_reversed_colormap_name(label)
    return str(reversed_callable.name), reversed_callable
