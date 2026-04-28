# Classes for storing user input data
import numpy as np
import pandas as pd
from . import orbit_maker
from . import point_sizes
import copy

class Trace:
    """
    Class for storing user input data of a star cluster.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data containing the following columns: x, y, z, U, V, W, name, age_myr.
    data_name : str
        Name of the star cluster data.
    min_size : int, optional
        Minimum size of the marker points (default is 1).
    max_size : int, optional
        Maximum size of the marker points (default is 5).
    color : str, optional
        Color of the marker points (default is 'gray').
    colormap : str, optional
        Plotly colorscale name for age-based coloring (default is None).
    cmin : float, optional
        Minimum value for the colorscale (default is None, auto-computed).
    cmax : float, optional
        Maximum value for the colorscale (default is None, auto-computed).
    opacity : float, optional
        Opacity of the marker points (default is 1.0).
    marker_style : str, optional
        Style of the marker points (default is 'circle').
    show_tracks : bool, optional
        Whether to show tracks (default is False).
    size_by_n_stars : bool, optional
        Whether to size by number of stars (default is False).
    fade_in_time : float, optional
        Duration of time for size easing (default is None, uses make_plot default).
    fade_in_and_out : bool, optional
        Whether to fade in and out (default is None, uses make_plot default).
    fade_in_and_disp : bool, optional
        Whether to fade in then drop to min size after disp_time (default is None, uses make_plot default).
    disp_time : float, optional
        Time after birth to keep max size before dropping to min size (default is None, uses make_plot default).

    Attributes
    ----------
    df : pandas.DataFrame
        Input data containing the following columns: x, y, z, U, V, W, name, age_myr.
    data_name : str
        Name of the star cluster data.
    cluster_int_coords : tuple
        Integrated coordinates of the star cluster.
    coordinates : numpy.ndarray
        Array of coordinates (x, y, z, U, V, W) of the star cluster.
    integrated : bool
        Whether the star cluster has been integrated.
    color : str
        Color of the marker points.
    colormap : str or None
        Plotly colorscale name for age-based coloring.
    cmin : float or None
        Minimum value for the colorscale.
    cmax : float or None
        Maximum value for the colorscale.
    opacity : float
        Opacity of the marker points.
    min_size : int
        Minimum size of the marker points.
    max_size : int
        Maximum size of the marker points.
    marker_style : str
        Style of the marker points.
    show_tracks : bool
        Whether to show tracks.
    size_by_n_stars : bool
        Whether to size by number of stars.
    fade_in_time : float or None
        Duration of time for size easing (instance-specific).
    fade_in_and_out : bool or None
        Whether to fade in and out (instance-specific).
    fade_in_and_disp : bool or None
        Whether to fade in then drop to min size after disp_time (instance-specific).
    disp_time : float or None
        Time after birth to keep max size before dropping to min size (instance-specific).

    Methods
    -------
    __init__(self, df, data_name, min_size=1, max_size=5, color='gray', opacity=1.0, marker_style='circle', show_tracks=False, size_by_n_stars=False, fade_in_time=None, fade_in_and_out=None, colormap=None, cmin=None, cmax=None, fade_in_and_disp=None, disp_time=None)
        Initializes the Trace object.
    create_integrated_dataframe(self, time)
        Creates an integrated DataFrame of the star cluster.
    set_age_based_sizes(self, fade_in_time=None, fade_in_and_out=None)
        Set the point sizes for clusters based on the number of stars.
    integrate_orbits(self, time, reference_frame_center=None)
        Integrates the orbits of the star cluster.
    limit_cluster_age(self, age_min, age_max)
        Limit the age of the star cluster.
    limit_cluster_by_name(self, names)
        Limit the star cluster by name.
    copy(self)
        Returns a copy of the Trace object.
    """

    def __init__(
        self,
        df,
        data_name,
        min_size=1,
        max_size=5,
        color='gray',
        opacity=1.0,
        marker_style='circle',
        show_tracks=False,
        size_by_n_stars=False,
        shifted_rf=None,
        fade_in_time=None,
        fade_in_and_out=None,
        colormap=None,
        cmin=None,
        cmax=None,
        fade_in_and_disp=None,
        disp_time=None,
    ):
        self.df = df
        self.data_name = data_name
        self.cluster_int_coords = None
        self.coordinates = None
        self.df_int = None
        self.integrated = False
        self.sizes_set = False
        self.color = color
        self.opacity = opacity
        self.marker_style = marker_style
        self.min_size = min_size
        self.max_size = max_size
        self.show_tracks = show_tracks
        self.size_by_n_stars = size_by_n_stars
        self.shifted_rf = shifted_rf
        self.fade_in_time = fade_in_time
        self.fade_in_and_out = fade_in_and_out
        self.fade_in_and_disp = fade_in_and_disp
        self.disp_time = disp_time
        self.colormap = colormap
        self.cmin = cmin
        self.cmax = cmax

        necessary_columns = ['x', 'y', 'z', 'U', 'V', 'W', 'name', 'age_myr']
        if self.size_by_n_stars:
            necessary_columns.append('n_stars')

        if not all(column in self.df.columns for column in necessary_columns):
            raise ValueError('User input data must contain the following columns: x, y, z, U, V, W, name, age_myr (and n_stars if size_by_n_stars=True)')

        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

    @staticmethod
    def _filter_integrated_component(component, member_mask):
        arr = np.asarray(component)

        if arr.ndim == 0:
            return arr

        if arr.shape[0] == member_mask.size:
            return arr[member_mask]

        if arr.ndim == 1 and member_mask.size and arr.size % member_mask.size == 0:
            reshaped = arr.reshape(member_mask.size, -1)
            return reshaped[member_mask]

        if member_mask.size == 1:
            return arr if member_mask[0] else arr[:0]

        raise ValueError('Integrated coordinate component shape is inconsistent with the current cluster membership mask.')

    def _filter_integrated_state(self, member_mask):
        if not self.integrated or self.cluster_int_coords is None:
            return

        self.cluster_int_coords = tuple(
            tuple(
                self._filter_integrated_component(component, member_mask)
                for component in coord_group
            )
            for coord_group in self.cluster_int_coords
        )

    def create_integrated_dataframe(self, time):
        """
        Create an integrated DataFrame of the star cluster.

        Parameters
        ----------
        time : array
            Array of time points.

        Returns
        -------
        df_int : pandas.DataFrame
            Integrated DataFrame of the star cluster.
        """
        if self.cluster_int_coords is None:
            print('Clusters have not yet been integrated')
            return

        xint, yint, zint = self.cluster_int_coords[0]
        xint_helio, yint_helio, zint_helio = self.cluster_int_coords[1]
        xint_gc, yint_gc, zint_gc = self.cluster_int_coords[2]
        rint_gc, phiint_gc, zint_gc = self.cluster_int_coords[3]

        df_int = pd.DataFrame({
            'x': xint.flatten(), 
            'y': yint.flatten(), 
            'z': zint.flatten(),
            'x_helio': xint_helio.flatten(),
            'y_helio': yint_helio.flatten(),
            'z_helio': zint_helio.flatten(),
            'x_gc': xint_gc.flatten(),
            'y_gc': yint_gc.flatten(),
            'z_gc': zint_gc.flatten(),
            'r_gc': rint_gc.flatten(),
            'phi_gc': np.rad2deg(phiint_gc.flatten()),
            'z_gc_cyl': zint_gc.flatten()
        })

        
        if 'n_stars' in self.df.columns:
            df_int['n_stars'] = np.repeat(self.df['n_stars'].values, len(time))

        df_int['age_myr'] = np.repeat(self.df['age_myr'].values, len(time))
        df_int['name'] = np.repeat(self.df['name'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))

        df_int.reset_index(drop=True, inplace=True)

        df_int = orbit_maker.coordFIX_to_coordROT(df_int)
        return df_int
    
    def set_age_based_sizes(self, fade_in_time=None, fade_in_and_out=None, fade_in_and_disp=None, disp_time=None):
        """
        Set the point sizes for clusters based on the number of stars.

        Parameters
        ----------
        fade_in_time : float, optional
            Duration of time for size easing. If None, uses instance value.
        fade_in_and_out : bool, optional
            Whether to fade in and out. If None, uses instance value.
        fade_in_and_disp : bool, optional
            Whether to fade in then drop to min size after disp_time. If None, uses instance value.
        disp_time : float, optional
            Time after birth to keep max size before dropping to min size. If None, uses instance value.
        """
        # Use instance values if parameters are None
        actual_fade_in_time = fade_in_time if fade_in_time is not None else self.fade_in_time
        actual_fade_in_and_out = fade_in_and_out if fade_in_and_out is not None else self.fade_in_and_out
        actual_fade_in_and_disp = fade_in_and_disp if fade_in_and_disp is not None else self.fade_in_and_disp
        actual_disp_time = disp_time if disp_time is not None else self.disp_time
        
        df_int_new_sizes = point_sizes.set_cluster_point_sizes(
            self.df_int, 
            min_size=self.min_size, 
            max_size=self.max_size, 
            fade_in_time=actual_fade_in_time, 
            fade_in_and_out=actual_fade_in_and_out, 
            size_by_n_stars=self.size_by_n_stars,
            fade_in_and_disp=actual_fade_in_and_disp,
            disp_time=actual_disp_time
        )
        self.df_int = df_int_new_sizes
        self.sizes_set = True

    def integrate_orbits(self, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
        """
        Integrate the orbits of the star cluster.

        Parameters
        ----------
        time : array
            Array of time points.
        reference_frame_center : tuple, optional
            Center of the reference frame.
        potential : galpy potential, optional
            Galactic potential to use. Defaults to MWPotential2014.
        vo : float, optional
            Circular velocity at the solar radius in km/s. Defaults to 236.
        ro : float, optional
            Solar radius in kpc. Defaults to 8.122.
        zo : float, optional
            Solar height above the galactic plane in kpc. Defaults to 0.0208.
        """

        if self.shifted_rf is not None:
            reference_frame_center = self.shifted_rf

        self.cluster_int_coords = orbit_maker.create_orbit(
            self.coordinates, time, 
            reference_frame_center=reference_frame_center,
            potential=potential,
            vo=vo, ro=ro, zo=zo
        )
        self.df_int = self.create_integrated_dataframe(time)
        self.integrated = True

    def limit_cluster_age(self, age_min, age_max):
        """
        Limit the age of the star cluster.

        Parameters
        ----------
        age_min : float
            Minimum age of the star cluster.
        age_max : float
            Maximum age of the star cluster.
        """
        member_mask = ((self.df['age_myr'] >= age_min) & (self.df['age_myr'] <= age_max)).to_numpy()
        self.df = self.df[member_mask]
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

        if self.integrated:
            self._filter_integrated_state(member_mask)
            self.df_int = self.df_int.loc[(self.df_int['age_myr'] >= age_min) & (self.df_int['age_myr'] <= age_max)]
    
    def limit_cluster_by_name(self, names):
        """
        Limit the star cluster by name.

        Parameters
        ----------
        names : list
            List of names of the star clusters to keep.
        """
        member_mask = self.df['name'].isin(names).to_numpy()
        self.df = self.df[member_mask]
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

        if self.integrated:
            self._filter_integrated_state(member_mask)
            self.df_int = self.df_int[self.df_int['name'].isin(names)]

    def copy(self):
        """
        Returns a copy of the Trace object.

        Returns
        -------
        Trace
            A copy of the Trace object.
        """
        return copy.deepcopy(self)

    @property
    def layer_name(self):
        return self.data_name

    @property
    def integrated_coords(self):
        return self.cluster_int_coords

    @property
    def has_time_varying_geometry(self):
        return True

    def limit_age(self, age_min, age_max):
        self.limit_cluster_age(age_min, age_max)

    def limit_by_name(self, names):
        self.limit_cluster_by_name(names)


class Layer(Trace):
    """
    General astronomy-facing spatial layer.

    `Layer` preserves the current `Trace` behavior but relaxes the dataframe
    requirements so callers can provide static 3D spatial data or full phase-space
    data for orbit tracing.
    """

    def __init__(
        self,
        df,
        layer_name,
        min_size=1,
        max_size=5,
        color='gray',
        opacity=1.0,
        marker_style='circle',
        show_tracks=False,
        size_by_n_stars=False,
        shifted_rf=None,
        fade_in_time=None,
        fade_in_and_out=None,
        colormap=None,
        cmin=None,
        cmax=None,
        fade_in_and_disp=None,
        disp_time=None,
        default_age_myr=0.0,
        assume_stationary=False,
    ):
        prepared_df, has_velocity_columns = self._normalize_layer_dataframe(
            df,
            layer_name=layer_name,
            size_by_n_stars=size_by_n_stars,
            default_age_myr=default_age_myr,
            assume_stationary=assume_stationary,
        )
        self._has_velocity_columns = has_velocity_columns

        super().__init__(
            prepared_df,
            data_name=layer_name,
            min_size=min_size,
            max_size=max_size,
            color=color,
            opacity=opacity,
            marker_style=marker_style,
            show_tracks=show_tracks,
            size_by_n_stars=size_by_n_stars,
            shifted_rf=shifted_rf,
            fade_in_time=fade_in_time,
            fade_in_and_out=fade_in_and_out,
            colormap=colormap,
            cmin=cmin,
            cmax=cmax,
            fade_in_and_disp=fade_in_and_disp,
            disp_time=disp_time,
        )

    @staticmethod
    def _normalize_layer_dataframe(df, *, layer_name, size_by_n_stars=False, default_age_myr=0.0, assume_stationary=False):
        if not isinstance(df, pd.DataFrame):
            raise TypeError('Layer input must be a pandas.DataFrame.')

        prepared_df = df.copy()
        missing_coordinates = [column for column in ('x', 'y', 'z') if column not in prepared_df.columns]
        if missing_coordinates:
            raise ValueError(
                'Layer input data must contain the following columns: x, y, z '
                f"(missing: {', '.join(missing_coordinates)})."
            )

        has_velocity_columns = all(column in prepared_df.columns for column in ('U', 'V', 'W'))
        if not has_velocity_columns:
            if not assume_stationary:
                raise ValueError(
                    'Layer input data must contain U, V, W to enable orbit tracing. '
                    'If this layer should stay static in time, pass assume_stationary=True.'
                )
            for column in ('U', 'V', 'W'):
                if column not in prepared_df.columns:
                    prepared_df[column] = 0.0

        if 'name' not in prepared_df.columns:
            prepared_df['name'] = str(layer_name)
        else:
            prepared_df['name'] = prepared_df['name'].fillna(str(layer_name))

        if 'age_myr' not in prepared_df.columns:
            prepared_df['age_myr'] = float(default_age_myr)
        else:
            prepared_df['age_myr'] = prepared_df['age_myr'].fillna(float(default_age_myr))

        if size_by_n_stars and 'n_stars' not in prepared_df.columns:
            raise ValueError('Layer input data must contain n_stars when size_by_n_stars=True.')

        return prepared_df, has_velocity_columns

    @classmethod
    def from_dataframe(cls, df, layer_name, **kwargs):
        return cls(df=df, layer_name=layer_name, **kwargs)

    @property
    def supports_orbit_tracing(self):
        return bool(self._has_velocity_columns)

    @property
    def has_time_varying_geometry(self):
        return self.supports_orbit_tracing

    def _static_cluster_int_coords(self, time):
        time = np.asarray(time, dtype=float)
        n_times = len(time)

        x = np.repeat(self.df['x'].to_numpy(dtype=float)[:, None], n_times, axis=1)
        y = np.repeat(self.df['y'].to_numpy(dtype=float)[:, None], n_times, axis=1)
        z = np.repeat(self.df['z'].to_numpy(dtype=float)[:, None], n_times, axis=1)
        r = np.hypot(x, y)
        phi = np.arctan2(y, x)

        return (
            (x, y, z),
            (x, y, z),
            (x, y, z),
            (r, phi, z),
        )

    def integrate_orbits(self, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
        if self.supports_orbit_tracing:
            return super().integrate_orbits(
                time,
                reference_frame_center=reference_frame_center,
                potential=potential,
                vo=vo,
                ro=ro,
                zo=zo,
            )

        self.cluster_int_coords = self._static_cluster_int_coords(time)
        self.df_int = self.create_integrated_dataframe(np.asarray(time, dtype=float))
        self.integrated = True

class TraceCollection:
    """
    Class for storing a collection of star clusters.

    Attributes
    ----------
    clusters : list
        List of Trace objects.

    Methods
    -------
    __init__(self, clusters=[])
        Initializes the TraceCollection object.
    add_cluster(self, cluster)
        Adds a star cluster to the collection.
    get_cluster(self, identifier)
        Retrieves a star cluster from the collection.
    get_all_clusters(self)
        Retrieves all star clusters from the collection.
    integrate_all_orbits(self, time, reference_frame_center=None)
        Integrates the orbits of all star clusters in the collection.
    set_all_cluster_sizes(self, fade_in_time, fade_in_and_out)
        Set the point sizes for all clusters in the collection.
    limit_all_cluster_ages(self, age_min, age_max)
        Limit the age of all star clusters in the collection.
    limit_all_cluster_names(self, names)
        Limit the star clusters in the collection by name.
    """

    def __init__(self, clusters=None):
        """
        Initialize the TraceCollection object.

        Parameters
        ----------
        clusters : list, optional
            List of Trace objects.

        Raises
        ------
        ValueError
            If any element in clusters is not a Trace object.
        """
        self.clusters = []
        self.time = None
        self.potential = None
        self.vo = None
        self.ro = None
        self.zo = None
        self.add_clusters(clusters or [])

    def add_clusters(self, clusters):
        """
        Add multiple star clusters to the collection.

        Parameters
        ----------
        clusters : list
            List of Trace objects.

        Raises
        ------
        ValueError
            If any element in clusters is not a Trace object.
        """
        for cluster in clusters:
            self.add_cluster(cluster)

    def add_cluster(self, cluster):
        """
        Add a star cluster to the collection.

        Parameters
        ----------
        cluster : Trace
            Trace object to be added.

        Raises
        ------
        ValueError
            If cluster is not a Trace object.
        """
        if not isinstance(cluster, Trace):
            raise ValueError('Input must be an instance of Trace')
        self.clusters.append(cluster)

    def get_cluster(self, identifier):
        """
        Retrieve a star cluster from the collection.

        Parameters
        ----------
        identifier : str or int
            Name or index of the star cluster.

        Returns
        -------
        cluster : Trace
            Trace object.

        Raises
        ------
        ValueError
            If identifier is neither a string (name) nor an integer (index).
        """
        if isinstance(identifier, str):
            for cluster in self.clusters:
                if cluster.data_name == identifier:
                    return cluster
        elif isinstance(identifier, int):
            return self.clusters[identifier]
        else:
            raise ValueError('Identifier must be a string (name) or an integer (index)')

    def get_all_clusters(self):
        """
        Retrieve all star clusters from the collection.

        Returns
        -------
        clusters : list
            List of Trace objects.
        """
        return self.clusters

    def integrate_all_orbits(self, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
        """
        Integrate the orbits of all star clusters in the collection.

        Parameters
        ----------
        time : array
            Array of time points.
        reference_frame_center : tuple, optional
            Center of the reference frame.
        potential : galpy potential, optional
            Galactic potential to use. Defaults to MWPotential2014.
        vo : float, optional
            Circular velocity at the solar radius in km/s. Defaults to 236.
        ro : float, optional
            Solar radius in kpc. Defaults to 8.122.
        zo : float, optional
            Solar height above the galactic plane in kpc. Defaults to 0.0208.
        """
        self.time = time
        self.potential = potential
        self.vo = vo
        self.ro = ro
        self.zo = zo
        for cluster in self.clusters:
            cluster.integrate_orbits(
                self.time, 
                reference_frame_center=reference_frame_center,
                potential=potential,
                vo=vo, ro=ro, zo=zo
            )
    
    def set_all_cluster_sizes(self, fade_in_time, fade_in_and_out, fade_in_and_disp=False, disp_time=0):
        """
        Set the point sizes for all clusters in the collection.

        Parameters
        ----------
        fade_in_time : int
            Duration of time for size easing (used as default for clusters without instance values).
        fade_in_and_out : bool
            Whether to fade in and out (used as default for clusters without instance values).
        fade_in_and_disp : bool
            Whether to fade in then drop to min size after disp_time (used as default for clusters without instance values).
        disp_time : float
            Time after birth to keep max size before dropping to min size (used as default for clusters without instance values).
        """
        for cluster in self.clusters:
            if not cluster.sizes_set:
                # Use cluster's instance values if they exist, otherwise use the provided defaults
                cluster_fade_in_time = cluster.fade_in_time if cluster.fade_in_time is not None else fade_in_time
                cluster_fade_in_and_out = cluster.fade_in_and_out if cluster.fade_in_and_out is not None else fade_in_and_out
                cluster_fade_in_and_disp = cluster.fade_in_and_disp if cluster.fade_in_and_disp is not None else fade_in_and_disp
                cluster_disp_time = cluster.disp_time if cluster.disp_time is not None else disp_time
                cluster.set_age_based_sizes(
                    cluster_fade_in_time,
                    cluster_fade_in_and_out,
                    cluster_fade_in_and_disp,
                    cluster_disp_time
                )

    def limit_all_cluster_ages(self, age_min, age_max):
        """
        Limit the age of all star clusters in the collection.

        Parameters
        ----------
        age_min : float
            Minimum age of the star cluster.
        age_max : float
            Maximum age of the star cluster.
        """
        for cluster in self.clusters:
            cluster.limit_cluster_age(age_min, age_max)
    
    def limit_all_cluster_names(self, names):
        """
        Limit the star clusters in the collection by name.

        Parameters
        ----------
        names : list
            List of names of the star clusters to keep.
        """
        for cluster in self.clusters:
            cluster.limit_cluster_by_name(names)

    @property
    def layers(self):
        return self.clusters

    def add_layer(self, layer):
        self.add_cluster(layer)

    def get_layer(self, identifier):
        return self.get_cluster(identifier)

    def get_all_layers(self):
        return self.get_all_clusters()

    def integrate_all_layers(self, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
        self.integrate_all_orbits(
            time,
            reference_frame_center=reference_frame_center,
            potential=potential,
            vo=vo,
            ro=ro,
            zo=zo,
        )

    def set_all_layer_sizes(self, fade_in_time, fade_in_and_out, fade_in_and_disp=False, disp_time=0):
        self.set_all_cluster_sizes(
            fade_in_time,
            fade_in_and_out,
            fade_in_and_disp=fade_in_and_disp,
            disp_time=disp_time,
        )


class LayerCollection(TraceCollection):
    pass
