# Classes for storing user input data
import numpy as np
import pandas as pd
from . import orbit_maker
from . import point_sizes
import copy

class StarClusterData:
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
    opacity : float, optional
        Opacity of the marker points (default is 1.0).

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
    opacity : float
        Opacity of the marker points.
    min_size : int
        Minimum size of the marker points.
    max_size : int
        Maximum size of the marker points.

    Methods
    -------
    __init__(self, df, data_name, min_size=1, max_size=5, color='gray', opacity=1.0)
        Initializes the StarClusterData object.
    create_integrated_dataframe(self, time)
        Creates an integrated DataFrame of the star cluster.
    set_age_based_sizes(self, use_n_stars=False)
        Set the point sizes for clusters based on the number of stars.
    integrate_orbits(self, time)
        Integrates the orbits of the star cluster.
    """

    def __init__(self, df, data_name, min_size=1, max_size=5, color='gray', opacity=1.0, marker_style = 'circle'):
        """
        Initialize the StarClusterData object.

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
        opacity : float, optional
            Opacity of the marker points (default is 1.0).

        Raises
        ------
        ValueError
            If the input data does not contain the required columns.
        """
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

        # copies are made so that the original data is still stored with the limit_cluster_age method
        #self.df_original = self.df.copy() # Save a copy of the original DataFrame
        #self.df_int_original = self.df_int.copy()

        try:
            # User input must contain the following columns
            assert all(column in self.df.columns for column in ['x', 'y', 'z', 'U', 'V', 'W', 'name', 'age_myr'])
        except AssertionError:
            raise ValueError('User input data must contain the following columns: x, y, z, U, V, W, name, age_myr')

        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values
    

    def create_integrated_dataframe(self, time):
        """
        Create an integrated DataFrame of the star cluster.

        Parameters
        ----------
        time : float
            The time at which the orbits are integrated.

        Returns
        -------
        df_int : pandas.DataFrame
            Integrated DataFrame of the star cluster.
        """
        if self.cluster_int_coords is None:
            print('Clusters have not yet been integrated')
            return

        xint, yint, zint = self.cluster_int_coords
        df_int = pd.DataFrame({'x': xint.flatten(), 'y': yint.flatten(), 'z': zint.flatten()})
        df_int['age_myr'] = np.repeat(self.df['age_myr'].values, len(time))
        df_int['name'] = np.repeat(self.df['name'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))
        return df_int
    
    def set_age_based_sizes(self, fade_in_time):
        """
        Set the point sizes for clusters based on the number of stars.

        Parameters
        ----------
        use_n_stars : bool, optional
            Flag indicating whether to scale based on the number of stars (default is False).

        Returns
        -------
        df_int : pandas.DataFrame
            Updated dataframe with scaled point sizes.
        """
        df_int_new_sizes = point_sizes.set_cluster_point_sizes(self.df_int, min_size=self.min_size, max_size=self.max_size, fade_in_time = fade_in_time)
        self.df_int = df_int_new_sizes
        self.sizes_set = True
        


    def integrate_orbits(self, time, reference_frame_center = None):
        """
        Integrate the orbits of the star cluster.

        Parameters
        ----------
        time : float
            The time at which the orbits are integrated.
        """
        self.cluster_int_coords = orbit_maker.create_orbit(self.coordinates, time, reference_frame_center = reference_frame_center)
        self.df_int = self.create_integrated_dataframe(time) # Create the integrated DataFrame
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
        #self.df = self.df_original.copy()
        self.df = self.df[(self.df['age_myr'] >= age_min) & (self.df['age_myr'] <= age_max)]
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

        if self.integrated:
            self.df_int = self.df_int.loc[(self.df_int['age_myr'] >= age_min) & (self.df_int['age_myr'] <= age_max)]
            self.cluster_int_coords = (self.df_int['x'].values, self.df_int['y'].values, self.df_int['z'].values)
    
    def limit_cluster_by_name(self, names):
        """
        Limit the star cluster by name.

        Parameters
        ----------
        names : list
            List of names of the star clusters to keep.
        """
        self.df = self.df[self.df['name'].isin(names)]
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

        if self.integrated:
            self.df_int = self.df_int[self.df_int['name'].isin(names)]
            self.cluster_int_coords = (self.df_int['x'].values, self.df_int['y'].values, self.df_int['z'].values)


    def copy(self):
        '''
        Returns a copy of the StarClusterData object
        '''
        return copy.deepcopy(self)

class StarClusterCollection:
    """
    Class for storing a collection of star clusters.

    Attributes
    ----------
    clusters : list
        List of StarClusterData objects.

    Methods
    -------
    __init__(self, clusters=[])
        Initializes the StarClusterCollection object.
    add_cluster(self, cluster)
        Adds a star cluster to the collection.
    get_cluster(self, identifier)
        Retrieves a star cluster from the collection.
    get_all_clusters(self)
        Retrieves all star clusters from the collection.
    """

    def __init__(self, clusters=[]):
        """
        Initialize the StarClusterCollection object.

        Parameters
        ----------
        clusters : list, optional
            List of StarClusterData objects.

        Raises
        ------
        ValueError
            If any element in clusters is not a StarClusterData object.
        """
        if len(clusters) != 0:
            for cluster in clusters:
                if not isinstance(cluster, StarClusterData):
                    raise ValueError('Input must be a list of StarClusterData objects')
        self.clusters = clusters
        self.time = None

    def add_cluster(self, cluster):
        """
        Add a star cluster to the collection.

        Parameters
        ----------
        cluster : StarClusterData
            StarClusterData object to be added.

        Raises
        ------
        ValueError
            If cluster is not a StarClusterData object.
        """
        if not isinstance(cluster, StarClusterData):
            raise ValueError('Input must be an instance of StarClusterData')
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
        cluster : StarClusterData
            StarClusterData object.

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
            List of StarClusterData objects.
        """
        return self.clusters

    def integrate_all_orbits(self, time, reference_frame_center = None):
        """
        Integrate the orbits of all star clusters in the collection.

        Parameters
        ----------
        time : float
            The time at which the orbits are integrated.
        """
        self.time = time # Makes sure each cluster group is integraetd with the same time array
        for cluster in self.clusters:
            cluster.integrate_orbits(self.time, reference_frame_center = reference_frame_center)
    
    def set_all_cluster_sizes(self, fade_in_time):
        """
        Set the point sizes for all clusters in the collection.
        """
        for cluster in self.clusters:
            if cluster.sizes_set:
                # don't reset the sizes if they've already been set
                continue
            else:
                cluster.set_age_based_sizes(fade_in_time)

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