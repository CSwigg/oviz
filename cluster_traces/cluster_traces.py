# Classes for storing user input data
import numpy as np
import pandas as pd
from cluster_traces import orbit_maker

class StarClusterData:
    """
    Class for storing user input data of a star cluster.

    Attributes:
    - df (pandas.DataFrame): Input data containing the following columns: x, y, z, U, V, W, name, age_myr.
    - data_name (str): Name of the star cluster data.
    - cluster_int_coords (tuple): Integrated coordinates of the star cluster.
    - coordinates (numpy.ndarray): Array of coordinates (x, y, z, U, V, W) of the star cluster.

    Methods:
    - __init__(self, df, data_name): Initializes the StarClusterData object.
    - create_integrated_dataframe(self, time): Creates an integrated DataFrame of the star cluster.
    - integrate_orbits(self, time): Integrates the orbits of the star cluster.
    """

    def __init__(self, df, data_name):
        """
        Initialize the StarClusterData object.

        Parameters:
        - df (pandas.DataFrame): Input data containing the following columns: x, y, z, U, V, W, name, age_myr.
        - data_name (str): Name of the star cluster data.

        Raises:
        - ValueError: If the input data does not contain the required columns.
        """
        self.df = df
        self.data_name = data_name
        self.cluster_int_coords = None
        self.coordinates = None
        self.df_int = None

        try:
            # User input must contain the following columns
            assert all(column in self.df.columns for column in ['x', 'y', 'z', 'U', 'V', 'W', 'name', 'age_myr'])
        except AssertionError:
            raise ValueError('User input data must contain the following columns: x, y, z, U, V, W, name, age_myr')

        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values
    

    def create_integrated_dataframe(self, time):
        """
        Create an integrated DataFrame of the star cluster.

        Parameters:
        - time (float): The time at which the orbits are integrated.

        Returns:
        - df_int (pandas.DataFrame): Integrated DataFrame of the star cluster.
        """
        if self.cluster_int_coords is None:
            print('Clusters have not yet been integrated')
            return

        xint, yint, zint = self.cluster_int_coords
        df_int = pd.DataFrame({'x': xint.flatten(), 'y': yint.flatten(), 'z': zint.flatten()})
        df_int['name'] = np.tile(self.df['name'].values, len(time))
        df_int['age_myr'] = np.tile(self.df['age_myr'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))
        return df_int


    def integrate_orbits(self, time):
        """
        Integrate the orbits of the star cluster.

        Parameters:
        - time (float): The time at which the orbits are integrated.
        """
        self.cluster_int_coords = orbit_maker.create_orbit(self.coordinates, time)
        self.df_int = self.create_integrated_dataframe(time) # Create the integrated DataFrame


class StarClusterCollection:
    """
    Class for storing a collection of star clusters.

    Attributes:
    - clusters (list): List of StarClusterData objects.

    Methods:
    - __init__(self, clusters=[]): Initializes the StarClusterCollection object.
    - add_cluster(self, cluster): Adds a star cluster to the collection.
    - get_cluster(self, identifier): Retrieves a star cluster from the collection.
    - get_all_clusters(self): Retrieves all star clusters from the collection.
    """

    def __init__(self, clusters=[]):
        """
        Initialize the StarClusterCollection object.

        Parameters:
        - clusters (list): List of StarClusterData objects.

        Raises:
        - ValueError: If any element in clusters is not a StarClusterData object.
        """
        if len(clusters) != 0:
            for cluster in clusters:
                if not isinstance(cluster, StarClusterData):
                    raise ValueError('Input must be a list of StarClusterData objects')
        self.clusters = clusters

    def add_cluster(self, cluster):
        """
        Add a star cluster to the collection.

        Parameters:
        - cluster (StarClusterData): StarClusterData object to be added.

        Raises:
        - ValueError: If cluster is not a StarClusterData object.
        """
        if not isinstance(cluster, StarClusterData):
            raise ValueError('Input must be an instance of StarClusterData')
        self.clusters.append(cluster)

    def get_cluster(self, identifier):
        """
        Retrieve a star cluster from the collection.

        Parameters:
        - identifier (str or int): Name or index of the star cluster.

        Returns:
        - cluster (StarClusterData): StarClusterData object.

        Raises:
        - ValueError: If identifier is neither a string (name) nor an integer (index).
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

        Returns:
        - clusters (list): List of StarClusterData objects.
        """
        return self.clusters


    def integrate_all_orbits(self, time):
        """
        Integrate the orbits of all star clusters in the collection.

        Parameters:
        - time (float): The time at which the orbits are integrated.
        """
        for cluster in self.clusters:
            cluster.integrate_orbits(time)
    
    def integrate_all_orbits(self, time):
        for cluster in self.clusters:
            cluster.integrate_orbits(time)