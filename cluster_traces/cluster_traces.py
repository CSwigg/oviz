# Classes for storing user input data
import numpy as np
import pandas as pd
#from orbit_maker import create_orbit
#import orbit_maker
from cluster_traces import orbit_maker

class StarClusterData:

    def __init__(self, df):
        """
        Initialize the ClusterTraces object.

        Parameters:
        df (pandas.DataFrame): Input data containing the following columns: x, y, z, vx, vy, vz, name, age_myr.

        Raises:
        ValueError: If the input data does not contain the required columns.

        """
        self.df = df
        self.cluster_int_coords = None
        self.coordinates = None
        self.df_int = None

        try:
            # User input must contain the following columns
            assert all(column in self.df.columns for column in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'name', 'age_myr'])
        except AssertionError:
            raise ValueError('User input data must contain the following columns: x, y, z, vx, vy, vz, name, age_myr')

        self.coordinates = self.df[['x', 'y', 'z', 'vx', 'vy', 'vz']].T.values
    

    def create_integrated_dataframe(self, time):
        if self.cluster_int_coords is None:
            print('Clusters have not yet been integrated')
            return

        xint, yint, zint = self.cluster_int_coords
        df_int = pd.DataFrame({'x': xint.flatten(), 'y': yint.flatten(), 'z': zint.flatten()})
        df_int['name'] = np.tile(self.df['name'].values, len(time))
        df_int['age_myr'] = np.tile(self.df['age_myr'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))
        # df_int['age_myr'] = np.tile(self.df['age_myr'].values, (len(time), 1)).T
        # df_int['time'] = np.tile(time, (len(self.df), 1))
        return df_int


    def integrate_orbits(self, time):
        """
        Integrate the orbits of the star clusters.

        Parameters:
        time (float): The time at which the orbits are integrated.
        """
        self.cluster_int_coords = orbit_maker.create_orbit(self.coordinates, time)
        self.df_int = self.create_integrated_dataframe(time) # Create the integrated DataFrame