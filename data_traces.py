# Classes for storing user input data
import pandas as pd


class StarClusterData:

    def __init__(self, df):
        self.df = df
        try:
            # User input must contain the following columns
            assert all(column in self.df.columns for column in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'name', 'age_myr'])
        except AssertionError:
            raise ValueError('User input data must contain the following columns: x, y, z, vx, vy, vz, name, age_myr')



