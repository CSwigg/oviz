import numpy as np
import pandas as pd

def set_cluster_point_sizes(df_int, min_size, max_size):
    def compute_size(group):
        age = group['age_myr'].mean()
        # size = np.where(group['time'] >= -age, initial_size, min_size)
        # size = np.where((group['time'] < -age) & (group['time'] >= -age - 3), size * np.exp(group['time'] + age), size)

        size_small = min_size
        size_large = max_size

        group['size'] = size_large
        group.loc[group['time'] < -age, 'size'] = size_small


        return group

    df_int = df_int.groupby('name').apply(compute_size).sort_index()
    return df_int