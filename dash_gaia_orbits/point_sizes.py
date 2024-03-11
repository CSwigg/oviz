import numpy as np
import pandas as pd

def size_easing(c, min_size, max_size, time_duration=5):
    """
    Apply size easing to the sizes of a cluster over time.

    Args:
        c (pd.DataFrame): Data frame containing the data.
        min_size (float): Minimum size value.
        max_size (float): Maximum size value.
        time_duration (int, optional): Duration of time for size easing. Defaults to 5.

    Returns:
        pd.DataFrame: Data frame with updated 'size' column.
    """

    t = -1 * c['time'].values
    age = c['age_myr'].values[0]

    a = min_size
    b = max_size

    sizes = np.array([max_size] * len(c))

    t_older = t[np.where((t > age) & (t <= age + time_duration))]
    sizes_younger = sizes[np.where(t <= age)]
    sizes_older = sizes[np.where((t > age) & (t <= (age + time_duration)))]
    sizes_oldest = sizes[np.where(t > age + time_duration)]
    sizes_oldest = [a] * len(sizes_oldest)

    w = 0.5
    D = np.linspace(2, 0, len(t_older))
    sigmaD = 1.0 / (1.0 + np.exp(-(1 - D) / w))
    sizes_older = a + (b - a) * (1 - sigmaD)

    c['size'] = np.concatenate([sizes_younger, sizes_older, sizes_oldest])

    return c


def set_cluster_point_sizes(df_int, min_size, max_size):
    """
    Set point sizes for clusters in the given data frame.

    Args:
        df_int (pd.DataFrame): Data frame containing the data.
        min_size (float): Minimum size value.
        max_size (float): Maximum size value.

    Returns:
        pd.DataFrame: Data frame with updated point sizes.
    """

    if len(df_int) == 0:
        return df_int

    def compute_size(group):
        age = group['age_myr'].mean()
        group = size_easing(group, min_size, max_size)
        return group

    df_int = df_int.groupby('name').apply(compute_size).sort_index()
    return df_int