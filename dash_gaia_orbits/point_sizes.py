import numpy as np
import pandas as pd

def size_easing(c, min_size, max_size, size_by_n_stars, fade_in_time=5, fade_in_and_out=False):
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
    

    if size_by_n_stars:
        n_stars = c['n_stars'].values[0]
        weight = n_stars/500
        if weight > 1:
            weight = 1
        if weight < 0.1:
            weight = 0.1
        min_size = min_size*weight
        max_size = max_size*weight

    a = min_size
    b = max_size

    sizes = np.array([max_size] * len(c))

    t_older = t[np.where((t > age) & (t <= age + fade_in_time))]
    t_younger = t[np.where(t <= age)]
    sizes_younger = sizes[np.where(t <= age)]
    sizes_older = sizes[np.where((t > age) & (t <= (age + fade_in_time)))]
    sizes_oldest = sizes[np.where(t > age + fade_in_time)]
    sizes_oldest = [a] * len(sizes_oldest)

    w = 0.5
    D = np.linspace(2, 0, len(t_older))
    sigmaD = 1.0 / (1.0 + np.exp(-(1 - D) / w))
    sizes_older = a + (b - a) * (1 - sigmaD)

    if fade_in_and_out:
        sizes_oldest = [0] * len(sizes_oldest)
        t_younger = t[np.where((t <= age) & (t >= age - fade_in_time))]
        t_youngest = t[np.where(t < age - fade_in_time)]
        D2 = np.linspace(0, 2, len(t_younger))
        sigmaD2 = 1.0 / (1.0 + np.exp(-(1 - D2) / w))
        sizes_younger = a + (b - a) * (1 - sigmaD2)
        sizes_youngest = [a]*len(t_youngest)
        c['size'] = np.concatenate([sizes_youngest, sizes_younger, sizes_older, sizes_oldest])


    else:
        c['size'] = np.concatenate([sizes_younger, sizes_older, sizes_oldest])

    return c


def set_cluster_point_sizes(df_int, min_size, max_size, fade_in_time, fade_in_and_out, size_by_n_stars):
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


    def compute_size(group, size_by_n_stars):
        age = group['age_myr'].mean()
        group = size_easing(c = group, min_size = min_size, max_size = max_size, fade_in_time=fade_in_time, fade_in_and_out = fade_in_and_out, size_by_n_stars = size_by_n_stars)
        return group

    df_int = df_int.groupby('name').apply(compute_size, size_by_n_stars).sort_index()

    return df_int