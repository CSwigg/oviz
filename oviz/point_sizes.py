import numpy as np
import pandas as pd

def size_easing(c, min_size, max_size, size_by_n_stars, fade_in_time=5, fade_in_and_out=False):
    """
    Apply size easing to the sizes of a cluster over time.

    Parameters:
    - c (pd.DataFrame): Data frame containing the data.
    - min_size (float): Minimum size value.
    - max_size (float): Maximum size value.
    - size_by_n_stars (bool): Whether to size by number of stars.
    - fade_in_time (int, optional): Duration of time for size easing. Defaults to 5.
    - fade_in_and_out (bool, optional): Whether to fade in and out. Defaults to False.

    Returns:
    pd.DataFrame: Data frame with updated 'size' column.
    """
    #t = -1 * c['time'].values
    t = c['time'].values
    age = -c['age_myr'].values[0]
    

    if size_by_n_stars:
        n_stars = c['n_stars'].values[0]
        weight = n_stars / 500
        weight = min(max(weight, 0.1), 1)
        min_size *= weight
        max_size *= weight

    a = min_size
    b = max_size

    #sizes = np.array([max_size] * len(c))
    ind_before_birth = np.where(t < age - fade_in_time)
    ind_birth = np.where((t >= age - fade_in_time) & (t <= age))
    ind_after_birth = np.where((t > age) & (t <= age + fade_in_time))
    ind_future = np.where((t > age + fade_in_time))


    times_before= t[ind_before_birth]
    times_birth = t[ind_birth]
    times_after = t[ind_after_birth]
    times_future = t[ind_future]

    sizes_before = [a] * len(times_before)
    sizes_future = [b] * len(times_future)

    w = 0.5
    if all(t <= 0):
        D = np.linspace(2, 0, len(times_birth))
    else:
        D = np.linspace(0, 2, len(times_birth))
    sigmaD = 1.0 / (1.0 + np.exp(-(1 - D) / w))
    sizes_birth = a + (b - a) * (1 - sigmaD)

    if fade_in_and_out:
        D2 = np.linspace(2, 0, len(times_after))
        sigmaD2 = 1.0 / (1.0 + np.exp(-(1 - D2) / w))
        sizes_after = a + (b - a) * (1 - sigmaD2)
        sizes_future = [a] * len(times_future)
    else:
        sizes_after = [b] * len(times_after)

    if all(t <= 0):
        c['size'] = np.concatenate([sizes_future, sizes_after, sizes_birth, sizes_before])
    else:
        c['size'] = np.concatenate([sizes_before, sizes_birth, sizes_after, sizes_future])
    return c

def set_cluster_point_sizes(df_int, min_size, max_size, fade_in_time, fade_in_and_out, size_by_n_stars):
    """
    Set point sizes for clusters in the given data frame.

    Parameters:
    - df_int (pd.DataFrame): Data frame containing the data.
    - min_size (float): Minimum size value.
    - max_size (float): Maximum size value.
    - fade_in_time (int): Duration of time for size easing.
    - fade_in_and_out (bool): Whether to fade in and out.
    - size_by_n_stars (bool): Whether to size by number of stars.

    Returns:
    pd.DataFrame: Data frame with updated point sizes.
    """
    if len(df_int) == 0:
        return df_int

    def compute_size(group, size_by_n_stars):
        group = size_easing(group, min_size, max_size, size_by_n_stars, fade_in_time, fade_in_and_out)
        return group

    df_int = df_int.groupby('name').apply(compute_size, size_by_n_stars).sort_index()
    return df_int