"""Compute time-dependent marker sizes and opacities for cluster traces."""

import numpy as np


def _clamp01(values):
    return np.clip(values, 0.0, 1.0)


def _ease(t, birth, low, high, fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time):
    """One object's values at times ``t``: ``low`` before ``birth``, ``high`` after.

    The change is a sigmoid over the ``fade_in_time`` Myr ending at ``birth``.
    ``fade_in_and_out`` eases back down over the next ``fade_in_time`` Myr;
    ``fade_in_and_disp`` holds ``high`` for ``disp_time`` Myr, then drops to
    ``low``. Values are ordered for a descending grid (0, -1, ...) when no
    time is positive, and for an ascending grid otherwise.
    """
    if fade_in_and_out and fade_in_and_disp:
        raise ValueError("fade_in_and_out and fade_in_and_disp cannot both be True")

    hold_time = disp_time if fade_in_and_disp else fade_in_time
    n_before = np.count_nonzero(t < birth - fade_in_time)
    n_birth = np.count_nonzero((t >= birth - fade_in_time) & (t <= birth))
    n_after = np.count_nonzero((t > birth) & (t <= birth + hold_time))
    n_future = np.count_nonzero(t > birth + hold_time)

    values_before = [low] * n_before
    values_future = [high] * n_future

    w = 0.5
    past_only = bool(np.all(t <= 0))
    D = np.linspace(2, 0, n_birth) if past_only else np.linspace(0, 2, n_birth)
    sigmaD = 1.0 / (1.0 + np.exp(-(1 - D) / w))
    values_birth = low + (high - low) * (1 - sigmaD)

    if fade_in_and_out:
        D2 = np.linspace(2, 0, n_after)
        sigmaD2 = 1.0 / (1.0 + np.exp(-(1 - D2) / w))
        values_after = low + (high - low) * (1 - sigmaD2)
        values_future = [low] * n_future
    elif fade_in_and_disp:
        values_after = [high] * n_after
        values_future = [low] * n_future
    else:
        values_after = [high] * n_after

    if past_only:
        if fade_in_and_out:
            values_after = np.flip(values_after)
        return np.concatenate([values_future, values_after, values_birth, values_before])
    return np.concatenate([values_before, values_birth, values_after, values_future])


def _star_count_weight(n_stars):
    """Size weight of a cluster with ``n_stars`` members: n / 500, kept within [0.1, 1]."""
    return min(max(n_stars / 500, 0.1), 1)


def size_easing(
    c,
    min_size,
    max_size,
    size_by_n_stars,
    fade_in_time=5,
    fade_in_and_out=False,
    fade_in_and_disp=False,
    disp_time=0
):
    """Set the ``size`` column of one cluster's rows (``time``, ``age_myr``).

    Sizes ease from ``min_size`` before the cluster's birth to ``max_size``
    after it. ``fade_in_and_out`` shrinks them again over the next
    ``fade_in_time`` Myr; ``fade_in_and_disp`` keeps ``max_size`` for
    ``disp_time`` Myr and then drops to ``min_size``. With ``size_by_n_stars``
    both bounds scale with ``n_stars`` (see the first row). Returns ``c``.
    """
    if size_by_n_stars:
        weight = _star_count_weight(c['n_stars'].values[0])
        min_size *= weight
        max_size *= weight
    c['size'] = _ease(
        c['time'].values, -c['age_myr'].values[0], min_size, max_size,
        fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time,
    )
    return c


def opacity_easing(
    c,
    min_opacity,
    max_opacity,
    fade_in_time=5,
    fade_in_and_out=False,
    fade_in_and_disp=False,
    disp_time=0
):
    """Set the ``opacity`` column of one cluster's rows, eased like :func:`size_easing`.

    Opacities are kept within [0, 1]. Returns ``c``.
    """
    low = float(_clamp01(min_opacity))
    high = float(_clamp01(max_opacity))
    c['opacity'] = _ease(
        c['time'].values, -c['age_myr'].values[0], low, high,
        fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time,
    )
    c['opacity'] = _clamp01(c['opacity'].values)
    return c


def _per_cluster(df_int, column, values_for, needed_columns):
    """A copy of ``df_int`` with ``column`` computed for each ``name`` group.

    ``values_for(rows, *arrays)`` gets one group's row positions and the
    ``needed_columns`` as arrays. Rows without a name are dropped, as by
    ``groupby('name')``; the others keep index order.
    """
    if len(df_int) == 0:
        return df_int
    if 'name' in df_int.index.names:
        df_int = df_int.reset_index(drop=True)

    groups = df_int.groupby('name', sort=False).indices
    arrays = [df_int[name].values for name in needed_columns]
    values = np.empty(len(df_int))
    grouped = np.zeros(len(df_int), dtype=bool)
    for rows in groups.values():
        values[rows] = values_for(rows, *arrays)
        grouped[rows] = True

    result = df_int.copy() if grouped.all() else df_int[grouped]
    if groups:
        result[column] = values[grouped]
    return result.sort_index()


def set_cluster_point_sizes(
    df_int,
    min_size,
    max_size,
    fade_in_time,
    fade_in_and_out,
    size_by_n_stars,
    fade_in_and_disp=False,
    disp_time=0
):
    """Return ``df_int`` with a ``size`` column eased per cluster (see :func:`size_easing`).

    ``df_int`` holds the rows of every cluster (``name``, ``time``,
    ``age_myr`` and, with ``size_by_n_stars``, ``n_stars``), each cluster's
    rows in timeline order.
    """
    def sizes(rows, time, age, n_stars=None):
        low, high = min_size, max_size
        if size_by_n_stars:
            weight = _star_count_weight(n_stars[rows[0]])
            low *= weight
            high *= weight
        return _ease(time[rows], -age[rows[0]], low, high, fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time)

    needed = ('time', 'age_myr', 'n_stars') if size_by_n_stars else ('time', 'age_myr')
    return _per_cluster(df_int, 'size', sizes, needed)


def set_cluster_point_opacities(
    df_int,
    min_opacity,
    max_opacity,
    fade_in_time,
    fade_in_and_out,
    fade_in_and_disp=False,
    disp_time=0
):
    """Return ``df_int`` with an ``opacity`` column eased per cluster (see :func:`opacity_easing`)."""
    def opacities(rows, time, age):
        low = float(_clamp01(min_opacity))
        high = float(_clamp01(max_opacity))
        return _clamp01(_ease(time[rows], -age[rows[0]], low, high, fade_in_time, fade_in_and_out, fade_in_and_disp, disp_time))

    return _per_cluster(df_int, 'opacity', opacities, ('time', 'age_myr'))
