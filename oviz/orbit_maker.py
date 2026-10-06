"""Integrate Galactic orbits and convert them to Oviz coordinate frames."""

import math
import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014

warnings.filterwarnings("ignore")


def normalize_time_grid(time, *, require_zero=True):
    """Return ``time`` as a float array, in the caller's order, after validation.

    The grid must be one-dimensional, finite and free of duplicates, and must
    contain 0 (the present) unless ``require_zero`` is False.
    """
    if time is None:
        raise ValueError("time must be a non-empty one-dimensional array-like.")
    try:
        time_array = np.asarray(time, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("time must contain only numeric values.") from exc
    if time_array.ndim != 1 or time_array.size == 0:
        raise ValueError("time must be a non-empty one-dimensional array-like.")
    if not np.all(np.isfinite(time_array)):
        raise ValueError("time must contain only finite values.")
    if np.unique(time_array).size != time_array.size:
        raise ValueError("time must not contain duplicate values.")
    if require_zero and not np.any(time_array == 0.0):
        raise ValueError("The time array must include 0 for the present day.")
    return np.array(time_array, dtype=float, copy=True)


def _integrate_orbit_in_requested_order(orbit, time, potential):
    """Integrate outward from zero, backward and forward separately.

    Returns ``values(accessor)``, which applies ``accessor`` to each branch's
    SkyCoord and reassembles the result (time on the last axis) in the
    caller's order.
    """
    requested_time = normalize_time_grid(time)
    negative_time = np.sort(requested_time[requested_time < 0.0])
    positive_time = np.sort(requested_time[requested_time > 0.0])
    branches = []

    if negative_time.size:
        integration_time = np.concatenate(([0.0], negative_time[::-1]))
        negative_orbit = deepcopy(orbit)
        negative_orbit.integrate(integration_time * u.Myr, potential)
        negative_coords = negative_orbit.SkyCoord(integration_time * u.Myr)
        branches.append((np.concatenate((negative_time, [0.0])), negative_coords, True, False))

    if positive_time.size:
        integration_time = np.concatenate(([0.0], positive_time))
        positive_orbit = deepcopy(orbit)
        positive_orbit.integrate(integration_time * u.Myr, potential)
        positive_coords = positive_orbit.SkyCoord(integration_time * u.Myr)
        branches.append((integration_time, positive_coords, False, bool(negative_time.size)))
    elif not negative_time.size:
        present_time = np.array([0.0])
        branches.append((present_time, orbit.SkyCoord(present_time * u.Myr), False, False))

    canonical_time = np.concatenate([branch[0][1:] if branch[3] else branch[0] for branch in branches])
    canonical_index = {float(value): idx for idx, value in enumerate(canonical_time)}
    requested_index = np.array([canonical_index[float(value)] for value in requested_time], dtype=int)

    def values(accessor):
        pieces = []
        for _branch_time, coords, reverse, drop_zero in branches:
            value = np.atleast_1d(np.asarray(accessor(coords)))
            if reverse:
                value = value[..., ::-1]
            if drop_zero:
                value = value[..., 1:]
            pieces.append(value)
        canonical_values = np.concatenate(pieces, axis=-1)
        return canonical_values[..., requested_index]

    return values


def reference_potential(potential, reference_frame_center=None):
    """Return the potential a reference-frame orbit is integrated in.

    The default frame follows the Local Standard of Rest, a circular orbit,
    which only the axisymmetric part of the potential has: non-axisymmetric
    terms (spiral arms, a bar) are left out of it, so the LSR stays at the
    origin of a frame whose Galactic centre circles it at R0. A custom frame
    centre (a focus group's orbit) is an object and feels the whole potential.
    """
    if potential is None:
        return MWPotential2014
    if reference_frame_center is not None:
        return potential
    parts = list(potential) if isinstance(potential, (list, tuple)) else [potential]
    axisymmetric = [part for part in parts if not getattr(part, "isNonAxi", False)]
    if not axisymmetric or len(axisymmetric) == len(parts):
        return potential
    return axisymmetric


def get_center_orbit_coords(time, reference_frame_center, potential=None, vo=236., ro=8.122, zo=0.0208):
    """Heliocentric Galactic x, y, z (pc) of the reference-frame orbit at each time.

    Parameters
    ----------
    time : array-like
        Timeline in Myr; must include 0.
    reference_frame_center : sequence of 6 floats or None
        ``x, y, z`` (pc) and ``U, V, W`` (km/s) of the frame's centre today.
        ``None`` follows the Local Standard of Rest.
    potential : galpy potential, optional
        Defaults to ``MWPotential2014``; the LSR frame uses only its
        axisymmetric part (see :func:`reference_potential`).
    vo, ro, zo : float
        Circular velocity (km/s), solar radius (kpc) and solar height (kpc).
    """
    time = normalize_time_grid(time)
    potential = reference_potential(potential, reference_frame_center)

    if reference_frame_center is None:
        rf_coords = [0, 0, 0, -11.1, -12.24, -7.25]
    else:
        rf_coords = reference_frame_center

    rf_sc = SkyCoord(
        u=rf_coords[0]*u.pc, v=rf_coords[1]*u.pc, w=rf_coords[2]*u.pc,
        U=rf_coords[3]*u.km/u.s, V=rf_coords[4]*u.km/u.s, W=rf_coords[5]*u.km/u.s,
        frame='galactic', representation_type='cartesian', differential_type='cartesian'
    )
    rf_orbit = Orbit(vxvv=rf_sc, solarmotion='schoenrich', ro=ro, vo=vo, zo=zo)

    integrated_values = _integrate_orbit_in_requested_order(rf_orbit, time, potential)
    # One frame transform per branch gives x, y and z together.
    x_rf_int, y_rf_int, z_rf_int = integrated_values(lambda sc: sc.galactic.cartesian.xyz.value) * 1000

    return (x_rf_int, y_rf_int, z_rf_int)


def center_orbit(coordinates_int, time, reference_frame_center, potential=None, vo=236., ro=8.122, zo=0.0208):
    """Subtract the reference-frame orbit from heliocentric ``(x, y, z)`` orbits.

    Takes the same arguments as :func:`get_center_orbit_coords`, plus the
    orbits to recentre; returns the frame-centred ``(x, y, z)`` in pc.
    """
    x_int, y_int, z_int = coordinates_int
    x_rf_int, y_rf_int, z_rf_int = get_center_orbit_coords(time, reference_frame_center, potential, vo, ro, zo)

    x_centered = x_int - x_rf_int
    y_centered = y_int - y_rf_int
    z_centered = z_int - z_rf_int

    return (x_centered, y_centered, z_centered)


def create_orbit(coordinates, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
    """Integrate phase-space positions through ``time``.

    Parameters
    ----------
    coordinates : sequence of arrays
        ``x, y, z`` (pc) and ``U, V, W`` (km/s), heliocentric Galactic.
    time, reference_frame_center, potential, vo, ro, zo
        As in :func:`get_center_orbit_coords`.

    Returns
    -------
    tuple
        Frame-centred ``(x, y, z)``, heliocentric ``(x, y, z)`` and
        Galactocentric ``(x, y, z)`` in pc, and frame-centred Galactocentric
        cylindrical ``(R, phi, z)``; each array is shaped (objects, times).
    """
    time = normalize_time_grid(time)
    if potential is None:
        potential = MWPotential2014

    x, y, z, U, V, W = coordinates
    sc = SkyCoord(
        u=x*u.pc, v=y*u.pc, w=z*u.pc, U=U*u.km/u.s, V=V*u.km/u.s, W=W*u.km/u.s,
        frame='galactic', representation_type='cartesian', differential_type='cartesian'
    )
    orbit = Orbit(vxvv=sc, ro=ro, vo=vo, zo=zo, solarmotion='schoenrich')

    integrated_values = _integrate_orbit_in_requested_order(orbit, time, potential)
    # One frame transform per branch and frame gives x, y and z together.
    helio_coords = tuple(integrated_values(lambda sc: sc.galactic.cartesian.xyz.value) * 1000)
    galactocentric_coords = tuple(integrated_values(lambda sc: sc.galactocentric.cartesian.xyz.value) * 1000)
    x_int_c, y_int_c, z_int_c = center_orbit(helio_coords, time, reference_frame_center, potential, vo, ro, zo)
    centered_coords = (x_int_c, y_int_c, z_int_c)

    # Galactocentric cylindrical coordinates of the frame-centred positions.
    sc_centered = SkyCoord(
        u=x_int_c*u.pc, v=y_int_c*u.pc, w=z_int_c*u.pc, frame='galactic', representation_type='cartesian')
    cylindrical = sc_centered.galactocentric.represent_as('cylindrical')
    gc_cylindcrical_centered_coords = (cylindrical.rho.value, cylindrical.phi.value, cylindrical.z.value)

    return (centered_coords, helio_coords, galactocentric_coords, gc_cylindcrical_centered_coords)


def _rotating_frame_xyz(x_gc, y_gc, z_gc, time_myr, r_sun=8.122, v_sun=236):
    """Galactocentric pc to the frame co-rotating with the Sun (Galactic centre at +R0 on x)."""
    w1 = (v_sun / r_sun) / 10
    t1 = time_myr * 0.01022
    r = np.sqrt(x_gc**2 + y_gc**2)
    theta = np.arctan2(x_gc, y_gc)
    x_rot = r_sun * 1000 - r * np.cos(theta - w1 * t1 + math.pi / 2)
    y_rot = r * np.sin(theta - w1 * t1 + math.pi / 2)
    return x_rot, y_rot, z_gc


def coordFIX_to_coordROT(df_gc, r_sun=8.122, v_sun=236):
    """Add co-rotating ``x_rot``, ``y_rot``, ``z_rot`` columns (pc) to ``df_gc``.

    ``df_gc`` needs ``time`` (Myr) and Galactocentric ``x_gc``, ``y_gc``,
    ``z_gc`` (pc); ``r_sun`` is in kpc and ``v_sun`` in km/s. The frame is
    modified in place and returned.
    """
    df_gc['x_rot'], df_gc['y_rot'], df_gc['z_rot'] = _rotating_frame_xyz(
        df_gc['x_gc'], df_gc['y_gc'], df_gc['z_gc'], df_gc['time'], r_sun=r_sun, v_sun=v_sun
    )
    return df_gc
