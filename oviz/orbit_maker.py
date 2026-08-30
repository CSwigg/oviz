"""Integrate Galactic orbits and convert them to Oviz coordinate frames."""

from astropy.coordinates import SkyCoord
import astropy.units as u
from copy import deepcopy
import math
import numpy as np
from galpy.orbit import Orbit 
import warnings 
warnings.filterwarnings("ignore")
#from galpy.potential import MWPotential2014, SpiralArmsPotential
from galpy.potential import MWPotential2014, DehnenBarPotential, SpiralArmsPotential


def normalize_time_grid(time, *, require_zero=True):
    """Return a validated one-dimensional timeline in caller-supplied order."""

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
    """Integrate outward from zero and return SkyCoords in requested order."""

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

def get_center_orbit_coords(time, reference_frame_center, potential=None, vo=236., ro=8.122, zo=0.0208):
    """
    Get the coordinates of the center orbit.

    Parameters:
    - time (array): Array of time points.
    - reference_frame_center (tuple): Center of the reference frame.
    - potential (galpy potential, optional): Galactic potential to use. Defaults to MWPotential2014.
    - vo (float, optional): Circular velocity at the solar radius in km/s. Defaults to 236.
    - ro (float, optional): Solar radius in kpc. Defaults to 8.122.
    - zo (float, optional): Solar height above the galactic plane in kpc. Defaults to 0.0208.

    Returns:
    tuple: Coordinates of the center orbit (x, y, z).
    """
    time = normalize_time_grid(time)
    if potential is None:
        potential = MWPotential2014
        
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
    x_rf_int = integrated_values(lambda sc: sc.galactic.cartesian.x.value) * 1000
    y_rf_int = integrated_values(lambda sc: sc.galactic.cartesian.y.value) * 1000
    z_rf_int = integrated_values(lambda sc: sc.galactic.cartesian.z.value) * 1000

    return (x_rf_int, y_rf_int, z_rf_int)

def center_orbit(coordinates_int, time, reference_frame_center, potential=None, vo=236., ro=8.122, zo=0.0208):
    """
    Centers the given coordinates with respect to the Local Standard of Rest (LSR) at a specific time.

    Parameters:
    - coordinates_int (tuple): The initial coordinates (x, y, z) in the galactic frame.
    - time (array): Array of time points.
    - reference_frame_center (tuple): Center of the reference frame.
    - potential (galpy potential, optional): Galactic potential to use. Defaults to MWPotential2014.
    - vo (float, optional): Circular velocity at the solar radius in km/s. Defaults to 236.
    - ro (float, optional): Solar radius in kpc. Defaults to 8.122.
    - zo (float, optional): Solar height above the galactic plane in kpc. Defaults to 0.0208.

    Returns:
    tuple: The centered coordinates (x_centered, y_centered, z_centered) in the galactic frame.
    """
    x_int, y_int, z_int = coordinates_int
    x_rf_int, y_rf_int, z_rf_int = get_center_orbit_coords(time, reference_frame_center, potential, vo, ro, zo)

    x_centered = x_int - x_rf_int
    y_centered = y_int - y_rf_int
    z_centered = z_int - z_rf_int

    return (x_centered, y_centered, z_centered)

def create_orbit(coordinates, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
    """
    Create orbits for star cluster(s).

    Parameters:
    - coordinates (tuple): A tuple containing the x, y, z coordinates and U, V, W velocities of the star cluster.
    - time (array): Array of time points.
    - reference_frame_center (tuple, optional): Center of the reference frame.
    - potential (galpy potential, optional): Galactic potential to use. Defaults to MWPotential2014.
    - vo (float, optional): Circular velocity at the solar radius in km/s. Defaults to 236.
    - ro (float, optional): Solar radius in kpc. Defaults to 8.122.
    - zo (float, optional): Solar height above the galactic plane in kpc. Defaults to 0.0208.

    Returns:
    tuple: Centered, heliocentric, and galactocentric coordinates.
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
    x_helio_int = integrated_values(lambda sc: sc.galactic.cartesian.x.value) * 1000
    y_helio_int = integrated_values(lambda sc: sc.galactic.cartesian.y.value) * 1000
    z_helio_int = integrated_values(lambda sc: sc.galactic.cartesian.z.value) * 1000
    x_gc_int = integrated_values(lambda sc: sc.galactocentric.cartesian.x.value) * 1000
    y_gc_int = integrated_values(lambda sc: sc.galactocentric.cartesian.y.value) * 1000
    z_gc_int = integrated_values(lambda sc: sc.galactocentric.cartesian.z.value) * 1000
    
    helio_coords = (x_helio_int, y_helio_int, z_helio_int)
    galactocentric_coords = (x_gc_int, y_gc_int, z_gc_int)
    x_int_c, y_int_c, z_int_c = center_orbit(helio_coords, time, reference_frame_center, potential, vo, ro, zo)
    centered_coords = (x_int_c, y_int_c, z_int_c)

    # Construct galactocentric cylindrical coordinates, centered
    sc_centered = SkyCoord(
        u=x_int_c*u.pc, v=y_int_c*u.pc, w=z_int_c*u.pc, frame='galactic', representation_type='cartesian')
    R_gc_int_c, phi_gc_int_c, z_gc_int_c = sc_centered.galactocentric.represent_as('cylindrical').rho.value, sc_centered.galactocentric.represent_as('cylindrical').phi.value, sc_centered.galactocentric.represent_as('cylindrical').z.value
    gc_cylindcrical_centered_coords = (R_gc_int_c, phi_gc_int_c, z_gc_int_c)

    return (centered_coords, helio_coords, galactocentric_coords, gc_cylindcrical_centered_coords)



def coordFIX_to_coordROT(df_gc, r_sun=8.122, v_sun=236):
    """
    Convert fixed coordinates to rotating coordinates.

    Parameters:
    - df_gc (pd.DataFrame): Data frame containing galactocentric coordinates.
    - r_sun (float, optional): Distance from the Sun to the galactic center in kpc. Defaults to 8.122.
    - v_sun (float, optional): Circular velocity of the Sun in km/s. Defaults to 236.

    Returns:
    pd.DataFrame: Data frame with rotating coordinates.
    """
    w0 = v_sun / r_sun
    r_sun = r_sun * 1000  # in pc!
    w1 = w0 / 10
    t1 = df_gc['time'] * 0.01022
    r = np.sqrt(df_gc['x_gc']**2 + df_gc['y_gc']**2)
    theta = np.arctan2(df_gc['x_gc'], df_gc['y_gc'])
    x_gc_rot = r_sun - r * np.cos(theta - w1 * t1 + math.pi / 2)
    y_gc_rot = r * np.sin(theta - w1 * t1 + math.pi / 2)
    z_gc_rot = df_gc['z_gc']

    df_gc['x_rot'] = x_gc_rot
    df_gc['y_rot'] = y_gc_rot
    df_gc['z_rot'] = z_gc_rot

    return df_gc
