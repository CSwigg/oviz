# Description: This file contains the functions to create the orbits of star clusters
from astropy.coordinates import SkyCoord
import astropy.units as u
import math
import numpy as np
from galpy.orbit import Orbit 
import warnings 
warnings.filterwarnings("ignore")
#from galpy.potential import MWPotential2014, SpiralArmsPotential
from galpy.potential import MWPotential2014, DehnenBarPotential, SpiralArmsPotential

def get_center_orbit_coords(time, reference_frame_center):
    """
    Get the coordinates of the center orbit.

    Parameters:
    - time (array): Array of time points.
    - reference_frame_center (tuple): Center of the reference frame.

    Returns:
    tuple: Coordinates of the center orbit (x, y, z).
    """
    if reference_frame_center is None:
        rf_coords = [0, 0, 0, -11.1, -12.24, -7.25]
    else:
        rf_coords = reference_frame_center 

    rf_sc = SkyCoord(
        u=rf_coords[0]*u.pc, v=rf_coords[1]*u.pc, w=rf_coords[2]*u.pc,
        U=rf_coords[3]*u.km/u.s, V=rf_coords[4]*u.km/u.s, W=rf_coords[5]*u.km/u.s,
        frame='galactic', representation_type='cartesian', differential_type='cartesian'
    )
    rf_orbit = Orbit(vxvv=rf_sc, solarmotion='schoenrich', ro=8.122, vo=236., zo=0.0208)
    rf_orbit.integrate(time*u.Myr, MWPotential2014)
    x_rf_int = rf_orbit.SkyCoord(time*u.Myr).galactic.cartesian.x.value * 1000
    y_rf_int = rf_orbit.SkyCoord(time*u.Myr).galactic.cartesian.y.value * 1000
    z_rf_int = rf_orbit.SkyCoord(time*u.Myr).galactic.cartesian.z.value * 1000
    return (x_rf_int, y_rf_int, z_rf_int)

def center_orbit(coordinates_int, time, reference_frame_center):
    """
    Centers the given coordinates with respect to the Local Standard of Rest (LSR) at a specific time.

    Parameters:
    - coordinates_int (tuple): The initial coordinates (x, y, z) in the galactic frame.
    - time (array): Array of time points.
    - reference_frame_center (tuple): Center of the reference frame.

    Returns:
    tuple: The centered coordinates (x_centered, y_centered, z_centered) in the galactic frame.
    """
    x_int, y_int, z_int = coordinates_int
    x_rf_int, y_rf_int, z_rf_int = get_center_orbit_coords(time, reference_frame_center)

    x_centered = x_int - x_rf_int
    y_centered = y_int - y_rf_int
    z_centered = z_int - z_rf_int

    return (x_centered, y_centered, z_centered)

def create_orbit(coordinates, time, reference_frame_center=None):
    """
    Create orbits for star cluster(s).

    Parameters:
    - coordinates (tuple): A tuple containing the x, y, z coordinates and U, V, W velocities of the star cluster.
    - time (array): Array of time points.
    - reference_frame_center (tuple, optional): Center of the reference frame.

    Returns:
    tuple: Centered, heliocentric, and galactocentric coordinates.
    """
    x, y, z, U, V, W = coordinates
    sc = SkyCoord(
        u=x*u.pc, v=y*u.pc, w=z*u.pc, U=U*u.km/u.s, V=V*u.km/u.s, W=W*u.km/u.s,
        frame='galactic', representation_type='cartesian', differential_type='cartesian'
    )
    orbit = Orbit(vxvv=sc, ro=8.122, vo=236, zo=0.0208, solarmotion='schoenrich')
    #pot = MWPotential2014 + SpiralArmsPotential(omega=13.5*u.km/u.s/u.kpc, amp=0.15)


    # Define the bar potential
    # bar_potential = DehnenBarPotential(
    #     omegab=37 * u.km/u.s/u.kpc,  # Pattern speed
    #     rb=12 * u.kpc,               # Bar radius
    #     barphi=28 * u.deg,           # Bar orientation angle
    #     alpha=0.01                   # Relative bar strength
    # )

    # Define the m=2 spiral arms potential
    spiral2_potential = SpiralArmsPotential(
        amp=0.249,                    # Amplitude
        N=2,                         # Number of spiral arms
        alpha=8.1 * u.deg,           # Pitch angle
        r_ref=6.6 * u.kpc,           # Reference radius
        phi_ref=47.8 * u.deg,        # Reference angle
        Rs=26.4 * u.kpc,             # Radial scale length
        H=0.13 * u.kpc,              # Scale height
        omega=13.1 * u.km/u.s/u.kpc, # Pattern speed
        Cs=[1]                       # Coefficients for harmonic terms
    )

    # Define the m=3 spiral arms potential
    spiral3_potential = SpiralArmsPotential(
        amp=0.093,                    # Amplitude
        N=3,                         # Number of spiral arms
        alpha=13.7 * u.deg,          # Pitch angle
        r_ref=8.0 * u.kpc,           # Reference radius
        phi_ref=81.7 * u.deg,        # Reference angle
        Rs=19.6 * u.kpc,             # Radial scale length
        H=0.13 * u.kpc,              # Scale height
        omega=16.4 * u.km/u.s/u.kpc, # Pattern speed
        Cs=[1]                       # Coefficients for harmonic terms
    )

    # Combine with MWPotential2014
    #pot = MWPotential2014 + [spiral2_potential, spiral3_potential]
    pot = MWPotential2014
    orbit.integrate(time*u.Myr, pot)
    sc_int = orbit.SkyCoord(time*u.Myr)

    x_helio_int = sc_int.galactic.cartesian.x.value * 1000
    y_helio_int = sc_int.galactic.cartesian.y.value * 1000
    z_helio_int = sc_int.galactic.cartesian.z.value * 1000
    helio_coords = (x_helio_int, y_helio_int, z_helio_int)

    x_gc_int = sc_int.galactocentric.cartesian.x.value * 1000
    y_gc_int = sc_int.galactocentric.cartesian.y.value * 1000
    z_gc_int = sc_int.galactocentric.cartesian.z.value * 1000
    galactocentric_coords = (x_gc_int, y_gc_int, z_gc_int)

    x_int_c, y_int_c, z_int_c = center_orbit((x_helio_int, y_helio_int, z_helio_int), time, reference_frame_center)
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