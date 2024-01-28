# Description: This file contains the functions to create the orbits of star clusters
from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit 
from galpy.potential import MWPotential2014

def center_orbit(coordinates_int, time):
    """
    Centers the given coordinates with respect to the Local Standard of Rest (LSR) at a specific time.

    Parameters:
    - coordinates_int (tuple): The initial coordinates (x, y, z) in the galactic frame.
    - time (float): The time at which to calculate the LSR coordinates.

    Returns:
    - tuple: The centered coordinates (x_centered, y_centered, z_centered) in the galactic frame.
    """
     
    x_int, y_int, z_int = coordinates_int

    # schoenrich lsr
    u_lsr = -11.1
    v_lsr = -12.24
    w_lsr = -7.25

    lsr_velocities = [u_lsr, v_lsr, w_lsr]
    solarmotion = [u_lsr, -1*v_lsr, -1*w_lsr]

    sc_orbit = SkyCoord(u = 0.*u.pc, v = 0*u.pc, w = 0.*u.pc, U = lsr_velocities[0]*u.km/u.s, V = lsr_velocities[1]*u.km/u.s, W = lsr_velocities[2]*u.km/u.s, frame = 'galactic', representation_type = 'cartesian', differential_type = 'cartesian')
    lsr_orbit = Orbit(vxvv = sc_orbit, solarmotion=solarmotion, ro = 8., vo = 220., zo = 0.0208)
    lsr_orbit.integrate(time*u.Myr, MWPotential2014)

    x_lsr_int = lsr_orbit.SkyCoord(time*u.Myr).galactic.cartesian.x.value*1000
    y_lsr_int = lsr_orbit.SkyCoord(time*u.Myr).galactic.cartesian.y.value*1000
    z_lsr_int = lsr_orbit.SkyCoord(time*u.Myr).galactic.cartesian.z.value*1000

    x_centered = x_int - x_lsr_int
    y_centered = y_int - y_lsr_int
    z_centered = z_int - z_lsr_int

    return (x_centered, y_centered, z_centered)

def create_orbit(coordinates, time):
    """
    Create orbits for star cluster(s).

    Parameters:
    coordinates (tuple): A tuple containing the x, y, z coordinates and U, V, W velocities of the star cluster.
    time (float): The time at which the orbit is created. Units should be in years.

    Returns:
    Orbit: An Orbit object representing the orbit of the star cluster.
    """
    x, y, z, U, V, W = coordinates
    sc = SkyCoord(u=x*u.pc, v=y*u.pc, w=z*u.pc, U=U*u.km/u.s, V=V*u.km/u.s, W=W*u.km/u.s, frame='galactic', representation_type='cartesian', differential_type='cartesian')
    orbit = Orbit(vxvv=sc, ro=8., vo=220., zo=0.0208, solarmotion='schoenrich')
    orbit.integrate(time*u.Myr, MWPotential2014)
    sc_int = orbit.SkyCoord(time*u.Myr)
    x_int = sc_int.galactic.cartesian.x.value*1000
    y_int = sc_int.galactic.cartesian.y.value*1000
    z_int = sc_int.galactic.cartesian.z.value*1000

    x_int_c, y_int_c, z_int_c = center_orbit((x_int, y_int, z_int), time)
    return (x_int_c, y_int_c, z_int_c) 
