# Description: This file contains the functions to create the orbits of star clusters
from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit 
from galpy.potential import MWPotential2014

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
    orbit = Orbit(vxvv=sc, ro=8., vo=220., zo=0., solarmotion='schoenrich')
    orbit.integrate(time*u.yr, MWPotential2014)
    sc_int = orbit.SkyCoord(time*u.yr)
    x_int = sc_int.cartesian.x.value
    y_int = sc_int.cartesian.y.value
    z_int = sc_int.cartesian.z.value
    return (x_int, y_int, z_int) 
