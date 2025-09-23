# Classes and functions for handling spiral arm visualization in oviz
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit
from galpy.potential import SpiralArmsPotential
from copy import deepcopy
import warnings
warnings.filterwarnings("ignore")

from . import orbit_maker

class SpiralArmTrace:
    """
    Class for creating animated spiral arm traces from galpy SpiralArmsPotential objects.
    
    This class extracts the spiral arm locus from SpiralArmsPotential and creates
    a trace that can be animated in oviz visualizations, showing how the spiral 
    pattern moves over time.
    
    Parameters
    ----------
    spiral_potential : galpy.potential.SpiralArmsPotential
        A galpy SpiralArmsPotential object defining the spiral arm
    arm_name : str
        Name for this spiral arm (e.g., "Carina-Sagittarius")
    r_range : tuple, optional
        Galactocentric radius range (min_r, max_r) in kpc. Default (4.0, 12.0)
    n_points : int, optional
        Number of points along the spiral arm. Default 200
    color : str, optional
        Color for the spiral arm trace. Default 'gray'
    opacity : float, optional
        Opacity of the spiral arm trace. Default 0.6
    line_width : float, optional
        Width of the spiral arm line. Default 2.0
    
    Attributes
    ----------
    arm_name : str
        Name of the spiral arm
    spiral_potential : galpy.potential.SpiralArmsPotential
        The spiral arm potential object
    df : pandas.DataFrame
        DataFrame containing spiral arm coordinates in oviz format
    coordinates : numpy.ndarray
        Array of coordinates (x, y, z, U, V, W) of spiral arm points
    color : str
        Color of the spiral arm
    opacity : float
        Opacity of the spiral arm
    line_width : float
        Width of the spiral arm line
    integrated : bool
        Whether the spiral arm has been orbit-integrated
    cluster_int_coords : tuple or None
        Integrated coordinates if orbit integration has been performed
    df_int : pandas.DataFrame or None
        Integrated DataFrame if orbit integration has been performed
    """
    
    def __init__(self, spiral_potential, arm_name, r_range=(4.0, 12.0), n_points=200, 
                 color='gray', opacity=0.6, line_width=2.0, visible=True):
        
        if not isinstance(spiral_potential, SpiralArmsPotential):
            raise ValueError("spiral_potential must be a galpy SpiralArmsPotential object")
            
        self.spiral_potential = spiral_potential
        self.arm_name = arm_name
        self.r_range = r_range
        self.n_points = n_points
        self.color = color
        self.opacity = opacity
        self.line_width = line_width
        self.visible = visible
        
        # Initialize as not integrated
        self.integrated = False
        self.cluster_int_coords = None
        self.df_int = None
        
        # Extract spiral arm parameters
        self._extract_arm_parameters()
        
        # Generate spiral arm locus
        self.df = self._generate_spiral_locus()
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values
        
    def _extract_arm_parameters(self):
        """Extract parameters from the SpiralArmsPotential object."""
        # Get the spiral arm parameters
        self.amp = self.spiral_potential._amp
        self.N = self.spiral_potential._N  # Number of arms
        self.alpha = self.spiral_potential._alpha  # Pitch angle
        self.r_ref = self.spiral_potential._r_ref  # Reference radius
        self.phi_ref = self.spiral_potential._phi_ref  # Reference azimuth
        self.Rs = self.spiral_potential._Rs  # Scale radius
        self.H = self.spiral_potential._H   # Scale height
        
        # Pattern speed (if available)
        if hasattr(self.spiral_potential, '_omega'):
            self.omega = self.spiral_potential._omega
        else:
            self.omega = 0.0  # Static spiral
            
    def _generate_spiral_locus(self):
        """
        Generate the spiral arm locus coordinates.
        
        Returns
        -------
        pandas.DataFrame
            DataFrame with spiral arm coordinates in oviz format
        """
        # Generate radial range
        r_spiral = np.linspace(self.r_range[0], self.r_range[1], self.n_points)
        
        # Calculate azimuthal angles from spiral equation
        # ln(r/r_ref) = -(phi - phi_ref) * tan(alpha)
        # Solving for phi: phi = phi_ref - ln(r/r_ref) / tan(alpha)
        phi_spiral = self.phi_ref - np.log(r_spiral / self.r_ref) / np.tan(self.alpha)
        
        # Convert to Cartesian galactocentric coordinates
        x_gc = r_spiral * np.cos(phi_spiral)  # kpc
        y_gc = r_spiral * np.sin(phi_spiral)  # kpc
        z_gc = np.zeros_like(r_spiral)  # Assume spiral is in galactic plane
        
        # Convert from galactocentric to heliocentric coordinates
        # (oviz uses heliocentric coordinates)
        R_sun = 8.122  # kpc, solar galactocentric distance
        
        # Transform: heliocentric = galactocentric + solar position
        x_helio = x_gc - R_sun  # kpc
        y_helio = y_gc          # kpc
        z_helio = z_gc + 0.027  # kpc, add solar height above plane
        
        # Convert to parsecs (oviz units)
        x_pc = x_helio * 1000  # pc
        y_pc = y_helio * 1000  # pc
        z_pc = z_helio * 1000  # pc
        
        # Assign zero velocities initially (will be updated during orbit integration)
        U = np.zeros_like(x_pc)  # km/s
        V = np.zeros_like(x_pc)  # km/s  
        W = np.zeros_like(x_pc)  # km/s
        
        # Create DataFrame in oviz format
        df_spiral = pd.DataFrame({
            'x': x_pc,
            'y': y_pc, 
            'z': z_pc,
            'U': U,
            'V': V,
            'W': W,
            'name': [f'{self.arm_name}_pt_{i:03d}' for i in range(len(x_pc))],
            'age_myr': np.full(len(x_pc), 0.0),  # Spiral arms are "ageless"
            'n_stars': np.ones(len(x_pc))  # Dummy value for compatibility
        })
        
        return df_spiral
        
    def integrate_orbits(self, time, reference_frame_center=None, potential=None, vo=236., ro=8.122, zo=0.0208):
        """
        Integrate spiral arm motion over time.
        
        For spiral arms, this accounts for:
        1. Pattern rotation (if omega != 0)
        2. Orbital motion of individual spiral arm "particles"
        
        Parameters
        ----------
        time : array
            Array of time points in Myr
        reference_frame_center : tuple, optional
            Center of the reference frame
        potential : galpy potential, optional
            Galactic potential to use for orbit integration
        vo : float, optional
            Circular velocity in km/s
        ro : float, optional
            Solar radius in kpc
        zo : float, optional
            Solar height in kpc
        """
        # Calculate pattern rotation effect
        omega_pattern = self.omega  # Pattern speed in km/s/kpc
        
        # For each time step, rotate the spiral pattern
        time_array = np.array(time)
        
        # Initialize storage for time-evolved coordinates
        all_coords = []
        
        for t in time_array:
            # Calculate rotation angle for this time
            # Convert time from Myr to years, then to appropriate units
            t_years = t * 1e6  # Myr to years
            t_seconds = t_years * 3.15576e7  # years to seconds
            
            # Pattern rotation: Δφ = ω * Δt
            if omega_pattern != 0:
                # omega is in km/s/kpc, need to convert to rad/s
                omega_rad_per_s = omega_pattern / ro  # Convert to 1/s
                delta_phi = omega_rad_per_s * t_seconds  # radians
            else:
                delta_phi = 0
                
            # Rotate spiral pattern
            coords_t = self._rotate_spiral_pattern(delta_phi)
            all_coords.append(coords_t)
            
        # Now perform standard orbit integration
        # (This handles motion of spiral arm material through the galaxy)
        self.cluster_int_coords = orbit_maker.create_orbit(
            self.coordinates, time,
            reference_frame_center=reference_frame_center,
            potential=potential,
            vo=vo, ro=ro, zo=zo
        )
        
        # Create integrated DataFrame
        self.df_int = self._create_integrated_dataframe(time)
        self.integrated = True
        
    def _rotate_spiral_pattern(self, delta_phi):
        """
        Rotate the spiral pattern by angle delta_phi.
        
        Parameters
        ----------
        delta_phi : float
            Rotation angle in radians
            
        Returns
        -------
        numpy.ndarray
            Rotated coordinates
        """
        # Get current coordinates
        x, y, z, U, V, W = self.coordinates
        
        # Convert to polar coordinates in the galactic plane
        r = np.sqrt((x/1000 + 8.122)**2 + (y/1000)**2)  # Convert pc to kpc, add solar position
        phi = np.arctan2(y/1000, x/1000 + 8.122)
        
        # Apply pattern rotation
        phi_new = phi + delta_phi
        
        # Convert back to Cartesian
        x_gc_new = r * np.cos(phi_new)
        y_gc_new = r * np.sin(phi_new)
        
        # Convert back to heliocentric coordinates
        x_new = (x_gc_new - 8.122) * 1000  # kpc to pc, subtract solar position
        y_new = y_gc_new * 1000
        z_new = z  # No change in z
        
        return np.array([x_new, y_new, z_new, U, V, W])
        
    def _create_integrated_dataframe(self, time):
        """
        Create integrated DataFrame for spiral arm motion.
        
        Parameters
        ----------
        time : array
            Array of time points
            
        Returns
        -------
        pandas.DataFrame
            Integrated DataFrame
        """
        if self.cluster_int_coords is None:
            raise ValueError("Must integrate orbits before creating integrated DataFrame")
            
        xint, yint, zint = self.cluster_int_coords[0]
        xint_helio, yint_helio, zint_helio = self.cluster_int_coords[1]
        xint_gc, yint_gc, zint_gc = self.cluster_int_coords[2]
        rint_gc, phiint_gc, zint_gc = self.cluster_int_coords[3]

        df_int = pd.DataFrame({
            'x': xint.flatten(), 
            'y': yint.flatten(), 
            'z': zint.flatten(),
            'x_helio': xint_helio.flatten(),
            'y_helio': yint_helio.flatten(),
            'z_helio': zint_helio.flatten(),
            'x_gc': xint_gc.flatten(),
            'y_gc': yint_gc.flatten(),
            'z_gc': zint_gc.flatten(),
            'r_gc': rint_gc.flatten(),
            'phi_gc': np.rad2deg(phiint_gc.flatten()),
            'z_gc_cyl': zint_gc.flatten()
        })

        # Add spiral arm specific columns
        df_int['n_stars'] = np.ones(len(df_int))  # Dummy values
        df_int['age_myr'] = np.zeros(len(df_int))  # Spiral arms are "ageless"
        df_int['name'] = np.repeat(self.df['name'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))
        df_int['spiral_arm'] = self.arm_name

        df_int.reset_index(drop=True, inplace=True)
        df_int = orbit_maker.coordFIX_to_coordROT(df_int)
        
        return df_int
        
    def copy(self):
        """
        Create a copy of the SpiralArmTrace.
        
        Returns
        -------
        SpiralArmTrace
            A copy of this spiral arm trace
        """
        return deepcopy(self)


def extract_spiral_arms_from_potential(potential, arm_names=None, **trace_kwargs):
    """
    Extract SpiralArmTrace objects from a galpy potential containing spiral arms.
    
    Parameters
    ----------
    potential : list or galpy potential
        Galpy potential object or list of potential objects
    arm_names : list of str, optional
        Names for the spiral arms. If None, will use generic names
    **trace_kwargs
        Additional keyword arguments passed to SpiralArmTrace constructor
        
    Returns
    -------
    list of SpiralArmTrace
        List of SpiralArmTrace objects, one for each spiral arm component
    """
    # Ensure potential is a list
    if not isinstance(potential, list):
        potential = [potential]
        
    spiral_traces = []
    spiral_count = 0
    
    for pot_component in potential:
        if isinstance(pot_component, SpiralArmsPotential):
            # Determine arm name
            if arm_names is not None and spiral_count < len(arm_names):
                arm_name = arm_names[spiral_count]
            else:
                arm_name = f"Spiral_Arm_{spiral_count + 1}"
                
            # Create SpiralArmTrace
            spiral_trace = SpiralArmTrace(
                spiral_potential=pot_component,
                arm_name=arm_name,
                **trace_kwargs
            )
            
            spiral_traces.append(spiral_trace)
            spiral_count += 1
            
    return spiral_traces