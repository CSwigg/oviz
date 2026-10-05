"""Create time-dependent spiral-arm traces from galpy potentials."""

import warnings
from copy import deepcopy

import numpy as np
import pandas as pd
from galpy.potential import SpiralArmsPotential

from . import orbit_maker

warnings.filterwarnings("ignore")


class SpiralArmTrace:
    """Points along one galpy ``SpiralArmsPotential`` arm, integrated like clusters.

    The locus ``ln(R / r_ref) = -(phi - phi_ref) tan(alpha)`` is sampled
    between ``r_range`` and each point is then integrated as a test particle
    through the timeline. For the published arm models drawn by
    ``make_plot(spiral_arm_models=...)`` see :mod:`oviz.spiral_models`.

    Parameters
    ----------
    spiral_potential : galpy.potential.SpiralArmsPotential
        The arm's potential term.
    arm_name : str
        Label, e.g. ``"Carina-Sagittarius"``.
    r_range : tuple of float
        Galactocentric radii (kpc) the locus spans.
    n_points : int
        Number of locus points.
    color, opacity, line_width, visible
        Display style.

    After :meth:`integrate_orbits`, ``df_int`` holds the integrated table.
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

        self.integrated = False
        self.cluster_int_coords = None
        self.df_int = None

        self._extract_arm_parameters()
        self.df = self._generate_spiral_locus()
        self.coordinates = self.df[['x', 'y', 'z', 'U', 'V', 'W']].T.values

    def _extract_arm_parameters(self):
        """Copy the arm's parameters from the SpiralArmsPotential."""
        self.amp = self.spiral_potential._amp
        self.N = self.spiral_potential._N  # Number of arms
        self.alpha = self.spiral_potential._alpha  # Pitch angle
        self.r_ref = self.spiral_potential._r_ref  # Reference radius
        self.phi_ref = self.spiral_potential._phi_ref  # Reference azimuth
        self.Rs = self.spiral_potential._Rs  # Scale radius
        self.H = self.spiral_potential._H   # Scale height
        # Pattern speed; 0 for a static spiral.
        self.omega = getattr(self.spiral_potential, '_omega', 0.0)

    def _generate_spiral_locus(self):
        """The arm's locus as an Oviz phase-space table (pc, zero velocities)."""
        r_spiral = np.linspace(self.r_range[0], self.r_range[1], self.n_points)
        # ln(r / r_ref) = -(phi - phi_ref) tan(alpha), solved for phi.
        phi_spiral = self.phi_ref - np.log(r_spiral / self.r_ref) / np.tan(self.alpha)

        # Galactocentric kpc, in the plane.
        x_gc = r_spiral * np.cos(phi_spiral)
        y_gc = r_spiral * np.sin(phi_spiral)
        z_gc = np.zeros_like(r_spiral)

        # Heliocentric: shift by the solar radius (8.122 kpc) and height (27 pc).
        R_sun = 8.122
        x_pc = (x_gc - R_sun) * 1000
        y_pc = y_gc * 1000
        z_pc = (z_gc + 0.027) * 1000

        U = np.zeros_like(x_pc)
        V = np.zeros_like(x_pc)
        W = np.zeros_like(x_pc)

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
        """Integrate each locus point as a test particle through ``time`` (Myr).

        The pattern speed is not applied: points move with the Galactic
        potential. Arguments are those of :meth:`oviz.Trace.integrate_orbits`.
        """
        time = orbit_maker.normalize_time_grid(time)
        self.cluster_int_coords = orbit_maker.create_orbit(
            self.coordinates, time,
            reference_frame_center=reference_frame_center,
            potential=potential,
            vo=vo, ro=ro, zo=zo
        )
        self.df_int = self._create_integrated_dataframe(time, ro=ro, vo=vo)
        self.integrated = True

    def _create_integrated_dataframe(self, time, *, ro=8.122, vo=236.):
        """One row per locus point and time, like :meth:`oviz.Trace.create_integrated_dataframe`."""
        if self.cluster_int_coords is None:
            raise ValueError("Must integrate orbits before creating integrated DataFrame")

        xint, yint, zint = self.cluster_int_coords[0]
        xint_helio, yint_helio, zint_helio = self.cluster_int_coords[1]
        xint_gc, yint_gc, zint_gc = self.cluster_int_coords[2]
        rint_gc, phiint_gc, zint_gc_cyl = self.cluster_int_coords[3]

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
            'z_gc_cyl': zint_gc_cyl.flatten()
        })

        # Add spiral arm specific columns
        df_int['n_stars'] = np.ones(len(df_int))  # Dummy values
        df_int['age_myr'] = np.zeros(len(df_int))  # Spiral arms are "ageless"
        df_int['name'] = np.repeat(self.df['name'].values, len(time))
        df_int['time'] = np.tile(time, len(self.df))
        df_int['spiral_arm'] = self.arm_name

        df_int.reset_index(drop=True, inplace=True)
        df_int = orbit_maker.coordFIX_to_coordROT(df_int, r_sun=ro, v_sun=vo)
        return df_int

    def copy(self):
        """Return a deep copy of the arm trace."""
        return deepcopy(self)


def extract_spiral_arms_from_potential(potential, arm_names=None, **trace_kwargs):
    """One :class:`SpiralArmTrace` per ``SpiralArmsPotential`` term in ``potential``.

    ``arm_names`` names them in order (``Spiral_Arm_<n>`` otherwise); other
    keyword arguments go to :class:`SpiralArmTrace`.
    """
    if not isinstance(potential, list):
        potential = [potential]

    spiral_traces = []
    for pot_component in potential:
        if isinstance(pot_component, SpiralArmsPotential):
            index = len(spiral_traces)
            if arm_names is not None and index < len(arm_names):
                arm_name = arm_names[index]
            else:
                arm_name = f"Spiral_Arm_{index + 1}"
            spiral_traces.append(SpiralArmTrace(spiral_potential=pot_component, arm_name=arm_name, **trace_kwargs))
    return spiral_traces
