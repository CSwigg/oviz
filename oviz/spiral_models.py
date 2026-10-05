"""Published Milky Way spiral arms, and the Khalil et al. (2025) spiral potential.

Each arm model is drawn as one trace holding all of its arms, turning at its
published pattern speeds as the figure's time moves (pass the models to
``make_plot(spiral_arm_models=...)``; the traces start hidden):

* :data:`KHALIL2025_ARMS`: the potential troughs of the two-armed (m = 2) and
  three-armed (m = 3) spiral modes that Khalil et al. (2025, A&A,
  arXiv:2411.12800) fitted to Gaia DR3 disk kinematics, with the values of
  their released SPIBACK code. The m = 2 mode follows the Crux-Scutum, Local
  and Outer arm segments, the m = 3 mode the Carina-Sagittarius and Perseus
  segments; they turn at 13.1 and 16.4 km/s/kpc.
* :data:`CASTRO_GINARD2021_ARMS`: the Perseus, Local, Sagittarius and Scutum
  segments that Castro-Ginard et al. (2021, A&A 652, A162) fitted to young
  open clusters and masers (their Table 1), each turning at the pattern speed
  they measured from cluster birthplaces (Table 2, 10-50 Myr).

:func:`khalil2025_potential` is galpy's MWPotential2014 plus the two Khalil
spiral modes as ``SpiralArmsPotential`` terms (the Cox & Gomez form SPIBACK
uses). In the plane it reproduces SPIBACK's spiral potential exactly, with
the same amplitudes, pitches, phases and pattern speeds; it adds galpy's
vertical profile and leaves out the paper's own axisymmetric background, its
bar and the modes' hard radial cutoffs.

Azimuths follow galpy: the Sun at phi = 0 and R = ro, phi growing in the
direction of Galactic rotation, time in Myr (negative in the past).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import astropy.units as u
import numpy as np
from astropy import constants as const

# Pattern speed to angular rate: 1 km/s/kpc is 1.0227e-3 rad/Myr.
KMS_KPC_TO_RAD_MYR = float((u.km / u.s / u.kpc).to(1 / u.Myr))

KHALIL2025_REFERENCE = "Khalil et al. 2025, A&A (arXiv:2411.12800); SPIBACK f270b32"
# The paper's solar radius: the spiral phases and amplitudes refer to it.
KHALIL2025_SOLAR_RADIUS_KPC = 8.275
KHALIL2025_SPIRAL_HEIGHT_KPC = 0.130
# SPIBACK's spirals do not decay with radius; galpy's term does, so its scale
# length is long enough to stay flat within 1.4% over 4-18 kpc.
KHALIL2025_SPIRAL_SCALE_LENGTH_KPC = 1000.0
# SPIBACK's released fiducial values: amplitude K (km/s)^2, pattern speed
# (km/s/kpc), pitch and present-day phase at the Sun (rad).
KHALIL2025_SPIRAL_MODES = (
    {"m": 2, "amplitude_kms2": 2818.0911260799921, "pattern_speed": 13.070958619442987,
     "pitch_rad": 0.14147012374228302, "phase_rad": 3.9764854092486122},
    {"m": 3, "amplitude_kms2": 1047.2810628832635, "pattern_speed": 16.406153217300506,
     "pitch_rad": 0.23862412058190050, "phase_rad": 1.4261355010802388},
)

CASTRO_GINARD2021_REFERENCE = "Castro-Ginard et al. 2021, A&A 652, A162 (Tables 1 and 2)"
# The solar radius their arm fits assumed (Gravity Collaboration 2019).
CASTRO_GINARD2021_SOLAR_RADIUS_KPC = 8.178
# Eq. 1, ln(R/R_ref) = -(theta - theta_ref) tan(pitch), with theta = 0 on the
# Sun-Galactic centre line, growing in the direction of rotation:
# (arm, theta_ref deg, theta range deg, R_ref kpc, pitch deg, pattern speed km/s/kpc).
CASTRO_GINARD2021_ARM_TABLE = (
    ("Perseus", -13.0, (-20.9, 88.2), 10.88, 9.8, 17.82),
    ("Local", -2.3, (-26.9, 26.6), 8.69, 8.9, 33.76),
    ("Sagittarius", 3.5, (-39.3, 67.7), 7.10, 10.6, 26.10),
    ("Scutum", -4.8, (-32.7, 100.9), 6.02, 14.9, 49.81),
)


def khalil2025_spiral_potentials(*, ro: float = 8.122, vo: float = 236.0) -> list:
    """The two Khalil et al. (2025) spiral modes as galpy ``SpiralArmsPotential`` terms.

    SPIBACK's mode is ``-K / (R0 k D) cos[m (phi - phase - Omega t +
    ln(R/R0) / tan p)]`` with ``k = m / (R sin p)`` and the Cox & Gomez
    ``D(kH)``. galpy's term is ``-G amp H / (k D) cos(gamma)`` for a density
    amplitude ``amp``, so ``amp = K / (G H R0)``. galpy turns its frame
    left-handed by negating N and alpha internally: a positive alpha gives
    the trailing arms SPIBACK has (a negative one flips the troughs to crests).
    """
    from galpy.potential import SpiralArmsPotential

    out = []
    for mode in KHALIL2025_SPIRAL_MODES:
        pitch = mode["pitch_rad"]
        amp = (
            mode["amplitude_kms2"] * (u.km / u.s) ** 2
            / (const.G * KHALIL2025_SPIRAL_HEIGHT_KPC * u.kpc * KHALIL2025_SOLAR_RADIUS_KPC * u.kpc)
        ).to(u.Msun / u.pc**3)
        # Referred to galpy's r_ref = ro, the phase keeps the troughs where
        # the paper puts them in Galactocentric radius.
        phi_ref = mode["phase_rad"] - np.log(ro / KHALIL2025_SOLAR_RADIUS_KPC) / np.tan(pitch)
        out.append(SpiralArmsPotential(
            amp=amp, amp_units="density", N=mode["m"], alpha=pitch * u.rad,
            r_ref=ro * u.kpc, phi_ref=phi_ref * u.rad,
            Rs=KHALIL2025_SPIRAL_SCALE_LENGTH_KPC * u.kpc, H=KHALIL2025_SPIRAL_HEIGHT_KPC * u.kpc,
            omega=mode["pattern_speed"] * u.km / u.s / u.kpc, Cs=[1], ro=ro, vo=vo,
        ))
    return out


def khalil2025_potential(*, ro: float = 8.122, vo: float = 236.0) -> list:
    """MWPotential2014 plus the Khalil et al. (2025) m = 2 and m = 3 spirals.

    The axisymmetric part is Oviz's usual MWPotential2014 (scaled by ``ro``
    and ``vo`` like every Oviz orbit); the reference frame's LSR keeps its
    circular orbit in it (see :func:`oviz.orbit_maker.reference_potential`).
    """
    from galpy.potential import MWPotential2014

    return [*MWPotential2014, *khalil2025_spiral_potentials(ro=ro, vo=vo)]


@dataclass(frozen=True)
class SpiralArmModel:
    """A published set of spiral arms, drawn as one line trace."""

    name: str
    reference: str
    color: str
    dash: str = "solid"
    width: float = 2.0
    opacity: float = 0.9

    def curves(self, time_myr: float, *, ro: float = 8.122) -> list[tuple[str, np.ndarray, np.ndarray]]:
        """``(arm, R kpc, phi rad)`` per arm at ``time_myr``, in galpy's Galactocentric frame."""
        raise NotImplementedError

    def pattern_speeds(self) -> dict[str, float]:
        """Pattern speed (km/s/kpc) of each arm."""
        raise NotImplementedError

    def provenance(self) -> dict:
        return {"name": self.name, "reference": self.reference, "pattern_speeds_kms_kpc": self.pattern_speeds()}


@dataclass(frozen=True)
class KhalilSpiralArms(SpiralArmModel):
    """The troughs of the Khalil et al. (2025) spiral potential: 2 + 3 arms.

    A trough is where the mode's potential is deepest at each radius,
    ``phi = phase - ln(R/R0) / tan p + Omega t + 2 pi k / m``. The arms are
    drawn over ``radius_range_kpc`` (the potential has no cutoffs) every
    ``step_deg`` of azimuth along them.
    """

    radius_range_kpc: tuple[float, float] = (4.0, 18.0)
    step_deg: float = 1.5

    def curves(self, time_myr, *, ro=8.122):
        r0, r1 = self.radius_range_kpc
        out = []
        for mode in KHALIL2025_SPIRAL_MODES:
            pitch, m = mode["pitch_rad"], mode["m"]
            n = int(np.ceil(np.log(r1 / r0) / np.tan(pitch) / np.radians(self.step_deg))) + 1
            radius = np.geomspace(r0, r1, n)
            base = (mode["phase_rad"] - np.log(radius / KHALIL2025_SOLAR_RADIUS_KPC) / np.tan(pitch)
                    + mode["pattern_speed"] * KMS_KPC_TO_RAD_MYR * float(time_myr))
            for k in range(m):
                out.append((f"m={m} arm {k + 1}", radius, base + 2.0 * np.pi * k / m))
        return out

    def pattern_speeds(self):
        return {f"m={mode['m']}": mode["pattern_speed"] for mode in KHALIL2025_SPIRAL_MODES}

    def provenance(self):
        return {
            **super().provenance(),
            "meaning": "Troughs (deepest potential at each radius) of the m = 2 and m = 3 modes; rigid rotation",
            "radius_range_kpc": list(self.radius_range_kpc),
            "solar_radius_kpc": KHALIL2025_SOLAR_RADIUS_KPC,
            "modes": [dict(mode) for mode in KHALIL2025_SPIRAL_MODES],
        }


@dataclass(frozen=True)
class CastroGinardSpiralArms(SpiralArmModel):
    """The Castro-Ginard et al. (2021) arm segments, each at its own pattern speed.

    The segments are drawn where the paper fitted them (Table 1 azimuth
    ranges), placed about the Sun as in the paper (R_sun = 8.178 kpc), and
    turned about the Galactic centre at their Table 2 (10-50 Myr) speeds.
    """

    step_deg: float = 0.5

    def curves(self, time_myr, *, ro=8.122):
        out = []
        r_sun = CASTRO_GINARD2021_SOLAR_RADIUS_KPC
        for arm, theta_ref, (t0, t1), r_ref, pitch_deg, speed in CASTRO_GINARD2021_ARM_TABLE:
            theta = np.radians(np.linspace(t0, t1, int(np.ceil((t1 - t0) / self.step_deg)) + 1))
            radius = r_ref * np.exp(-(theta - np.radians(theta_ref)) * np.tan(np.radians(pitch_deg)))
            # Heliocentric in the paper's frame (X to the centre, Y along
            # rotation), then about the figure's Galactic centre at ro.
            x_helio = r_sun - radius * np.cos(theta)
            y_helio = radius * np.sin(theta)
            phi = np.arctan2(y_helio, ro - x_helio) + speed * KMS_KPC_TO_RAD_MYR * float(time_myr)
            out.append((arm, np.hypot(ro - x_helio, y_helio), phi))
        return out

    def pattern_speeds(self):
        return {arm: speed for arm, *_rest, speed in CASTRO_GINARD2021_ARM_TABLE}

    def provenance(self):
        return {
            **super().provenance(),
            "meaning": "Present-day arm segments fitted to young open clusters and HMSFRs, each rotated rigidly at its own pattern speed",
            "solar_radius_kpc": CASTRO_GINARD2021_SOLAR_RADIUS_KPC,
            "arms": [
                {"arm": arm, "theta_ref_deg": theta_ref, "theta_range_deg": list(theta_range),
                 "r_ref_kpc": r_ref, "pitch_deg": pitch, "pattern_speed_kms_kpc": speed}
                for arm, theta_ref, theta_range, r_ref, pitch, speed in CASTRO_GINARD2021_ARM_TABLE
            ],
        }


KHALIL2025_ARMS = KhalilSpiralArms(
    name="Spiral arms (Khalil+ 2025)", reference=KHALIL2025_REFERENCE, color="#6ee7b7",
)
CASTRO_GINARD2021_ARMS = CastroGinardSpiralArms(
    name="Spiral arms (Castro-Ginard+ 2021)", reference=CASTRO_GINARD2021_REFERENCE,
    color="#fcd34d", dash="dash",
)


@lru_cache(maxsize=8)
def _plane_transforms(ro: float, zo: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Affine maps from galpy's Galactocentric Cartesian (kpc) to heliocentric
    Galactic and to astropy Galactocentric Cartesian (pc).

    They come from galpy's own SkyCoord conversion (solar height and tilt
    included), the one Oviz uses for every orbit, so arms and clusters share
    one frame exactly.
    """
    from galpy.orbit import Orbit

    # R, phi, z (kpc, rad, kpc): four points that span the space.
    probe = np.array([[ro, 0.0, 0.0], [ro + 1.0, 0.0, 0.0], [ro, 0.25, 0.0], [ro, 0.0, 0.5]])
    vxvv = np.column_stack([probe[:, 0] / ro, np.zeros(4), np.ones(4), probe[:, 2] / ro, np.zeros(4), probe[:, 1]])
    sky = Orbit(vxvv, ro=ro, vo=236.0, zo=zo, solarmotion="schoenrich").SkyCoord()
    src = np.column_stack([probe[:, 0] * np.cos(probe[:, 1]), probe[:, 0] * np.sin(probe[:, 1]), probe[:, 2], np.ones(4)])
    maps = []
    for xyz in (sky.galactic.cartesian.xyz.to_value(u.pc), sky.galactocentric.cartesian.xyz.to_value(u.pc)):
        solution = np.linalg.solve(src, xyz.T)  # rows: A^T (3) then b
        maps.extend([solution[:3].T, solution[3]])
    return tuple(maps)


def spiral_arm_coordinates(model: SpiralArmModel, time_myr: float, *, ro: float = 8.122,
                           zo: float = 0.0208) -> tuple[np.ndarray, np.ndarray]:
    """One arm model at ``time_myr`` in the Galactic plane.

    Returns heliocentric Galactic and astropy Galactocentric Cartesian
    positions (pc), each shaped (3, N), with a NaN column between arms. The
    number of points never changes with time, so every frame's line has the
    same vertices and the viewer interpolates between frames smoothly.
    """
    columns = []
    for _arm, radius, phi in model.curves(time_myr, ro=ro):
        if columns:
            columns.append(np.full((3, 1), np.nan))
        columns.append(np.vstack([radius * np.cos(phi), radius * np.sin(phi), np.zeros_like(radius)]))
    plane = np.hstack(columns)
    helio_a, helio_b, gc_a, gc_b = _plane_transforms(float(ro), float(zo))
    return helio_a @ plane + helio_b[:, None], gc_a @ plane + gc_b[:, None]


def provenance(models=(), *, ro: float = 8.122, vo: float = 236.0) -> dict:
    """What a figure using :func:`khalil2025_potential` and these arm models shows."""
    return {
        "potential": "MWPotential2014 + Khalil et al. (2025) m = 2 and m = 3 spirals (galpy SpiralArmsPotential)",
        "reference": KHALIL2025_REFERENCE,
        "ro_kpc": ro,
        "vo_kms": vo,
        "spiral_modes": [dict(mode) for mode in KHALIL2025_SPIRAL_MODES],
        "spiral_height_kpc": KHALIL2025_SPIRAL_HEIGHT_KPC,
        "spiral_scale_length_kpc": KHALIL2025_SPIRAL_SCALE_LENGTH_KPC,
        "approximation": (
            "Matches SPIBACK's spiral potential in the plane; adds galpy's Cox & Gomez vertical profile; "
            "no bar, no hard radial cutoffs, MWPotential2014 instead of the paper's AGAMA background; "
            "spirals fully grown and steady over the timeline"
        ),
        "reference_frame": "The LSR orbit is integrated in the axisymmetric part only, so it stays circular",
        "arm_models": [model.provenance() for model in models],
    }


__all__ = [
    "CASTRO_GINARD2021_ARMS",
    "CASTRO_GINARD2021_ARM_TABLE",
    "CastroGinardSpiralArms",
    "KHALIL2025_ARMS",
    "KHALIL2025_SPIRAL_MODES",
    "KhalilSpiralArms",
    "SpiralArmModel",
    "khalil2025_potential",
    "khalil2025_spiral_potentials",
    "provenance",
    "spiral_arm_coordinates",
]
