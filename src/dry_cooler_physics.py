"""
dry_cooler_physics.py
======================
Cooler geometry: tube/fin dimensions and every derived area, ratio and
pitch the heat-transfer correlations need. `Geometry` is a plain data
container -- all the "raw input" fields are the actual design knobs;
everything under "intermediate" / "outputs" is computed from those.

`get_geometry()` is the single source of truth for the current cooler
design -- change the numbers there to try a different geometry.
"""

from dataclasses import dataclass
from math import pi
import CoolProp.CoolProp as CP
import numpy as np


@dataclass
class Geometry:

    # -------------------------------------------------------------------
    # 1. Raw input -- the actual design variables
    # -------------------------------------------------------------------
    # OBSOLETE (safe to delete): `D: float` used to live here -- outer
    # diameter of a circular fin. Irrelevant since the switch to
    # continuous (rectangular) fins; A_R below was reworked accordingly.
    d: float                       # Outer diameter of the tube [m]
    s: float                       # Fin thickness [m]
    a: float                       # Fin spacing [m]
    d_i: float                     # Inner diameter of the tube [m]
    material: str                  # Fin material
    n_tubes: int                   # Number of tubes -- vary this to scale the array
    n_rows: int                    # Number of tube rows
    t_q: float                     # Tube pitch (center-to-center) [m] -- "quer" (transverse)
    t_l: float                     # Tube pitch (center-to-center) [m] -- "längs" (longitudinal)
    height: float                  # Height of the cooler [m] == length of a single tube pass

    # -------------------------------------------------------------------
    # 2. Intermediate properties
    # -------------------------------------------------------------------

    @property
    def fin_density(self) -> float:
        """Number of fins per meter (fixed at 9 fins/inch)."""
        return 9 * 2.54 * 100      # 9 fins/inch * 2.54 cm/inch * 100 cm/m

    @property
    def fins_per_pipe(self) -> int:
        """Number of fins per tube, derived from fin density and cooler height."""
        return int(self.fin_density * self.height)

    @property
    def t_R(self) -> float:
        """Fin pitch [m] (left edge to next left edge).
        NOTE: not currently referenced by any correlation in this codebase
        (checked July 2026) -- kept because it's a natural companion value
        to fin_density/fins_per_pipe. Safe to delete if you're sure you
        won't need fin pitch itself somewhere."""
        return self.a + self.s

    @property
    def width(self) -> float:
        """Width of the cooler [m]. Scales directly with n_tubes since t_q
        (tube pitch) is held fixed -- this is how the array grows/shrinks
        when varying the number of tubes."""
        return self.n_tubes * self.t_q

    @property
    def inflow_cross_section(self) -> float:
        """Cross-sectional area for the inflow [m²]. Derived from width
        (which scales with n_tubes at fixed t_q) and the fixed cooler
        height, so varying n_tubes changes the frontal area while tube
        pitch and tube-pass length (height) stay constant."""
        return self.width * self.height

    @property
    def lambda_R(self) -> float:
        """Thermal conductivity of the fin material at room temperature [W/m-K]."""
        # CoolProp only does fluids! We use a static lookup for solid metals.
        solid_conductivities = {
            'Aluminum':     237.0,     # W/m-K (standard pure aluminum)
            'Copper':       401.0,     # W/m-K
            'Carbon Steel':  50.0,     # W/m-K
        }
        # Look up the material, default to 237.0 (Aluminum) if there's a typo
        return solid_conductivities.get(self.material, 237.0)

    # -------------------------------------------------------------------
    # 3. Outputs
    # -------------------------------------------------------------------

    @property
    def A_R(self) -> float:
        """Outer surface area of fin segments on one tube [m²].
        Rectangular fin area (not the classic circular-fin formula)."""
        return 2 * (self.t_q * self.t_l - (np.pi * self.d**2) / 4) * self.fins_per_pipe

    @property
    def A_G(self) -> float:
        """Exposed base tube area between fins [m²]."""
        return (self.fins_per_pipe + 1) * np.pi * self.d * self.a

    @property
    def A(self) -> float:
        """Total outer surface area product [m²]."""
        return self.A_R + self.A_G

    @property
    def A_i(self) -> float:
        """Inner surface area of one tube [m²]. Uses fixed height."""
        return self.height * self.d_i * np.pi

    @property
    def A_flow_coolant(self) -> float:
        """Inner (tube-side) cross-sectional flow area across all tubes,
        assuming a 1-pass tube arrangement [m²]. Scales with n_tubes, which
        is exactly the lever you want: more tubes at the same total mass
        flow rate lowers velocity/Re per tube -- fewer tubes raises it."""
        return self.n_tubes * (np.pi / 4) * (self.d_i ** 2)

    @property
    def Ao_Ae_ratio(self) -> float:
        """Ratio of total to minimum cross-sectional area for airflow [-].
        Adapted for continuous plate fins (Blocklamellen). Depends only on
        t_q (fixed) and tube/fin dimensions.
        Source: [not given in original notes -- TODO: find & fill in]."""
        numerator = self.t_q * (self.a + self.s)
        denominator = (self.t_q - self.d) * self.a
        return numerator / denominator

    @property
    def A_Go(self) -> float:
        """Bare tube surface area per element [m²].
        Source: VDI Heat Atlas, Section M1, p. [not recorded -- TODO], A_Go = pi * d * height."""
        return np.pi * self.d * self.height

    @property
    def l(self) -> float:
        """Length of a single tube, accounting for n_rows passes through
        the cooler. Ignores bend radii for now."""
        return self.n_rows * self.height


def get_geometry() -> Geometry:
    """Single source of truth for the current cooler design. Edit the
    numbers below to try a different geometry."""
    return Geometry(
        # OBSOLETE (safe to delete): `D=0.056,` used to sit here -- default
        # for the old circular-fin diameter field, see note on Geometry above.
        d=0.012,                   # 0.0254 in an earlier iteration
        s=0.00012,
        a=0.0023,
        d_i=0.011,                 # 12 mm OD with 0.5 mm wall thickness (0.021 in an earlier iteration)
        material='Aluminum',
        n_tubes=11,                # 17 in an earlier iteration
        n_rows=6,
        t_q=0.03,                  # 0.06 in an earlier iteration
        t_l=0.03,                  # added for fin-efficiency calculation
        height=1.0,
    )