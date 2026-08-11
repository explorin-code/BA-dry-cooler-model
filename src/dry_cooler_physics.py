"""
dry_cooler_physics.py
======================
Cooler geometry: tube/fin dimensions and every derived area, ratio and
pitch the heat-transfer correlations need. Raw inputs come from
parameters.py; everything else here is derived from those.
"""

from dataclasses import dataclass
import numpy as np

from src.param_loader import from_parameters


@dataclass
class Geometry:
    # --- Raw input -----------------------------------------------------
    d: float                       # tube outer diameter [m]
    s: float                       # fin thickness [m]
    a: float                       # fin spacing [m]
    d_i: float                     # tube inner diameter [m]
    n_tubes: int                   # number of tubes
    n_rows: int                    # number of tube rows
    t_q: float                     # tube pitch, transverse [m]
    t_l: float                     # tube pitch, longitudinal [m]
    height: float                  # cooler height == single tube-pass length [m]
    material: str = 'Aluminum'     # fin material

    # --- Intermediate ----------------------------------------------------
    @property
    def fin_density(self) -> float:
        """Fins per meter, derived from fin pitch (spacing + thickness)."""
        return 1.0 / self.t_R

    @property
    def fins_per_pipe(self) -> int:
        return int(self.fin_density * self.height)

    @property
    def t_R(self) -> float:
        """Fin pitch [m]."""
        return self.a + self.s

    @property
    def width(self) -> float:
        return self.n_tubes * self.t_q

    @property
    def inflow_cross_section(self) -> float:
        return self.width * self.height

    @property
    def lambda_R(self) -> float:
        """Fin material thermal conductivity at room temp [W/m-K]."""
        solid_conductivities = {
            'Aluminum':     237.0,
            'Copper':       401.0,
            'Carbon Steel':  50.0,
        }
        return solid_conductivities.get(self.material, 237.0)

    # --- Outputs -----------------------------------------------------------
    @property
    def A_R(self) -> float:
        """Fin surface area on one tube [m²]."""
        return 2 * (self.t_q * self.t_l - (np.pi * self.d**2) / 4) * self.fins_per_pipe

    @property
    def A_G(self) -> float:
        """Exposed base tube area between fins [m²]."""
        return (self.fins_per_pipe + 1) * np.pi * self.d * self.a

    @property
    def A(self) -> float:
        """Total outer surface area [m²]."""
        return self.A_R + self.A_G

    @property
    def A_i(self) -> float:
        """Inner surface area of one tube [m²]."""
        return self.height * self.d_i * np.pi

    @property
    def A_flow_coolant(self) -> float:
        """Tube-side flow area across all tubes, 1-pass [m²]."""
        return self.n_tubes * (np.pi / 4) * (self.d_i ** 2)

    @property
    def Ao_Ae_ratio(self) -> float:
        """Ratio of total to minimum airflow cross-section [-].
        Source: [not given in original notes -- TODO: find & fill in]."""
        numerator = self.t_q * (self.a + self.s)
        denominator = (self.t_q - self.d) * self.a
        return numerator / denominator

    @property
    def A_Go(self) -> float:
        """Bare tube surface area per element [m²].
        Source: VDI Heat Atlas, Section M1, A_Go = pi * d * height."""
        return np.pi * self.d * self.height

    @property
    def l(self) -> float:
        """Single tube length across all n_rows passes [m]."""
        return self.n_rows * self.height


def get_geometry() -> Geometry:
    """Current cooler design -- edit parameters.py to change the numbers."""
    return from_parameters(Geometry, {
        'd': 'D_TUBE_OUTER',
        's': 'FIN_THICKNESS',
        'a': 'FIN_SPACING',
        'd_i': 'D_TUBE_INNER',
        'n_tubes': 'N_TUBES',
        'n_rows': 'N_ROWS',
        't_q': 'T_Q',
        't_l': 'T_L',
        'height': 'HEIGHT',
        'material': 'MATERIAL',
    })