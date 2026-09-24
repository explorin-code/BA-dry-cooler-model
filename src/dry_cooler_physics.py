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


# Solid thermal conductivities at room temp [W/m-K], shared by fin and pipe.
SOLID_CONDUCTIVITIES = {
    'Aluminum':     237.0,
    'Copper':       401.0,
    'Carbon Steel':  50.0,
}


def _solid_conductivity(material: str) -> float:
    """Looks up a material's conductivity; raises on unknown names instead
    of silently falling back to a default."""
    try:
        return SOLID_CONDUCTIVITIES[material]
    except KeyError:
        raise ValueError(f"Unknown material {material!r} -- known: {list(SOLID_CONDUCTIVITIES)}") from None


@dataclass
class Geometry:
    # --- Raw input -----------------------------------------------------
    d_a: float                     # tube outer diameter [m]
    delta_R: float                 # fin thickness [m]
    t_s: float                     # fin spacing [m]
    d_i: float                     # tube inner diameter [m]
    n_tubes: int                   # number of tubes
    n_rows: int                    # number of tube rows
    s_1: float                     # tube pitch, transverse [m]
    s_2: float                     # tube pitch, longitudinal [m]
    height: float                  # cooler height == single tube-pass length [m]
    pipe_material: str = 'Copper'     # pipe material
    fin_material: str = 'Aluminum'      # fin material
    h_collar: float = 0.0          # collar height [m] -- no data available yet

    # --- Intermediate ----------------------------------------------------
    # AI-REVIEW: this replaced a hardcoded "9 fins/inch" constant that was
    # numerically inconsistent with fin spacing (t_s) + fin thickness
    # (delta_R). Confirm which one -- the original constant, or this
    # derived formula -- actually matches the physical hardware. See
    # CLAUDE.md.
    @property
    def n_R(self) -> float:
        """Fins per meter, derived from fin pitch (spacing + thickness)."""
        return 1.0 / self.t_R

    @property
    def fins_per_pipe(self) -> int:
        return int(self.n_R * self.height)

    @property
    def t_R(self) -> float:
        """Fin pitch [m]."""
        return self.t_s + self.delta_R

    @property
    def width(self) -> float:
        return self.n_tubes * self.s_1

    @property
    def inflow_cross_section(self) -> float:
        return self.width * self.height

    @property
    def lambda_f(self) -> float:
        """Fin material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.fin_material)

    @property
    def lambda_p(self) -> float:
        """Pipe (tube wall) material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.pipe_material)

    # --- Outputs -----------------------------------------------------------
    @property
    def A_R(self) -> float:
        """Fin surface area on one tube [m²]."""
        return 2 * (self.s_1 * self.s_2 - (np.pi * self.d_a**2) / 4) * self.fins_per_pipe

    @property
    def A_G(self) -> float:
        """Exposed base tube area between fins [m²]."""
        return (self.fins_per_pipe + 1) * np.pi * self.d_a * self.t_s

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
        Source: [not given in original notes -- TODO: find & fill in].
        Checked 2026-08 against VDI Wärmeatlas Chapter M1 without finding
        it there; next leads are L1.5 or N4 (Kühltürme). Likely lead: this
        is the inverse of Kays & London's sigma = Ac/Afr for a plate-fin
        core (min. gap transverse), with d_c as the blocking diameter.
        Uses d_c (== d_a while h_collar = 0) since that is what actually
        blocks the flow between fins."""
        numerator = self.s_1 * (self.t_s + self.delta_R)
        denominator = (self.s_1 - self.d_c) * self.t_s
        return numerator / denominator

    @property
    def A_total(self) -> float:
        """Total outer (air-side) surface area across the whole array [m²]."""
        return self.A * self.n_tubes * self.n_rows

    @property
    def A_c_air(self) -> float:
        """Minimum free-flow area on the air side [m²]."""
        return self.inflow_cross_section / self.Ao_Ae_ratio

    @property
    def A_Go(self) -> float:
        """Bare tube surface area per element [m²].
        Source: VDI Heat Atlas, Section M1, A_Go = pi * d_a * height."""
        return np.pi * self.d_a * self.height

    @property
    def l(self) -> float:
        """Single tube length across all n_rows passes [m] -- this is the
        book's "L" (Rohrlänge) used as the characteristic length in the
        coolant-side Nu correlations; deliberately kept lowercase/spelled
        differently from the `height` field (a different quantity, the
        single-pass height) to avoid an L/l case-only collision."""
        return self.n_rows * self.height

    @property
    def d_c(self) -> float:
        """Fin-collar outside diameter [m]. Source: Wang, Chi & Chang,
        Part II, Nomenclature, p. 2694: Dc = fin collar outside diameter."""
        return self.d_a + 2.0 * self.h_collar


def get_geometry() -> Geometry:
    """Current cooler design -- edit parameters.py to change the numbers."""
    return from_parameters(Geometry, {
        'd_a': 'D_TUBE_OUTER',
        'delta_R': 'FIN_THICKNESS',
        't_s': 'FIN_SPACING',
        'd_i': 'D_TUBE_INNER',
        'n_tubes': 'N_TUBES',
        'n_rows': 'N_ROWS',
        's_1': 'S_1',
        's_2': 'S_2',
        'height': 'HEIGHT',
        'pipe_material': 'PIPE_MATERIAL',
        'fin_material': 'FIN_MATERIAL',
    })
