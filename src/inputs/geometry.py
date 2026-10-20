"""
geometry.py
======================
Cooler geometry: tube/fin dimensions and every derived area, ratio and
pitch the heat-transfer correlations need. Raw inputs come from
parameters.py; everything else here is derived from those.
"""

from dataclasses import dataclass
import numpy as np

from src.inputs.cooler_sets import CELL_SETS
from src.inputs.param_loader import MM


# Solid thermal conductivities [W/m-K], shared by fin and pipe.
# Source: VDI-Wärmeatlas, 12th ed. (2019), Chapter D6, Tab. 8 (design values
# from DIN EN ISO 10456): aluminium / aluminium alloys 160, copper 380,
# steel (structural / reinforcing) 50 -- technical alloys, not the
# pure-metal 401 / 237.
SOLID_CONDUCTIVITIES = {
    'Aluminum':     160.0,
    'Copper':       380.0,
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
    d: float                     # tube outer diameter [m]
    delta_f: float                 # fin thickness [m]
    s_f: float                     # fin pitch (spacing + thickness) [m] -- Wang's s_f
    d_i: float                     # tube inner diameter [m]
    N_t: int                   # number of tubes
    N_r: int                       # number of tube rows (= coolant passes N_p, one row per pass)
    s_t: float                     # tube pitch, transverse [m]
    s_l: float                     # tube pitch, longitudinal [m]
    l: float                       # length of one tube pass (= cooler height) [m]
    pipe_material: str = 'Copper'     # pipe material
    fin_material: str = 'Aluminum'      # fin material
    delta_c: float = 0.0           # fin-collar thickness [m] -- currently 0 (no collar, d_c = d);
                                   # only enters Wang's d_c (friction factor, narrowest gap)

    # --- Intermediate ----------------------------------------------------
    # The fin count follows from the fin pitch s_f (see docs/model-notes.md).
    @property
    def N_f(self) -> int:
        """Number of fins along one tube pass (whole fins only)."""
        return int(self.l / self.s_f)

    @property
    def delta_p(self) -> float:
        """Tube wall thickness [m]."""
        return (self.d - self.d_i) / 2

    # Fin rectangle around one tube for the fin efficiency (staggered bank,
    # the modelled layout). Source: VDI M1, p. 1687, Eq. (14); book notation
    # s_1/s_2 == s_t/s_l. (An in-line bank would use b_f = min(s_t, s_l),
    # l_f = max(s_t, s_l), Eq. (13) -- see calc_fin_efficiency_inline.)
    @property
    def b_f(self) -> float:
        """Fin rectangle width [m] (staggered)."""
        return self.s_t if self.s_l >= self.s_t / 2 else 2 * self.s_l

    @property
    def l_f(self) -> float:
        """Fin rectangle length [m] (staggered), l_f >= b_f."""
        return (self.s_l**2 + (self.s_t / 2)**2)**0.5

    @property
    def d_c(self) -> float:
        """Fin-collar outside diameter [m]. Source: Wang, Chi & Chang,
        Part II, Nomenclature, p. 2694: Dc = fin collar outside diameter."""
        return self.d + 2.0 * self.delta_c

    @property
    def b(self) -> float:
        """Cooler width (frontal, across the tubes) [m]."""
        return self.N_t * self.s_t

    @property
    def A_fr(self) -> float:
        """Frontal (face) area [m²]."""
        return self.b * self.l

    @property
    def l_tot(self) -> float:
        """Total tube length of one circuit across all N_r passes [m] -- the
        tube length l used as the characteristic length in the coolant-side
        Nu correlations (VDI G1)."""
        return self.N_r * self.l

    @property
    def lambda_f(self) -> float:
        """Fin material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.fin_material)

    @property
    def lambda_p(self) -> float:
        """Pipe (tube wall) material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.pipe_material)

    # --- Areas per tube pass -------------------------------------------------
    # One tube in one row over the pass length l, as in the VDI-Wärmeatlas
    # (2019) M1 worked example (p. 1689: A_f, A_p, A, A_p0, A_c), with the
    # circular fin's face replaced by the plate-fin share s_t*s_l - pi*d^2/4.
    @property
    def A_f(self) -> float:
        """Fin surface area per tube pass (both fin faces) [m²]."""
        return 2 * (self.s_t * self.s_l - (np.pi * self.d**2) / 4) * self.N_f

    @property
    def A_p(self) -> float:
        """Exposed tube surface between the fins per tube pass [m²]: tube
        length minus the fin feet."""
        return np.pi * self.d * (self.l - self.N_f * self.delta_f)

    @property
    def A(self) -> float:
        """Outer (air-side) surface area per tube pass [m²]."""
        return self.A_f + self.A_p

    @property
    def A_c(self) -> float:
        """Inner (coolant-side) surface area per tube pass [m²]."""
        return self.l * self.d_i * np.pi

    @property
    def A_p0(self) -> float:
        """Bare tube surface area per tube pass (as if unfinned) [m²].
        Source: VDI Heat Atlas, Section M1, A_p0 = pi * d * height."""
        return np.pi * self.d * self.l

    # --- Whole-exchanger areas -----------------------------------------------
    @property
    def A_tot(self) -> float:
        """Total outer (air-side) surface area across the whole array [m²]."""
        return self.A * self.N_t * self.N_r

    @property
    def A_cs_c(self) -> float:
        """Tube-side flow area across all tubes, 1-pass [m²]."""
        return self.N_t * (np.pi / 4) * (self.d_i ** 2)

    @property
    def sigma(self) -> float:
        """Ratio of narrowest to frontal airflow cross-section A_e/A_fr [-]
        (Kays & London's sigma). Source: VDI Waermeatlas, Chapter M1,
        p. 1689, worked example ("Verengter Stroemungsquerschnitt"), for
        circular fins (inverted):
            A_o/A_e = t_q (a + s) / ((t_q - d) a + (t_q - D) s)
        with VDI's t_q = s_t, a = clear fin spacing s_f - delta_f, s = fin
        thickness delta_f.
        Adapted to continuous plate fins: the fin spans the whole transverse
        pitch (D = t_q), so the fin-band term (t_q - D) s vanishes. Uses d_c
        (== d while delta_c = 0) as the blocking diameter. Only the transverse
        gap is considered -- for staggered banks the diagonal gap can be
        narrower at small s_l (not the case for the current geometry)."""
        return (self.s_t - self.d_c) * (self.s_f - self.delta_f) / (self.s_t * self.s_f)

    @property
    def A_e(self) -> float:
        """Minimum free-flow area on the air side [m²]."""
        return self.sigma * self.A_fr


# parameters.py constant -> (CellSet attribute, divisor): the elementary cell, preset + overrides
CELL_PARAMETERS = {
    'D_TUBE_OUTER_MM': 'd_mm', 'D_TUBE_INNER_MM': 'd_i_mm', 'FIN_THICKNESS_MM': 'delta_f_mm',
    'FIN_PITCH_MM': 's_f_mm', 'PITCH_TRANSVERSE_MM': 's_t_mm', 'PITCH_LONGITUDINAL_MM': 's_l_mm',
    'FIN_COLLAR_THICKNESS_MM': 'delta_c_mm', 'PIPE_MATERIAL': 'pipe_material', 'FIN_MATERIAL': 'fin_material',
}


def get_cell_values() -> dict:
    """The elementary cell of the current run: the preset cell (parameters.GEOMETRY_SET, src/inputs/cooler_sets.py)
    with every parameters.py constant that is not None replacing the preset's value. {CellSet attribute: value}."""
    import src.parameters as parameters
    from src.inputs.cooler_sets import CELL_SETS
    preset = CELL_SETS[parameters.GEOMETRY_SET] if parameters.GEOMETRY_SET else None
    values = {}
    for constant, attribute in CELL_PARAMETERS.items():
        override = getattr(parameters, constant)
        if override is not None:
            values[attribute] = override
        elif preset is not None:
            values[attribute] = getattr(preset, attribute)
        else:
            raise ValueError(f"GEOMETRY_SET is None, so {constant} must be given.")
    return values


def get_geometry() -> Geometry:
    """Current cooler: the elementary cell (preset + overrides, see get_cell_values) in a module of
    N_TUBES (None: the preset's reference count) x N_ROWS tubes of length HEIGHT_MM -- edit parameters.py."""
    import src.parameters as parameters
    c = get_cell_values()
    N_t = parameters.N_TUBES or CELL_SETS[parameters.GEOMETRY_SET].N_t_ref     # fixed per preset, not tied to the inputs
    return Geometry(
        d=c['d_mm'] / MM, delta_f=c['delta_f_mm'] / MM, s_f=c['s_f_mm'] / MM, d_i=c['d_i_mm'] / MM,
        N_t=N_t, N_r=parameters.N_ROWS, s_t=c['s_t_mm'] / MM, s_l=c['s_l_mm'] / MM,
        l=parameters.HEIGHT_MM / MM, pipe_material=c['pipe_material'], fin_material=c['fin_material'],
        delta_c=c['delta_c_mm'] / MM)
