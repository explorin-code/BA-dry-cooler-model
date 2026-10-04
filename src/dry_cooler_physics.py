"""
dry_cooler_physics.py
======================
Cooler geometry: tube/fin dimensions and every derived area, ratio and
pitch the heat-transfer correlations need. Raw inputs come from
parameters.py; everything else here is derived from those.
"""

from dataclasses import dataclass
import numpy as np

from src.param_loader import from_parameters, MM


# Solid thermal conductivities [W/m-K], shared by fin and pipe.
# Source: VDI-Wärmeatlas, 12th ed. (2019), Chapter D6 (Cu 380, Al 160 --
# technical copper / aluminium alloys, not the pure-metal 401 / 237).
# Carbon steel 50: NOT yet checked against D6.
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
    F_p: float                     # fin pitch (spacing + thickness) [m] -- Wang's F_p
    d_i: float                     # tube inner diameter [m]
    N_t: int                   # number of tubes
    N_r: int                       # number of tube rows (= coolant passes N_p, one row per pass)
    P_t: float                     # tube pitch, transverse [m]
    P_l: float                     # tube pitch, longitudinal [m]
    L_p: float                     # single tube-pass length (= cooler height) [m]
    pipe_material: str = 'Copper'     # pipe material
    fin_material: str = 'Aluminum'      # fin material
    h_collar: float = 0.0          # fin-collar thickness [m] -- deliberately 0: no fin collar is
                                   # modelled (D_c = d), a modelling decision, not missing data

    # --- Intermediate ----------------------------------------------------
    # AI-REVIEW: this replaced a hardcoded "9 fins/inch" constant that was
    # numerically inconsistent with fin pitch F_p (spacing t_s + fin thickness)
    # (delta_f). Confirm which one -- the original constant, or this
    # derived formula -- actually matches the physical hardware. See
    # CLAUDE.md.
    @property
    def n_R(self) -> float:
        """Fins per meter, derived from fin pitch (spacing + thickness)."""
        return 1.0 / self.F_p

    @property
    def fins_per_pipe(self) -> int:
        return int(self.n_R * self.L_p)

    @property
    def t_s(self) -> float:
        """Clear spacing between two fins [m] (fin pitch minus fin thickness)."""
        return self.F_p - self.delta_f

    @property
    def width(self) -> float:
        return self.N_t * self.P_t

    @property
    def A_fr(self) -> float:
        return self.width * self.L_p

    @property
    def lambda_f(self) -> float:
        """Fin material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.fin_material)

    @property
    def lambda_p(self) -> float:
        """Pipe (tube wall) material thermal conductivity at room temp [W/m-K]."""
        return _solid_conductivity(self.pipe_material)

    # --- Outputs -----------------------------------------------------------
    # Areas per tube and row, as in the VDI-Wärmeatlas (2019) M1 worked
    # example (p. 1689: A_f, A_p, A, A_p0, A_i), with the circular fin's
    # face replaced by the plate-fin share P_t*P_l - pi*d^2/4.
    @property
    def A_f(self) -> float:
        """Fin surface area on one tube [m²]."""
        return 2 * (self.P_t * self.P_l - (np.pi * self.d**2) / 4) * self.fins_per_pipe

    @property
    def A_p(self) -> float:
        """Exposed base tube area between fins [m²]."""
        return (self.fins_per_pipe + 1) * np.pi * self.d * self.t_s

    @property
    def A(self) -> float:
        """Total outer surface area [m²]."""
        return self.A_f + self.A_p

    @property
    def A_i(self) -> float:
        """Inner surface area of one tube [m²]."""
        return self.L_p * self.d_i * np.pi

    @property
    def A_cs_c(self) -> float:
        """Tube-side flow area across all tubes, 1-pass [m²]."""
        return self.N_t * (np.pi / 4) * (self.d_i ** 2)

    @property
    def Afr_Ae_ratio(self) -> float:
        """Ratio of frontal to narrowest airflow cross-section A_o/A_e [-].
        Source: VDI Waermeatlas, Chapter M1, p. 1689, worked example
        ("Verengter Stroemungsquerschnitt"), for circular fins:
            A_o/A_e = t_q (a + s) / ((t_q - d) a + (t_q - D) s)
        with t_q = P_t, a = fin spacing t_s, s = fin thickness delta_f.
        Adapted to continuous plate fins: the fin spans the whole transverse
        pitch (D = t_q), so the fin-band term (t_q - D) s vanishes. Equals
        Kays & London's 1/sigma (sigma = A_e/A_fr). Uses D_c (== d while
        h_collar = 0) as the blocking diameter. Only the transverse gap is
        considered -- for staggered banks the diagonal gap can be narrower
        at small P_l (not the case for the current geometry)."""
        numerator = self.P_t * self.F_p
        denominator = (self.P_t - self.D_c) * self.t_s
        return numerator / denominator

    @property
    def A_tot(self) -> float:
        """Total outer (air-side) surface area across the whole array [m²]."""
        return self.A * self.N_t * self.N_r

    @property
    def A_e(self) -> float:
        """Minimum free-flow area on the air side [m²]."""
        return self.A_fr / self.Afr_Ae_ratio

    @property
    def A_p0(self) -> float:
        """Bare tube surface area per element [m²].
        Source: VDI Heat Atlas, Section M1, A_p0 = pi * d * height."""
        return np.pi * self.d * self.L_p

    @property
    def l(self) -> float:
        """Single tube length across all N_r passes [m] -- this is the
        book's "L" (Rohrlänge) used as the characteristic length in the
        coolant-side Nu correlations; deliberately kept lowercase/spelled
        differently from the `height` field (a different quantity, the
        single-pass height) to avoid an L/l case-only collision."""
        return self.N_r * self.L_p

    @property
    def D_c(self) -> float:
        """Fin-collar outside diameter [m]. Source: Wang, Chi & Chang,
        Part II, Nomenclature, p. 2694: Dc = fin collar outside diameter."""
        return self.d + 2.0 * self.h_collar


def get_geometry() -> Geometry:
    """Current cooler design -- edit parameters.py to change the numbers."""
    return from_parameters(Geometry, {
        'd': ('D_TUBE_OUTER_MM', MM),
        'delta_f': ('FIN_THICKNESS_MM', MM),
        'F_p': ('FIN_PITCH_MM', MM),
        'd_i': ('D_TUBE_INNER_MM', MM),
        'N_t': 'N_TUBES',
        'N_r': 'N_ROWS',
        'P_t': ('PITCH_TRANSVERSE_MM', MM),
        'P_l': ('PITCH_LONGITUDINAL_MM', MM),
        'L_p': ('HEIGHT_MM', MM),
        'pipe_material': 'PIPE_MATERIAL',
        'fin_material': 'FIN_MATERIAL',
    })
