"""
config.py -- the settings of a cooler comparison, and its research constants.

Everything the user can change lives in parameters.py (groups "Comparison" / "Analyses", shown in the GUI's comparison
mode). `settings()` reads it ONCE per run into a `Settings` object that is handed to every function that needs it
(`cfg`) -- nothing in the comparison reads parameters.py behind the scenes, so a run is reproducible from its `cfg`.

What is NOT a setting is defined here: the literature ranges of the sensitivity study and the material densities.
The coolers themselves are src/cooler_sets.py.
"""

import dataclasses
from dataclasses import dataclass

import src.parameters as parameters
from src import weather
from src.cooler_sets import CELL_SETS
from src.solver_settings import SolverSettings, get_solver_settings

# The water target is met when T_c,o <= target + TARGET_TOL: the variable-speed fan lands on the target within its root-finder
# tolerance (face velocity to 1e-3 m/s = a few mK), far below the accuracy of the model. ONE definition for every flag and count.
TARGET_TOL = 0.01             # [K]

# Uniform inlet air: Cell 2D == 3D. 60 NTU elements: the NTU table's R = C_a/C_c range is 0.1 * N_e = 6 (Eq. 1 validity);
# the 1 kg/s designs sit at R ~ 2-5, beyond 6 (high air flows) the table falls back to the field solve (correct, ~100x slower).
NTU_ELEMENTS = 60


@dataclass(frozen=True)
class Settings:
    """One cooler comparison: coolers, water, operation, solver, prices (parameters.py, group Comparison + prices)."""
    set_keys: tuple                 # the compared coolers (src/cooler_sets.py keys)
    solver_key: str                 # one of solvers.SOLVER_LABELS
    solver_settings: SolverSettings  # tuning of the solvers (parameters.py, group Solver tuning), Cell 2D, 60 NTU elements
    m_dot_c: float                  # water per module [kg/s]
    T_c_i: float                    # water inlet [°C]
    T_c_target: float               # target water outlet [°C]
    n_rows: int
    n_tubes: int | None             # None: each cooler's fixed reference count (cooler_sets N_t_ref)
    l_tube: float                   # tube length [m]
    Re_c_design: float              # reference tube Reynolds number of the tubes / wall studies
    w_fr_design: float              # fixed fan speed = top of the variable range [m/s]
    w_fr_min: float                 # slowest fan setting [m/s]
    T_freeze_risk: float            # water outlet below this = freezing risk [°C]
    pad_in_dry_path: bool           # PAD_IN_DRY_AIR_PATH
    design_T_air: float             # design ambient [°C]
    design_phi: float               # design ambient relative humidity [0-1]
    year_stride: int                # year analyses: every n-th hour
    sensitivity_set: str
    electricity_eur_per_kwh: float
    water_eur_per_m3: float         # incl. the wastewater fee
    copper_eur_per_kg: float
    aluminium_eur_per_kg: float
    lifetime_years: float

    def cell(self, key):
        """Cooler `key` (any of src/cooler_sets.py, compared or not) as a module cell with the comparison's rows."""
        return dataclasses.replace(CELL_SETS[key], N_r=self.n_rows)

    @property
    def sets(self) -> dict:
        """The compared coolers {key: CellSet}."""
        return {key: self.cell(key) for key in self.set_keys}

    def design_points(self) -> dict:
        """Representative ambient states of the climate year + the design ambient."""
        return weather.design_points(weather.add_wet_bulb(weather.load_year()),
                                     design_T=self.design_T_air, design_phi=self.design_phi)

    def design_point(self) -> dict:
        return self.design_points()['design']

    def solver(self, table_ok=True) -> str:
        """The solver of this comparison. 'ntu_table' is fast for years, but needs a table per row count: analyses that
        change the rows or the cell on every point pass table_ok=False and get the table-free 'ntu_field' instead."""
        return 'ntu_field' if (self.solver_key == 'ntu_table' and not table_ok) else self.solver_key


def settings() -> Settings:
    """The cooler-comparison settings of parameters.py as they are NOW (call once per run)."""
    p = parameters
    return Settings(
        set_keys=tuple(p.COMPARE_SETS), solver_key=p.COMPARE_SOLVER,
        solver_settings=dataclasses.replace(get_solver_settings(), cell_2d=True, ntu_n_elements=NTU_ELEMENTS),
        m_dot_c=p.COMPARE_M_COOLANT, T_c_i=p.T_COOLANT_IN, T_c_target=p.T_COOLANT_TARGET_OUT,
        n_rows=int(p.COMPARE_N_ROWS), n_tubes=p.COMPARE_N_TUBES, l_tube=p.COMPARE_TUBE_LENGTH_MM / 1000.0,
        Re_c_design=p.COMPARE_RE_TUBE, w_fr_design=p.COMPARE_W_FR_DESIGN, w_fr_min=p.COMPARE_W_FR_MIN,
        T_freeze_risk=p.COMPARE_T_FREEZE_RISK, pad_in_dry_path=p.PAD_IN_DRY_AIR_PATH,
        design_T_air=p.COMPARE_DESIGN_T_AIR, design_phi=p.COMPARE_DESIGN_PHI_PERCENT / 100.0,
        year_stride=int(p.COMPARE_YEAR_STRIDE), sensitivity_set=p.COMPARE_SENSITIVITY_SET,
        electricity_eur_per_kwh=p.ELECTRICITY_EUR_PER_KWH, water_eur_per_m3=p.WATER_EUR_PER_M3 + p.WASTEWATER_EUR_PER_M3,
        copper_eur_per_kg=p.COPPER_EUR_PER_KG, aluminium_eur_per_kg=p.ALUMINIUM_EUR_PER_KG, lifetime_years=p.LIFETIME_YEARS)


# =============================================================================
# Research constants (not settings)
# =============================================================================

# Literature ranges of the one-at-a-time sensitivity study: name -> (lo, hi, unit, source)
RANGES = {
    'wall_mm':    (0.3, 2.0, "mm", "Kayansayan 0.5 (copper) ... heavy industrial tube ~2 mm; ASSUMED where the source gives none"),
    'd_mm':       (6.35, 16.3, "mm", "Wang II Tab. 1 (6.35-13.3, applicable to 12.7) / Kayansayan (9.52-16.3)"),
    's_f_mm':     (1.19, 8.7, "mm", "Wang II conclusion p. 2700 (Kayansayan 2.2-4.2)"),
    'delta_f_mm': (0.1, 0.3, "mm", "Wang I 0.115 / Kayansayan 0.2 (0.1-0.3 bracketing both)"),
    's_t_mm':     (17.7, 40.0, "mm", "Wang II 17.7-31.75 / Kayansayan 25.4-40"),
    's_l_mm':     (12.4, 34.67, "mm", "Wang II 12.4-27.5 / Kayansayan 22-34.67"),
    'N_r':        (1, 6, "-", "Wang II 1-6; Kroger: 3-6 rows typical, Towler: max 10"),
    'w_fr':       (0.5, 3.0, "m/s", "typical dry-cooler range; Kayansayan tunnel 0-15"),
}

# Material densities for the mass estimate [kg/m3] -- ASSUMPTION (typical handbook values;
# cite e.g. VDI Heat Atlas D6 / DIN EN ISO 10456 before use in the thesis)
DENSITY = {'Aluminum': 2700.0, 'Copper': 8960.0, 'Carbon Steel': 7850.0}
