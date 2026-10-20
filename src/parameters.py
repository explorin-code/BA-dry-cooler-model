"""
parameters.py
==============
Raw numerical inputs for a run. No logic -- see operating_conditions.py, dry_cooler_physics.py and
solver_settings.py (solver comparison) and comparison/config.py (cooler comparison) for how they are used.

The file is also the source of the GUI (gui.py builds its form from it): every section header
    # --- Name [mode]: sub-heading ---
puts the constants below it into the group "Name", shown in the GUI mode(s) in [brackets]:
    [solver]   solver comparison: one cooler, the selected solvers, ambient / precooled scenario, solver benchmarks
    [compare]  cooler comparison: the cooler sets raced against each other over the climate year
    [both]     used by both
Inline comments become the help texts. Leave a value as None to fall back to that field's default.
"""

# =====================================================================================================
# Both workflows
# =====================================================================================================

# --- Coolant [both]: temperatures [°C] -------------------------------------
T_COOLANT_IN = 37.0             # hot coolant entering
T_COOLANT_TARGET_OUT = 25.0     # target coolant outlet; precool / speed up the fan if exceeded

# --- Precooling [both]: dry-operation air path ----------------------------
PAD_IN_DRY_AIR_PATH = False     # dry operation: True = the air passes the dry pad (its pressure drop counts), False = bypass

# --- Prices [compare]: operating cost ---------------------------------------
ELECTRICITY_EUR_PER_KWH = 0.2726  # electricity [EUR/kWh] net: SWM M-Strom single tariff from 15 Sep 2026 (basic charge not included)
WATER_EUR_PER_M3 = 1.84         # Munich drinking water, net, from 2026-07-01 incl. "Wassercent"
                                 # (SWM press release 3.6.2026); 1.97 incl. 7 % VAT
WASTEWATER_EUR_PER_M3 = 0.0     # evaporated pad water normally does not reach the sewer; if charged, ~2 EUR/m3
COPPER_EUR_PER_KG = 13.10465    # copper [EUR/kg]: LME cash 9 Oct 2026 14,689 USD/t at 1.1209 USD/EUR (Westmetall); raw metal, no fabrication premium
ALUMINIUM_EUR_PER_KG = 2.72593  # aluminium [EUR/kg]: LME cash 9 Oct 2026 3,055.5 USD/t at 1.1209 USD/EUR (Westmetall); raw metal
LIFETIME_YEARS = 20             # service life of the straight-line material write-off [years]

# =====================================================================================================
# Solver comparison: one cooler, the selected solvers LMTD / NTU / Cell, solver benchmarks
# =====================================================================================================

# --- Operating conditions [solver]: air temperature [°C] -----------------
T_AIR_IN = 20.0                 # cold air entering

# --- Operating conditions [solver]: coolant flow -- exactly one nonzero --
W_COOLANT = 0.0                 # velocity inside tubes         [m/s]
M_COOLANT = 1.0                 # mass flow rate                [kg/s]  (the thesis' primary operating case: 0.28 kg/s)
V_COOLANT = 0.0                 # volumetric flow rate          [m3/s]

# --- Operating conditions [solver]: air flow -- exactly one nonzero ------
W_O = 2.5                       # air approach velocity [m/s] -- a normal speed; the design / top speed is 3 m/s
M_O = 0.0                       # air mass flow rate            [kg/s]
V_O = 0.0                       # air volumetric flow rate      [m3/s]

# --- Operating conditions [solver]: fallback values -- None uses the class default ---
COOLANT_TYPE = 'Water'          # CoolProp fluid, e.g. 'INCOMP::MEG-30%' (ethylene glycol, 30 % by mass);
                                 # GUI choices: fluid_properties.COOLANTS
P_COOLANT = None                # [Pa]
P_AIR = None                    # [Pa]
PHI_AIR_PERCENT = 60.0          # inlet relative humidity [%] -- give at most one of PHI_AIR_PERCENT/Y_AIR
Y_AIR = None                    # inlet humidity ratio Y [kg water/kg dry air]

# --- Cooler [solver]: elementary cell -------------------------------------
# GEOMETRY_SET picks the standard cell (src/cooler_sets.py); any value below that is not None replaces the
# preset's value for this run (in the GUI the fields show the preset and can be changed). Lengths in mm.
GEOMETRY_SET = 'P2'             # 'P1' dense (Wang) | 'P2' medium (Kayansayan 5) | 'P3' open (Kayansayan 3)
D_TUBE_OUTER_MM = None          # tube outer diameter d                      [mm]
D_TUBE_INNER_MM = None          # tube inner diameter d_i (= d - 2 wall)     [mm]
FIN_THICKNESS_MM = None         # fin thickness delta_f                      [mm]
FIN_PITCH_MM = None             # fin pitch s_f = spacing + thickness        [mm]
PITCH_TRANSVERSE_MM = None      # tube pitch s_t, transverse                 [mm]
PITCH_LONGITUDINAL_MM = None    # tube pitch s_l, longitudinal               [mm]
FIN_COLLAR_THICKNESS_MM = None  # fin-collar thickness delta_c [mm]; only enters Wang's d_c (air friction factor)
PIPE_MATERIAL = None            # tube wall material: 'Copper', 'Aluminum', 'Carbon Steel'
FIN_MATERIAL = None             # fin material

# --- Cooler [solver]: module (not part of the preset) ---------------------
N_TUBES = None                  # tubes per row; None = the cooler's fixed reference count (does not follow the inputs)
N_ROWS = 6                      # tube rows (= passes in series)
HEIGHT_MM = 2000.0              # cooler height == single tube-pass length l [mm] (same as COMPARE_TUBE_LENGTH_MM)

# --- Run [solver]: solvers and outputs (command-line flags override, see main.py) ---
SOLVERS = ['lmtd', 'ntu_field', 'cell_3d']   # solvers that run; the most detailed is the reference of the deviations
INSIGHT_MODE = False             # per-iteration solver progress + detailed results in the terminal
                                 # (off: one summary block per scenario only)
SOLVER_OUTPUTS = ['cooler_results', 'caching', 'resolution']   # figures and benchmarks (annual = minutes)

# --- Output settings [solver]: only used by the optional outputs and benchmarks ---
BENCHMARK_RESOLUTIONS = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100]
                                 # Cell segments / NTU elements swept by the resolution benchmark
BENCHMARK_TABLE_MAX_N = 100       # resolution benchmark: NTU tables are built and timed up to this N (the build grows steeply: ~13 s at N = 30)
ANNUAL_CELL_STRIDE = 4           # annual benchmark: Cell runs on every n-th hour only (its year total is extrapolated)
ANNUAL_TIME_BUDGET_S = 120.0     # annual benchmark: a solver exceeding this is stopped and extrapolated [s]; 0 = no limit
ANNUAL_WEATHER_CSV = None        # annual benchmark: weather file; None = the past year
PROFILE_T_COOLANT_HOT = 90.0     # profile grids: second water inlet temperature, next to the run's own [°C]
ACCURACY_T_COOLANT = [30, 40, 50, 60, 70, 80, 90]   # accuracy sweep: water inlet temperatures [°C]
BENCHMARK_REPEATS = 3            # timed repeats per caching measurement (median is reported)
BENCHMARK_SWEEP_REPEATS = 3      # timed repeats per resolution point (iteration counts are exact anyway)

# --- Solver tuning [solver]: shared by LMTD/NTU/Cell -----------------------
DT_HOT_IT_INIT = 30.0           # initial Delta_T_hot guess                  [K]
DT_COLD_IT_INIT = 30.0          # initial Delta_T_cold guess                 [K]
CENTRAL_OMEGA = 0.3             # relaxation of LMTD's outer loop (stable up to 0.4)
OUTER_MAX_ITER = 500            # iteration cap of LMTD/NTU's outer loop (raises if reached)
CELL_N_SEGMENTS = 20            # coolant-direction segments per tube pass (Cell only)
CELL_OMEGA = 1.0                # relaxation of Cell (1.0 converges fast at any resolution)
NTU_N_ELEMENTS = 20             # elements per tube per row (NTU field and table) -- independent of CELL_N_SEGMENTS
NTU_OMEGA = 1.0                 # relaxation of NTU's outer loop (converges monotonically)

# --- Solver tuning [solver]: convergence ----------------------------------
CONVERGENCE_THRESHOLD = 1e-3    # stop when both outlet temperatures change less than this [K]; keep >= 1e-4
NTU_FIELD_MAX_ITER = 2000
CELL_MAX_ITER = 1000            # iteration cap of Cell (raises if reached)

# =====================================================================================================
# Cooler comparison: the standard coolers raced against each other over the climate year
# =====================================================================================================

# --- Comparison [compare]: coolers and module -----------------------------
COMPARE_SETS = ['P1', 'P2', 'P3']   # the coolers compared (src/cooler_sets.py)
COMPARE_SOLVER = 'ntu_table'    # solver of every analysis (Cell is ~100x slower: use the year stride)
COMPARE_M_COOLANT = 1.0         # water mass flow per module [kg/s]; more flow = more tubes of the same kind
COMPARE_RE_TUBE = 5000          # reference tube Reynolds number of the tubes / wall studies
COMPARE_N_ROWS = 6              # tube rows of every module (Wang's tested limit; Kroger 3-6 typical)
COMPARE_N_TUBES = None          # tubes per row; None = each cooler's fixed reference count
COMPARE_TUBE_LENGTH_MM = 2000.0 # tube length (module height) of every module [mm] -- changes the behaviour, not only the size
COMPARE_SENSITIVITY_SET = 'P2'  # the cell the sensitivity study varies

# --- Comparison [compare]: operation --------------------------------------
COMPARE_DESIGN_T_AIR = 20.0     # design ambient for the single-point analyses (geometry, rows, tubes, wall) [°C]
COMPARE_DESIGN_PHI_PERCENT = 60.0   # design ambient relative humidity [%]
COMPARE_W_FR_DESIGN = 3.0       # fan speed [m/s]: the fixed speed and the top of the variable range
COMPARE_W_FR_MIN = 0.4          # slowest fan setting [m/s]
COMPARE_T_FREEZE_RISK = 2.0     # coolant outlet below this = freezing risk [°C] (the tube wall is colder than the bulk)

# --- Analyses [compare]: what to compute and plot (each writes CSVs and figures, see comparison/analyses.py) ---
COMPARE_ANALYSES = ['geometry']  # analyses to run (the year with fan control takes 2-3 min)
COMPARE_YEAR_STRIDE = 1         # year analyses: use every n-th hour (1 = all 8760 h; 24 = quick test)
