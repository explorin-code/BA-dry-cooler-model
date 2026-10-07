"""
parameters.py
==============
Raw numerical inputs for a solver run: operating conditions and cooler
geometry. No logic -- see operating_conditions.py and dry_cooler_physics.py
for how these are used/derived. Leave a value as None to fall back to
that field's class-level default instead.
"""

# --- Operating conditions: temperatures [°C] ----------------------------
T_COOLANT_IN = 37.0             # hot coolant entering -- Konrad: 37 °C (target output: see T_COOLANT_TARGET_OUT below)
T_AIR_IN = 20.0                 # cold air entering

# --- Operating conditions: coolant flow -- exactly one nonzero ----------
W_COOLANT = 0.0                 # velocity inside tubes         [m/s]
M_COOLANT = 0.28                # mass flow rate                [kg/s]  -- Konrad: 0.28 kg/s
V_COOLANT = 0.0                 # volumetric flow rate          [m3/s]

# --- Operating conditions: air flow -- exactly one nonzero --------------
W_O = 2.0                       # air approach velocity         [m/s]  -- sweet spot ~2-2.5 m/s
M_O = 0.0                       # air mass flow rate            [kg/s]
V_O = 0.0                       # air volumetric flow rate      [m3/s]

# --- Operating conditions: fallback values -- None uses the class default ---
COOLANT_TYPE = 'Water'          # CoolProp fluid, e.g. 'INCOMP::MEG-30%' (ethylene glycol, 30 % by mass);
                                 # GUI choices: fluid_properties.COOLANTS
P_COOLANT = None                # [Pa]
P_AIR = None                    # [Pa]
PHI_AIR_PERCENT = 30.0          # inlet relative humidity [%] -- give at most one of PHI_AIR_PERCENT/Y_AIR
Y_AIR = None                    # inlet humidity ratio Y [kg water/kg dry air]

# --- Geometry: tube/fin dimensions ----------------------------------------
# Lengths in mm here (and in the GUI); converted to m when read (param_loader.MM).
D_TUBE_OUTER_MM = 9.0           # tube outer diameter           [mm]
FIN_THICKNESS_MM = 0.12         # fin thickness                 [mm]
FIN_PITCH_MM = 1.62             # fin pitch s_f = spacing + thickness [mm]  -- 1.5 mm clear spacing
D_TUBE_INNER_MM = 8.0           # tube inner diameter           [mm]  -- 0.5 mm wall
N_TUBES = 17                    # number of tubes (parallel)
N_ROWS = 6                      # number of tube rows
PITCH_TRANSVERSE_MM = 30.0      # tube pitch s_t, transverse    [mm]
PITCH_LONGITUDINAL_MM = 30.0    # tube pitch s_l, longitudinal  [mm]
HEIGHT_MM = 1000.0              # cooler height == single tube-pass length l [mm]
FIN_COLLAR_THICKNESS_MM = 0.0   # fin-collar thickness delta_c [mm] -- 0 = no collar (d_c = d); only
                                 # enters Wang's d_c (air friction factor, narrowest gap)

# --- Geometry: fallback values -- None uses the class default ------------
PIPE_MATERIAL = 'Copper'        # tube wall material, e.g. 'Aluminum', 'Carbon Steel'
FIN_MATERIAL = 'Aluminum'       # fin material

# --- Run modes: each can be overridden per run by a command-line flag (see main.py) ---
INSIGHT_MODE = False             # per-iteration solver progress + detailed results in the terminal
                                 # (off: one summary block per scenario only)
PLOT_RESULTS = True              # cooler-results figures (inputs, results, profiles) per scenario
PLOT_CONVERGENCE = False         # convergence + iteration-error figure per scenario
BENCHMARK_MODE = True          # caching benchmark + solver-performance figure (adds ~20 s)
RESOLUTION_MODE = True          # resolution sweep + its figure (adds ~15 s)
CELL_2D = False                 # Cell: one representative tube (2D, ~N_TUBES x faster, no neighbour
                                 # averaging of air) instead of all tubes (3D)
BENCHMARK_RESOLUTIONS = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100]
                                 # Cell segments / NTU elements swept by RESOLUTION_MODE
ANNUAL_MODE = False             # annual benchmark: LMTD, NTU (field), NTU (table) and Cell over 8760
                                 # hourly operating points + its figure (several minutes)
ANNUAL_TIME_BUDGET_S = 120.0    # annual benchmark: a solver exceeding this is stopped and extrapolated [s]
ANNUAL_WEATHER_CSV = None       # annual benchmark: CSV path (columns T_air_C, phi_percent; 8760 rows);
                                 # None = synthetic Munich-like year
BENCHMARK_REPEATS = 3           # timed repeats per caching measurement (median is reported)
BENCHMARK_SWEEP_REPEATS = 1     # timed repeats per resolution point (iteration counts are exact anyway)

# --- Solver tuning: shared by LMTD/NTU/Cell -------------------------------
DT_HOT_IT_INIT = 30.0           # initial Delta_T_hot guess                  [K]
DT_COLD_IT_INIT = 30.0          # initial Delta_T_cold guess                 [K]
CENTRAL_OMEGA = 0.3             # under-relaxation factor of LMTD's outer loop (the loop is shared
                                 # with NTU, but NTU uses NTU_OMEGA and Cell CELL_OMEGA, see below);
                                 # stable up to 0.4, diverges from 0.5 (winter -10 °C) / 0.7 (all)
OUTER_MAX_ITER = 500            # iteration cap of LMTD/NTU's outer loop (raises if reached)
CELL_N_SEGMENTS = 20            # coolant-direction segments per tube pass (Cell only)
CELL_OMEGA = 1.0                # Cell's own relaxation factor -- NOT CENTRAL_OMEGA: at 0.2 Cell's
                                 # (then per-cell) stopping criterion fired ~0.02 K before
                                 # convergence (Q +0.17 %); at 1.0 it converges fast at any resolution
NTU_N_ELEMENTS = 20             # elements per tube per row (NTU only) -- independent of CELL_N_SEGMENTS
NTU_MODE = 'field'              # 'field': element field in every outer iteration (internal profile);
                                 # 'table': P(NTU, R) computed once per N_r/N_e, stored in cache/, looked up
                                 # -- for many operating points at constant row count (pays off from ~1000)
NTU_OMEGA = 1.0                 # NTU's own outer-loop (property) relaxation -- NOT CENTRAL_OMEGA: NTU
                                 # converges monotonically, undamped in ~5 instead of 14 iterations

# --- Solver tuning: convergence -------------------------------------------
CONVERGENCE_THRESHOLD = 1e-3    # stopping criterion of all three solvers: both outlet temperatures
                                 # change by less than this per iteration [K] (NTU's inner field:
                                 # 1/10 of it). Keep >= ~1e-4: the property cache rounds T to
                                 # 0.01 K, which leaves Cell a ~1.5e-5 K limit cycle
NTU_FIELD_MAX_ITER = 2000
CELL_MAX_ITER = 1000            # iteration cap of Cell (raises if reached)

# --- Precooling: decision target and dry-operation air path --------------
T_COOLANT_TARGET_OUT = 25.0     # precool if Cell's T_c_o exceeds this -- Konrad's target
PAD_IN_DRY_AIR_PATH = False      # dry operation (no precooling): True = air still passes the (dry) pad,
                                 # so its pressure drop counts; False = separate bypass inlet, no pad ΔP.
                                 # Precooled operation always passes the pad.
