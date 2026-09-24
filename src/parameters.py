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
COOLANT_TYPE = None             # e.g. 'Water'
P_COOLANT = None                # [Pa]
P_AIR = None                    # [Pa]
PHI_AIR = 0.3                  # inlet relative humidity [0-1] -- give at most one of PHI_AIR/X_AIR
X_AIR = None                    # inlet humidity ratio [kg water/kg dry air]

# --- Geometry: tube/fin dimensions ----------------------------------------
D_TUBE_OUTER = 0.009            # tube outer diameter           [m]
FIN_THICKNESS = 0.00012         # fin thickness                 [m]
FIN_SPACING = 0.0015            # fin spacing                   [m]
D_TUBE_INNER = 0.008            # tube inner diameter           [m]  -- 12 mm OD, 0.5 mm wall
N_TUBES = 17                    # number of tubes (parallel)
N_ROWS = 6                      # number of tube rows
S_1 = 0.03                      # tube pitch, transverse        [m]
S_2 = 0.03                      # tube pitch, longitudinal      [m]
HEIGHT = 1                    # cooler height == single tube-pass length [m]

# --- Geometry: fallback values -- None uses the class default ------------
PIPE_MATERIAL = 'Copper'        # tube wall material, e.g. 'Aluminum', 'Carbon Steel'
FIN_MATERIAL = 'Aluminum'       # fin material

# --- Solver tuning: shared by LMTD/NTU/Cell -------------------------------
DT_HOT_IT_INIT = 30             # initial dT_hot guess                  [K]
DT_COLD_IT_INIT = 30            # initial dT_cold guess                 [K]
DT_KICK = 5.0                   # forced initial dT_hot/dT_cold to kickstart the while loop [K]
CONVERGENCE_THRESHOLD = 1e-3    # convergence threshold on dT_hot/dT_cold change [K]
CENTRAL_OMEGA = 0.2             # under-relaxation factor -- same for all three solvers,
                                 # so their step sizes are directly comparable
CELL_N_SEGMENTS = 20            # coolant-direction segments per tube pass (Cell only)

# --- Solver tuning: Cell's two-stage relaxation (see solvers._relax_cell_grid) ---
CELL_STAGE1_THRESHOLD = 1e-1    # stage 1: coarse, fast propagation at a large omega
CELL_STAGE1_MAX_ITER = 50
CELL_STAGE1_MIN_ITER = 3
CELL_STAGE2_THRESHOLD = 1e-3    # stage 2: fine polish at the requested omega
CELL_STAGE2_MAX_ITER = 1000
CELL_STAGE2_MIN_ITER = 3

# --- Precooling: decision target ------------------------------------------
T_COOLANT_TARGET_OUT = 25.0     # precool if Cell's T_coolant_out exceeds this -- Konrad's target