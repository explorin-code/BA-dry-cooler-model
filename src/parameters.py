"""
parameters.py
==============
Raw numerical inputs for a solver run: operating conditions and cooler
geometry. No logic -- see operating_conditions.py and dry_cooler_physics.py
for how these are used/derived. Leave a value as None to fall back to
that field's class-level default instead.
"""

# --- Operating conditions: temperatures [°C] ----------------------------
T_COOLANT_IN = 37.0             # hot coolant entering -- Konrad: 37 °C, 25 °C target output
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
PHI_AIR = None                  # inlet relative humidity [0-1] -- give at most one of PHI_AIR/X_AIR
X_AIR = None                    # inlet humidity ratio [kg water/kg dry air]

# --- Geometry: tube/fin dimensions ----------------------------------------
D_TUBE_OUTER = 0.009            # tube outer diameter           [m]
FIN_THICKNESS = 0.00012         # fin thickness                 [m]
FIN_SPACING = 0.00023            # fin spacing                   [m]
D_TUBE_INNER = 0.008            # tube inner diameter           [m]  -- 12 mm OD, 0.5 mm wall
N_TUBES = 17                    # number of tubes (parallel)
N_ROWS = 6                      # number of tube rows
T_Q = 0.03                      # tube pitch, transverse        [m]
T_L = 0.03                      # tube pitch, longitudinal      [m]
HEIGHT = 1                    # cooler height == single tube-pass length [m]

# --- Geometry: fallback values -- None uses the class default ------------
MATERIAL = None                 # e.g. 'Aluminum'