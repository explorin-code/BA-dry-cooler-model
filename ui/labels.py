"""
labels.py -- how the GUI names and orders things: modes, group order per mode, the "exactly one nonzero" flow triples, and
the display labels of the parameters (thesis symbols in LaTeX, rendered by Streamlit).
"""

MODES = {"Solver comparison": "solver", "Cooler comparison": "compare"}

# Sidebar order of the groups per mode (parameters.py's file order is not the order of use); groups not listed follow in file
# order. The open ones are what one changes most often, the rest sits in collapsed expanders.
GROUP_ORDER = {"solver": ("Run", "Operating conditions", "Coolant", "Precooling", "Cooler", "Benchmark settings", "Solver tuning"),
               "compare": ("Analyses", "Comparison", "Coolant", "Precooling", "Prices")}
OPEN_GROUPS = {"solver": ("Run", "Operating conditions", "Cooler"), "compare": ("Analyses", "Comparison")}

# "Exactly one nonzero" triples of parameters.py: shown as ONE choice (what is given) + ONE number.
# (constant, plain name for the selectbox -- its options cannot render LaTeX, symbol label for the number field)
FLOW_TRIPLES = {
    "water flow": (("M_COOLANT", "mass flow [kg/s]", r"$\dot{m}_c$ water mass flow [kg/s]"),
                   ("W_COOLANT", "velocity in the tubes [m/s]", r"$w_c$ water velocity in the tubes [m/s]"),
                   ("V_COOLANT", "volume flow [m³/s]", r"$\dot{V}_c$ water volume flow [m³/s]")),
    "air flow": (("W_O", "face velocity [m/s]", r"$w_{fr}$ air face velocity [m/s]"),
                 ("M_O", "mass flow [kg/s]", r"$\dot{m}_a$ air mass flow [kg/s]"),
                 ("V_O", "volume flow [m³/s]", r"$\dot{V}_a$ air volume flow [m³/s]")),
}
FLOW_MEMBERS = {name: title for title, members in FLOW_TRIPLES.items() for name, *_ in members}

# Display labels: parameters.py constant -> thesis symbol (LaTeX, rendered by Streamlit) + short text + unit.
# Constants without a symbol get a plain-text label. The constant's own name stays visible in the help text.
LABELS = {
    # both
    'T_COOLANT_IN': r"$T_{c,i}$ coolant inlet [°C]", 'T_COOLANT_TARGET_OUT': r"$T_{c,o}$ target coolant outlet [°C]",
    'PAD_IN_DRY_AIR_PATH': r"pad in the dry air path ($\Delta p_{pad}$ counts)",
    'ELECTRICITY_EUR_PER_KWH': "electricity price [€/kWh]", 'WATER_EUR_PER_M3': "water price [€/m³]",
    'WASTEWATER_EUR_PER_M3': "wastewater fee [€/m³]",
    'COPPER_EUR_PER_KG': "copper price [€/kg]", 'ALUMINIUM_EUR_PER_KG': "aluminium price [€/kg]",
    'LIFETIME_YEARS': "service life for the material write-off [years]",
    # solver: operating conditions
    'T_AIR_IN': r"$T_{a,i}$ air inlet [°C]",
    'W_COOLANT': r"$w_c$ water velocity in the tubes [m/s]", 'M_COOLANT': r"$\dot{m}_c$ water mass flow [kg/s]",
    'V_COOLANT': r"$\dot{V}_c$ water volume flow [m³/s]",
    'W_O': r"$w_{fr}$ air face velocity [m/s]", 'M_O': r"$\dot{m}_a$ air mass flow [kg/s]",
    'V_O': r"$\dot{V}_a$ air volume flow [m³/s]",
    'COOLANT_TYPE': "coolant fluid", 'P_COOLANT': r"$p_c$ coolant pressure [Pa]", 'P_AIR': r"$p_a$ air pressure [Pa]",
    'PHI_AIR_PERCENT': r"$\varphi$ air relative humidity [%]", 'Y_AIR': r"$Y$ air humidity ratio [kg/kg]",
    # solver: cooler
    'D_TUBE_OUTER_MM': r"$d$ tube outer diameter [mm]", 'D_TUBE_INNER_MM': r"$d_i$ tube inner diameter [mm]",
    'FIN_THICKNESS_MM': r"$\delta_f$ fin thickness [mm]", 'FIN_PITCH_MM': r"$s_f$ fin pitch [mm]",
    'PITCH_TRANSVERSE_MM': r"$s_t$ transverse tube pitch [mm]", 'PITCH_LONGITUDINAL_MM': r"$s_l$ longitudinal tube pitch [mm]",
    'FIN_COLLAR_THICKNESS_MM': r"$\delta_c$ fin collar thickness [mm]",
    'PIPE_MATERIAL': r"tube material ($\lambda_p$)", 'FIN_MATERIAL': r"fin material ($\lambda_f$)",
    'N_TUBES': r"$N_t$ tubes per row", 'N_ROWS': r"$N_r$ tube rows", 'HEIGHT_MM': r"$l$ tube length = cooler height [mm]",
    # solver: run
    'SOLVERS': "solvers that run", 'INSIGHT_MODE': "insight mode (per-iteration progress)", 'SOLVER_OUTPUTS': "figures and benchmarks",
    'BENCHMARK_RESOLUTIONS': r"resolutions swept ($N_s$ / $N_e$)",
    'ANNUAL_CELL_STRIDE': "annual benchmark: Cell on every n-th hour",
    'ANNUAL_TIME_BUDGET_S': "annual benchmark: time budget per solver [s]", 'ANNUAL_WEATHER_CSV': "annual benchmark: own weather CSV",
    'BENCHMARK_REPEATS': "caching benchmark: repeats", 'BENCHMARK_SWEEP_REPEATS': "resolution benchmark: repeats",
    # solver: solver tuning
    'DT_HOT_IT_INIT': r"$\Delta T_{hot}$ initial guess [K]", 'DT_COLD_IT_INIT': r"$\Delta T_{cold}$ initial guess [K]",
    'CENTRAL_OMEGA': r"$\omega_{LMTD}$ relaxation", 'OUTER_MAX_ITER': "outer loop: maximum iterations",
    'CELL_N_SEGMENTS': r"$N_s$ Cell segments per tube pass", 'CELL_OMEGA': r"$\omega_{Cell}$ relaxation",
    'NTU_N_ELEMENTS': r"$N_e$ NTU elements per tube and row",
    'NTU_OMEGA': r"$\omega_{NTU}$ relaxation", 'CONVERGENCE_THRESHOLD': "convergence threshold, outlets [K]",
    'NTU_FIELD_MAX_ITER': "NTU field: maximum sweeps", 'CELL_MAX_ITER': "Cell: maximum iterations",
    # compare
    'COMPARE_SETS': "coolers compared", 'COMPARE_SOLVER': "solver", 'COMPARE_SENSITIVITY_SET': "cooler of the sensitivity study",
    'COMPARE_M_COOLANT': r"$\dot{m}_c$ water mass flow per module [kg/s]", 'COMPARE_RE_TUBE': r"$Re_c$ tube Reynolds number",
    'COMPARE_N_ROWS': r"$N_r$ tube rows", 'COMPARE_TUBE_LENGTH_MM': r"$l$ tube length [mm]",
    'COMPARE_DESIGN_T_AIR': r"design $T_{a,i}$ [°C]", 'COMPARE_DESIGN_PHI_PERCENT': r"design $\varphi$ [%]",
    'COMPARE_W_FR_DESIGN': r"$w_{fr,design}$ nominal face velocity [m/s]", 'COMPARE_W_FR_MIN': r"$w_{fr,min}$ slowest fan setting [m/s]",
'COMPARE_T_FREEZE_RISK': r"$T_{c,o}$ freezing-risk limit [°C]",
    'COMPARE_N_TUBES': r"$N_t$ tubes per row (0 = each cooler's own)", 'COMPARE_ANALYSES': "analyses", 'COMPARE_YEAR_STRIDE': "year analyses: every n-th hour",
}


def label(name):
    return LABELS.get(name, name)
