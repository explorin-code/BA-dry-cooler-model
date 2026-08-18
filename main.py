"""
main.py
=======
Runs the dry-cooler scenario twice: once on ambient conditions, once with
precooling.calc_precooler() applied -- only if precooling actually engages.
The precooling decision reuses the ambient scenario's already-solved Cell
result rather than re-solving. Both windows are shown together.
"""

import matplotlib.pyplot as plt

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import solve_scenario
from src.precooling import calc_precooler
from src.pressure_drop import calc_delta_p_total, calc_delta_p_bundle_air, calc_delta_p_air_total
from src.economics import calc_pump_power, calc_fan_power, calc_water_usage
from src.run_scenario import plot_scenario
from src.plot_profiles import plot_profiles

# ------------------------------------------------------------------
# Centralized under-relaxation factor / cell resolution.
# All three solvers get the SAME omega so the step size is directly
# comparable between them.
# ------------------------------------------------------------------
CENTRAL_OMEGA = 0.2
CELL_N_SEGMENTS = 20

ops_ambient = get_operating_conditions()
geo = get_geometry()

result_ambient = solve_scenario(ops_ambient, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS)

# Coolant-side pressure drop / pump power -- independent of precooling,
# since the coolant loop is unaffected by the air-side pad.
coolant_state = get_fluid_properties(ops_ambient, ops_ambient.T_coolant_in, ops_ambient.P_coolant)
delta_p_coolant = calc_delta_p_total(geo, ops_ambient, coolant_state)
W_pump = calc_pump_power(ops_ambient, coolant_state, delta_p_coolant)

# Air-side pressure drop / fan power -- depends on ops, so recomputed per scenario.
air_state_ambient = get_air_properties(ops_ambient, ops_ambient.T_air_in, ops_ambient.P_air)
delta_p_bundle_ambient = calc_delta_p_bundle_air(geo, ops_ambient, air_state_ambient)
delta_p_air_ambient = calc_delta_p_air_total(delta_p_bundle_ambient)
W_fan_ambient = calc_fan_power(ops_ambient, air_state_ambient, delta_p_air_ambient)

plot_scenario(result_ambient, ops_ambient, geo, label="Ambient (no precooling)",
              omega=CENTRAL_OMEGA, W_pump=W_pump, W_fan=W_fan_ambient)
plot_profiles(result_ambient, ops_ambient, geo, label="Ambient (no precooling)", n_segments=CELL_N_SEGMENTS)

ops_final, was_precooled = calc_precooler(ops_ambient, result_ambient)
if was_precooled:
    m_water = calc_water_usage(ops_ambient, ops_final)

    result_final = solve_scenario(ops_final, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS)

    air_state_final = get_air_properties(ops_final, ops_final.T_air_in, ops_final.P_air)
    delta_p_bundle_final = calc_delta_p_bundle_air(geo, ops_final, air_state_final)
    delta_p_air_final = calc_delta_p_air_total(delta_p_bundle_final)  # pad ΔP not yet modeled
    W_fan_final = calc_fan_power(ops_final, air_state_final, delta_p_air_final)

    plot_scenario(result_final, ops_final, geo, label="Precooled",
                  omega=CENTRAL_OMEGA, W_pump=W_pump, W_fan=W_fan_final, m_water=m_water)
    plot_profiles(result_final, ops_final, geo, label="Precooled", n_segments=CELL_N_SEGMENTS)

plt.show()
