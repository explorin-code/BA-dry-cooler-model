"""
scenario_pipeline.py
======================
Assembles and runs the dry-cooler scenario(s): solve ambient conditions,
compute pump/fan power, plot; then -- only if precooling.calc_precooler()
decides it's needed, reusing the ambient scenario's already-solved Cell
result rather than re-solving -- do the same for the precooled scenario.

This module exists so main.py doesn't have to import every subsystem
directly: main.py just calls run_scenarios() and decides when to show the
resulting figures.
"""

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import solve_scenario
from src.solver_settings import get_solver_settings
from src.precooling import calc_precooler
from src.pressure_drop import calc_delta_p_coolant, calc_delta_p_bundle_air, calc_delta_p_pad, calc_delta_p_air_total
from src.economics import calc_pump_power, calc_fan_power, calc_water_usage
from src.run_scenario import plot_scenario
from src.plot_profiles import plot_profiles


def _solve_air_side_power(geo, ops, T_air_out):
    """Air-side pressure drop -> fan power for one scenario (ambient or
    precooled) -- depends on ops, so called once per scenario. T_air_out
    comes from the Cell result, the same reference the precooling decider
    uses. The pad is physically installed in both scenarios, so its ΔP
    always counts."""
    air_in = get_air_properties(ops, ops.T_air_in, ops.P_air)
    air_out = get_air_properties(ops, T_air_out, ops.P_air)
    air_mean = get_air_properties(ops, (ops.T_air_in + T_air_out) / 2.0, ops.P_air)

    delta_p_bundle = calc_delta_p_bundle_air(geo, ops, air_in, air_out, air_mean)
    delta_p_air = calc_delta_p_air_total(delta_p_bundle, calc_delta_p_pad())
    return calc_fan_power(ops, air_in, air_out, delta_p_air)


def run_scenarios():
    """Builds and returns nothing -- solving and plotting are both
    side-effecting (plot_scenario/plot_profiles build matplotlib figures
    but don't show them). Caller decides when to call plt.show()."""
    # Central under-relaxation factor / cell resolution -- see parameters.py
    # (CENTRAL_OMEGA, CELL_N_SEGMENTS). All three solvers get the SAME omega
    # so the step size is directly comparable between them.
    settings = get_solver_settings()
    CENTRAL_OMEGA = settings.central_omega
    CELL_N_SEGMENTS = settings.cell_n_segments

    geo = get_geometry()
    ops_ambient = get_operating_conditions(geo=geo)

    result_ambient = solve_scenario(ops_ambient, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS)

    # Coolant-side pressure drop / pump power -- independent of precooling,
    # since the coolant loop is unaffected by the air-side pad.
    coolant_state = get_fluid_properties(ops_ambient, ops_ambient.T_coolant_in, ops_ambient.P_coolant)
    delta_p_coolant = calc_delta_p_coolant(geo, ops_ambient, coolant_state)
    P_p = calc_pump_power(ops_ambient, coolant_state, delta_p_coolant)

    P_f_ambient = _solve_air_side_power(geo, ops_ambient, result_ambient.cell.T_air_out)

    plot_scenario(result_ambient, ops_ambient, geo, label="Ambient (no precooling)",
                  omega=CENTRAL_OMEGA, P_p=P_p, P_f=P_f_ambient)
    plot_profiles(result_ambient, ops_ambient, geo, label="Ambient (no precooling)", n_segments=CELL_N_SEGMENTS)

    ops_final, was_precooled = calc_precooler(ops_ambient, result_ambient)
    if was_precooled:
        m_dot_w = calc_water_usage(ops_ambient, ops_final)

        result_final = solve_scenario(ops_final, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS)

        P_f_final = _solve_air_side_power(geo, ops_final, result_final.cell.T_air_out)

        plot_scenario(result_final, ops_final, geo, label="Precooled",
                      omega=CENTRAL_OMEGA, P_p=P_p, P_f=P_f_final, m_dot_w=m_dot_w)
        plot_profiles(result_final, ops_final, geo, label="Precooled", n_segments=CELL_N_SEGMENTS)
