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

Run modes (run_modes.py, set in parameters.py or by main.py's flags):
  insight_mode   -- per-iteration progress + detailed results (DEBUG log
                    level); a short summary per scenario is always logged
  plot_results     -- per scenario: cooler-results figure (inputs, economics,
                      per-solver results next to their spatial profiles)
  plot_convergence -- per scenario: convergence + iteration-error figure
  benchmark_mode   -- caching benchmark (ambient operating point) +
                      solver-performance figure
  resolution_mode  -- resolution sweep (ambient operating point) + figure
"""

import dataclasses
import logging

import src.parameters as parameters

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import solve_scenario
from src.solver_settings import get_solver_settings
from src.precooling import calc_precooler
from src.pressure_drop import calc_delta_p_coolant, calc_delta_p_bundle_air, calc_delta_p_pad, calc_delta_p_air_total
from src.economics import calc_pump_power, calc_fan_power, calc_water_usage
from src.run_scenario import (format_input_conditions, format_geometry_info,
                              format_output_conditions, format_economics, calc_pinch, calc_deviation)
from src.run_modes import get_run_modes

logger = logging.getLogger(__name__)


def _solve_air_side_power(geo, ops, theta_a_o, through_pad: bool):
    """Air-side pressure drop -> fan power for one scenario (ambient or
    precooled) -- depends on ops, so called once per scenario. theta_a_o
    comes from the Cell result, the same reference the precooling decider
    uses. through_pad: whether the air passes the adiabatic pad (always when
    precooling; in dry operation per parameters.PAD_IN_DRY_AIR_PATH, since
    some designs have a separate bypass inlet)."""
    air_in = get_air_properties(ops, ops.theta_a_i, ops.p_a)
    air_out = get_air_properties(ops, theta_a_o, ops.p_a)
    air_mean = get_air_properties(ops, (ops.theta_a_i + theta_a_o) / 2.0, ops.p_a)

    Delta_P_bundle = calc_delta_p_bundle_air(geo, ops, air_in, air_out, air_mean)
    Delta_P_a = calc_delta_p_air_total(Delta_P_bundle, calc_delta_p_pad() if through_pad else None)
    return calc_fan_power(ops, air_in, air_out, Delta_P_a)


def _log_scenario_summary(label, result, ops, geo, W_pump, W_fan, m_dot_ev=None):
    """Always (INFO): one line per solver -- Q, outlet temperatures, pinch,
    deviation from Cell, iterations, wall time -- plus power/water.
    Insight mode (DEBUG): additionally the full input/result text boxes of
    the results figure, so nothing is lost without plots."""
    cell = result.cell
    pinch_cell = calc_pinch(cell.theta_c_o, ops.theta_a_i)
    logger.info("[%s]", label)
    for name, r in (("LMTD", result.lmtd), ("NTU", result.ntu), ("Cell", cell)):
        pinch = calc_pinch(r.theta_c_o, ops.theta_a_i)
        vs_cell = ("" if r is cell else
                   f"  (vs. Cell: Q {calc_deviation(r.Q_dot, cell.Q_dot):+5.2f} %, "
                   f"pinch {calc_deviation(pinch, pinch_cell):+5.2f} %)")
        logger.info("  %-4s Q = %6.2f kW   T_coolant_out = %5.2f °C   T_air_out = %5.2f °C   "
                    "pinch = %5.2f K   k = %5.2f W/m²K   %3d it   %8.1f ms%s",
                    name, r.Q_dot / 1000, r.theta_c_o, r.theta_a_o, pinch, r.k,
                    len(r.history_hot), r.solve_time * 1e3, vs_cell)
    logger.info("  %s", format_economics(W_pump, W_fan, m_dot_ev))

    logger.debug("\n%s\n%s", format_input_conditions(ops), format_geometry_info(geo))
    for name, r in (("LMTD", result.lmtd), ("NTU", result.ntu), ("Cell", cell)):
        logger.debug("\n%s", format_output_conditions(name, r, ops.theta_a_i,
                                                       reference=None if r is cell else cell))


def run_scenarios(modes=None):
    """Solves the ambient scenario and, if the decider engages precooling,
    the precooled one. Figures (plot_results / benchmark_mode) are built
    but not shown -- the caller decides when to call plt.show().
    Returns a dict of the solved results for programmatic use."""
    if modes is None:
        modes = get_run_modes()

    # Central under-relaxation factor / Cell and NTU resolution -- see parameters.py
    # (CENTRAL_OMEGA, CELL_N_SEGMENTS, NTU_N_ELEMENTS). LMTD and NTU share CENTRAL_OMEGA;
    # Cell uses its own CELL_OMEGA (see solve_it_cell for why).
    settings = dataclasses.replace(get_solver_settings(), cell_2d=modes.cell_2d)
    CENTRAL_OMEGA = settings.central_omega
    CELL_N_SEGMENTS = settings.cell_n_segments
    NTU_N_ELEMENTS = settings.ntu_n_elements

    geo = get_geometry()
    ops_ambient = get_operating_conditions(geo=geo)

    result_ambient = solve_scenario(ops_ambient, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS, n_elements=NTU_N_ELEMENTS,
                                    settings=settings)

    # Coolant-side pressure drop / pump power -- independent of precooling,
    # since the coolant loop is unaffected by the air-side pad.
    coolant_state = get_fluid_properties(ops_ambient, ops_ambient.theta_c_i, ops_ambient.p_c)
    Delta_P_c = calc_delta_p_coolant(geo, ops_ambient, coolant_state)
    W_pump = calc_pump_power(ops_ambient, coolant_state, Delta_P_c)

    # Ambient = dry operation: pad ΔP only if the dry air path runs through the pad.
    W_fan_ambient = _solve_air_side_power(geo, ops_ambient, result_ambient.cell.theta_a_o,
                                        through_pad=parameters.PAD_IN_DRY_AIR_PATH)

    label_ambient = "Ambient (no precooling)"
    _log_scenario_summary(label_ambient, result_ambient, ops_ambient, geo, W_pump, W_fan_ambient)

    results = {'geo': geo, 'ambient': {'label': label_ambient, 'ops': ops_ambient, 'result': result_ambient,
                                       'W_pump': W_pump, 'W_fan': W_fan_ambient, 'm_dot_ev': None},
               'precooled': None}

    ops_final, was_precooled = calc_precooler(ops_ambient, result_ambient)
    if was_precooled:
        m_dot_ev = calc_water_usage(ops_ambient, ops_final)

        result_final = solve_scenario(ops_final, geo, omega=CENTRAL_OMEGA, n_segments=CELL_N_SEGMENTS, n_elements=NTU_N_ELEMENTS,
                                    settings=settings)

        W_fan_final = _solve_air_side_power(geo, ops_final, result_final.cell.theta_a_o, through_pad=True)

        _log_scenario_summary("Precooled", result_final, ops_final, geo, W_pump, W_fan_final, m_dot_ev)
        results['precooled'] = {'label': "Precooled", 'ops': ops_final, 'result': result_final,
                                'W_pump': W_pump, 'W_fan': W_fan_final, 'm_dot_ev': m_dot_ev}

    # Benchmarks last, on their own fresh runs (ambient operating point) --
    # imported lazily so normal runs don't pay for them.
    if modes.benchmark_mode or modes.resolution_mode:
        results['benchmark'] = {}
    if modes.benchmark_mode:
        from src.benchmark_solvers import run_cache_benchmark
        results['benchmark']['cache'] = run_cache_benchmark(ops_ambient, geo, settings, modes)
    if modes.resolution_mode:
        from src.benchmark_solvers import run_resolution_benchmark
        results['benchmark']['resolution'] = run_resolution_benchmark(ops_ambient, geo, settings, modes)

    # --- Figures (built but not shown; the caller calls plt.show()) ----------
    scenarios = [results['ambient']] + ([results['precooled']] if results['precooled'] else [])
    if modes.plot_results:
        from src.plot_profiles import plot_cooler_results
        for sc in scenarios:
            plot_cooler_results(sc['result'], sc['ops'], geo, sc['label'], n_segments=CELL_N_SEGMENTS,
                                n_elements=NTU_N_ELEMENTS, W_pump=sc['W_pump'], W_fan=sc['W_fan'], m_dot_ev=sc['m_dot_ev'])
    if modes.plot_convergence:
        from src.plot_solver import plot_convergence
        for sc in scenarios:
            plot_convergence(sc['result'], sc['ops'], sc['label'], settings)
    if modes.benchmark_mode:
        from src.plot_solver import plot_solver_performance
        plot_solver_performance(label_ambient, settings, results['benchmark']['cache'])
    if modes.resolution_mode:
        from src.benchmark_solvers import plot_resolution
        plot_resolution(results['benchmark']['resolution'], label_ambient)

    return results
