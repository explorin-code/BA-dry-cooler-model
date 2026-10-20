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
  outputs        -- which of src/workflow/outputs.py's catalogue to produce after the
                    scenarios: cooler-results / convergence figures and the
                    caching / resolution / annual benchmarks
"""

import logging

import src.parameters as parameters

from src.inputs.operating_conditions import get_operating_conditions
from src.inputs.geometry import get_geometry
from src.thermo.fluid_properties import get_fluid_properties
from src.solvers import SOLVER_FAMILY, scenario_keys, solve_scenario
from src.inputs.solver_settings import get_solver_settings
from src.thermo.precooling import calc_precooler
from src.auxiliaries.pressure_drop import calc_delta_p_coolant, calc_delta_p_air
from src.auxiliaries.economics import calc_pump_power, calc_fan_power, calc_water_usage
from src.workflow.run_scenario import (format_input_conditions, format_geometry_info,
                              format_output_conditions, format_economics, calc_pinch, calc_deviation)
from src.catalogue import run_items
from src.workflow.outputs import OUTPUTS, RunContext
from src.inputs.run_modes import get_run_modes

logger = logging.getLogger(__name__)


def _log_scenario_summary(label, result, ops, geo, P_pump, P_fan, T_amb, m_dot_ev=None):
    """Always (INFO): one line per solver -- Q, outlet temperatures, pinch, deviation from the reference
    (most detailed) solver, iterations, wall time -- plus power/water.
    T_amb: the AMBIENT air temperature (the pinch is that of the whole cooler, also in the precooled scenario).
    Insight mode (DEBUG): additionally the full input/result text boxes of
    the results figure, so nothing is lost without plots."""
    ref_key, ref = result.reference_key, result.reference
    ref_label = SOLVER_FAMILY[ref_key]
    logger.info("[%s]", label)
    for key, r in result.results.items():
        pinch = calc_pinch(r.T_c_o, T_amb)
        vs_ref = ("" if key == ref_key else
                  f"  (vs. {ref_label}: Q {calc_deviation(r.Q_dot, ref.Q_dot):+5.2f} %, "
                  f"pinch {r.T_c_o - ref.T_c_o:+5.2f} K)")
        logger.info("  %-11s Q = %6.2f kW   T_coolant_out = %5.2f °C   T_air_out = %5.2f °C   "
                    "pinch = %5.2f K   k = %5.2f W/m²K   %3d it   %8.1f ms%s",
                    SOLVER_FAMILY[key], r.Q_dot / 1000, r.T_c_o, r.T_a_o, pinch, r.k,
                    len(r.history_hot), r.solve_time * 1e3, vs_ref)
    logger.info("  %s", format_economics(P_pump, P_fan, m_dot_ev))

    logger.debug("\n%s\n%s", format_input_conditions(ops), format_geometry_info(geo))
    for key, r in result.results.items():
        logger.debug("\n%s", format_output_conditions(SOLVER_FAMILY[key], r, T_amb,
                                                       reference=None if key == ref_key else ref,
                                                       reference_label=ref_label))


def run_scenarios(modes=None):
    """Solves the ambient scenario and, if the decider engages precooling,
    the precooled one, then the selected outputs. Figures are built
    but not shown -- the caller decides when to call plt.show().
    Returns a dict of the solved results for programmatic use."""
    if modes is None:
        modes = get_run_modes()

    # Solver tuning (relaxation factors, resolutions): see parameters.py.
    settings = get_solver_settings()

    geo = get_geometry()
    ops_ambient = get_operating_conditions(geo=geo)

    result_ambient = solve_scenario(ops_ambient, geo, scenario_keys(modes.solvers), settings)

    # Coolant-side pressure drop / pump power -- independent of precooling,
    # since the coolant loop is unaffected by the air-side pad.
    coolant_state = get_fluid_properties(ops_ambient, ops_ambient.T_c_i, ops_ambient.p_c)
    Delta_p_c = calc_delta_p_coolant(geo, ops_ambient, coolant_state)
    P_pump = calc_pump_power(ops_ambient, coolant_state, Delta_p_c)

    # Ambient = dry operation: pad ΔP only if the dry air path runs through the pad.
    T_a_o = result_ambient.reference.T_a_o
    P_fan_ambient = calc_fan_power(ops_ambient, T_a_o, calc_delta_p_air(geo, ops_ambient, T_a_o, parameters.PAD_IN_DRY_AIR_PATH))

    label_ambient = "Dry"
    _log_scenario_summary(label_ambient, result_ambient, ops_ambient, geo, P_pump, P_fan_ambient, ops_ambient.T_a_i)

    results = {'geo': geo, 'ambient': {'label': label_ambient, 'ops': ops_ambient, 'result': result_ambient,
                                       'P_pump': P_pump, 'P_fan': P_fan_ambient, 'm_dot_ev': None},
               'precooled': None}

    ops_final, was_precooled = calc_precooler(ops_ambient, result_ambient, geo)
    if was_precooled:
        m_dot_ev = calc_water_usage(ops_ambient, ops_final)

        result_final = solve_scenario(ops_final, geo, scenario_keys(modes.solvers), settings)

        T_a_o = result_final.reference.T_a_o
        P_fan_final = calc_fan_power(ops_final, T_a_o, calc_delta_p_air(geo, ops_final, T_a_o, through_pad=True))

        _log_scenario_summary("Precooled", result_final, ops_final, geo, P_pump, P_fan_final, ops_ambient.T_a_i, m_dot_ev)
        results['precooled'] = {'label': "Precooled", 'ops': ops_final, 'result': result_final,
                                'P_pump': P_pump, 'P_fan': P_fan_final, 'm_dot_ev': m_dot_ev}

    # The selected outputs (figures, benchmarks -- src/workflow/outputs.py); the benchmarks run their own fresh solves.
    scenarios = [results['ambient']] + ([results['precooled']] if results['precooled'] else [])
    context = RunContext(geo=geo, settings=settings, modes=modes, scenarios=scenarios)
    results['outputs'] = run_items(OUTPUTS, modes.outputs, args=(context,), logger=logger)
    return results
