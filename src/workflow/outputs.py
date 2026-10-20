"""
outputs.py
===========
The catalogue of what the solver comparison can produce besides its terminal summary (parameters.SOLVER_OUTPUTS,
`main.py --outputs`, the GUI's multiselect). Each output is a catalogue.Item whose steps get the run context:

  cooler_results  the cooler-results figure per scenario (dry / precooled), identical axes
  profile_grids   profiles at the run's water temperature and at PROFILE_T_COOLANT_HOT: absolute / normalised grids,
                  the collapsed normalised graphs, a result table (the cooler_results figures stay as they are)
  convergence     outlet temperatures and iteration errors per scenario
  caching         caching benchmark (ambient operating point) + its figures (times, time split, speed-up, presentation sequence)
  resolution      resolution sweep (ambient operating point) + its figure
  annual          the selected solvers over the 8760 hours of the climate year + its figure
  accuracy        deviation of the solvers from the reference over the water inlet temperature + its figures (both graphs, and the pinch graph alone)

The benchmarks run their own fresh solves (benchmark_solvers.py, benchmark_annual.py); the figures are built, not shown.
"""

import logging
from dataclasses import dataclass

from src.catalogue import Item

logger = logging.getLogger(__name__)


@dataclass
class RunContext:
    """What an output may use: the solved scenarios of the run (dicts: label, ops, result, P_pump, P_fan, m_dot_ev;
    the dry one first, the precooled one if the decider engaged it), the geometry, solver settings and run modes."""
    geo: object
    settings: object
    modes: object
    scenarios: list

    @property
    def ops_ambient(self):
        return self.scenarios[0]['ops']


# --- figures of the scenarios --------------------------------------------------------------------------------

def _plot_cooler_results(_, ctx):
    from src.plots.plot_profiles import common_ticks, plot_cooler_results
    n_segments, n_elements = ctx.settings.cell_n_segments, ctx.settings.ntu_n_elements
    ticks = common_ticks([(sc['result'], sc['ops']) for sc in ctx.scenarios], ctx.geo, n_segments, n_elements)
    for sc in ctx.scenarios:                                  # identical axes in the dry and the precooled figure
        plot_cooler_results(sc['result'], sc['ops'], ctx.geo, sc['label'], n_segments=n_segments, n_elements=n_elements,
                            P_pump=sc['P_pump'], P_fan=sc['P_fan'], m_dot_ev=sc['m_dot_ev'], ticks=ticks,
                            T_amb=ctx.ops_ambient.T_a_i)


def _compute_profile_cases(ctx):
    from src.inputs.operating_conditions import ops_like
    from src.solvers import scenario_keys, solve_scenario
    ops_hot = ops_like(ctx.geo, ctx.ops_ambient, T_c_i=ctx.modes.profile_t_coolant)
    result_hot = solve_scenario(ops_hot, ctx.geo, scenario_keys(ctx.modes.solvers), ctx.settings)
    return [(ctx.ops_ambient, ctx.scenarios[0]['result']), (ops_hot, result_hot)]


def _plot_profile_grids(cases, ctx):
    from src.plots import plot_profiles
    n_segments, n_elements = ctx.settings.cell_n_segments, ctx.settings.ntu_n_elements
    plot_profiles.plot_profile_grid(cases, ctx.geo, n_segments, n_elements, normalised=False)
    plot_profiles.plot_profile_grid(cases, ctx.geo, n_segments, n_elements, normalised=True)
    plot_profiles.plot_profile_collapsed(cases, ctx.geo, n_segments, n_elements)
    plot_profiles.plot_result_table(cases, ctx.geo)


def _plot_convergence(_, ctx):
    from src.plots.plot_solver import plot_convergence
    for sc in ctx.scenarios:
        plot_convergence(sc['result'], sc['ops'], sc['label'], ctx.settings)


# --- benchmarks ----------------------------------------------------------------------------------------------

def _compute_caching(ctx):
    from src.benchmarks.benchmark_solvers import run_cache_benchmark
    return run_cache_benchmark(ctx.ops_ambient, ctx.geo, ctx.settings, ctx.modes)


def _plot_caching(data, ctx):
    from src.benchmarks.benchmark_solvers import counts_table_latex
    from src.plots.figures import SOLVER, save_text
    from src.plots.plot_solver import plot_caching
    plot_caching(data)
    logger.info("  counts table: %s", save_text(SOLVER, 'caching_counts.tex', counts_table_latex(data)))


def _compute_resolution(ctx):
    from src.benchmarks.benchmark_solvers import RESULTS_FILE, run_resolution_benchmark
    from src.plots.figures import SOLVER, save_data
    data = run_resolution_benchmark(ctx.ops_ambient, ctx.geo, ctx.settings, ctx.modes)
    save_data(SOLVER, RESULTS_FILE, data)                  # a long run is never lost to a plotting problem
    return data


def _plot_resolution(data, ctx):
    from src.plots.plot_solver import plot_resolution
    plot_resolution(data, ctx.settings)


def _compute_annual(ctx):
    from src.benchmarks.benchmark_annual import run_annual_benchmark
    return run_annual_benchmark(ctx.geo, ctx.settings, ctx.modes)


def _plot_annual(data, ctx):
    from src.plots.plot_solver import plot_annual
    plot_annual(data)


def _compute_accuracy(ctx):
    from src.benchmarks.benchmark_accuracy import run_accuracy_sweep
    return run_accuracy_sweep(ctx.ops_ambient, ctx.geo, ctx.settings, ctx.modes)


def _plot_accuracy(data, ctx):
    from src.benchmarks.benchmark_accuracy import plot_accuracy_sweep
    plot_accuracy_sweep(data)
    plot_accuracy_sweep(data, panels=('pinch',), name='accuracy_sweep_pinch')


def _log_accuracy(data, ctx):
    from src.benchmarks.benchmark_accuracy import summary_lines
    return summary_lines(data)


OUTPUTS = {item.key: item for item in (
    Item('cooler_results', "Cooler results", lambda ctx: None, _plot_cooler_results),
    Item('profile_grids', "Profile grids (two water temperatures)", _compute_profile_cases, _plot_profile_grids),
    Item('convergence', "Convergence", lambda ctx: None, _plot_convergence),
    Item('caching', "Caching benchmark", _compute_caching, _plot_caching),
    Item('resolution', "Resolution benchmark", _compute_resolution, _plot_resolution),
    Item('annual', "Annual benchmark (8760 h)", _compute_annual, _plot_annual),
    Item('accuracy', "Accuracy sweep", _compute_accuracy, _plot_accuracy, _log_accuracy),
)}
