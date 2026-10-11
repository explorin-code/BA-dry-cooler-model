"""
outputs.py
===========
The catalogue of what the solver comparison can produce besides its terminal summary (parameters.SOLVER_OUTPUTS,
`main.py --outputs`, the GUI's multiselect). Each output is a catalogue.Item whose steps get the run context:

  cooler_results  the cooler-results figure per scenario (dry / precooled), identical axes
  convergence     outlet temperatures and iteration errors per scenario
  caching         caching benchmark (ambient operating point) + its figures (linear and logarithmic time axes)
  resolution      resolution sweep (ambient operating point) + its figure
  annual          the selected solvers over the 8760 hours of the climate year + its figure

The benchmarks run their own fresh solves (benchmark_solvers.py, benchmark_annual.py); the figures are built, not shown.
"""

from dataclasses import dataclass

from src.catalogue import Item


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
    from src.plot_profiles import common_ticks, plot_cooler_results
    n_segments, n_elements = ctx.settings.cell_n_segments, ctx.settings.ntu_n_elements
    ticks = common_ticks([(sc['result'], sc['ops']) for sc in ctx.scenarios], ctx.geo, n_segments, n_elements)
    for sc in ctx.scenarios:                                  # identical axes in the dry and the precooled figure
        plot_cooler_results(sc['result'], sc['ops'], ctx.geo, sc['label'], n_segments=n_segments, n_elements=n_elements,
                            P_pump=sc['P_pump'], P_fan=sc['P_fan'], m_dot_ev=sc['m_dot_ev'], ticks=ticks)


def _plot_convergence(_, ctx):
    from src.plot_solver import plot_convergence
    for sc in ctx.scenarios:
        plot_convergence(sc['result'], sc['ops'], sc['label'], ctx.settings)


# --- benchmarks ----------------------------------------------------------------------------------------------

def _compute_caching(ctx):
    from src.benchmark_solvers import run_cache_benchmark
    return run_cache_benchmark(ctx.ops_ambient, ctx.geo, ctx.settings, ctx.modes)


def _plot_caching(data, ctx):
    from src.plot_solver import plot_solver_performance
    for log in (False, True):                                 # the same page on linear and on logarithmic time axes
        plot_solver_performance(ctx.settings, data, log=log)


def _compute_resolution(ctx):
    from src.benchmark_solvers import run_resolution_benchmark
    return run_resolution_benchmark(ctx.ops_ambient, ctx.geo, ctx.settings, ctx.modes)


def _plot_resolution(data, ctx):
    from src.benchmark_solvers import plot_resolution
    plot_resolution(data)


def _compute_annual(ctx):
    from src.benchmark_annual import run_annual_benchmark
    return run_annual_benchmark(ctx.geo, ctx.settings, ctx.modes)


def _plot_annual(data, ctx):
    from src.benchmark_annual import plot_annual_benchmark
    plot_annual_benchmark(data)


OUTPUTS = {item.key: item for item in (
    Item('cooler_results', "Cooler results", lambda ctx: None, _plot_cooler_results),
    Item('convergence', "Convergence", lambda ctx: None, _plot_convergence),
    Item('caching', "Caching benchmark", _compute_caching, _plot_caching),
    Item('resolution', "Resolution benchmark", _compute_resolution, _plot_resolution),
    Item('annual', "Annual benchmark (8760 h)", _compute_annual, _plot_annual),
)}
