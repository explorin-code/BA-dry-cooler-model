"""
run_modes.py
=============
Run-mode toggles (what a run prints, plots and measures) -- distinct from
parameters.py's physical inputs and solver_settings.py's algorithm tuning.
Defaults come from parameters.py; main.py's command-line flags override
them per run. The toggles are independent: any combination works.
"""

from dataclasses import dataclass

from src.param_loader import from_parameters


@dataclass
class RunModes:
    insight_mode: bool             # per-iteration progress + detailed results in the terminal
    plot_results: bool             # cooler-results figures (inputs, results, profiles)
    plot_convergence: bool         # convergence + iteration-error figures
    benchmark_mode: bool           # caching benchmark + solver-performance figure
    resolution_mode: bool          # resolution sweep + its figure
    cell_2d: bool                  # Cell solves one representative tube instead of all
    benchmark_resolutions: list    # Cell segments / NTU elements swept by the benchmark
    benchmark_repeats: int         # timed repeats per caching measurement (median reported)
    benchmark_sweep_repeats: int   # timed repeats per resolution point


def get_run_modes(**overrides) -> RunModes:
    """Current run modes -- parameters.py defaults, with any non-None
    keyword override (e.g. from main.py's flags) taking priority."""
    return from_parameters(RunModes, {
        'insight_mode': 'INSIGHT_MODE',
        'plot_results': 'PLOT_RESULTS',
        'plot_convergence': 'PLOT_CONVERGENCE',
        'benchmark_mode': 'BENCHMARK_MODE',
        'resolution_mode': 'RESOLUTION_MODE',
        'cell_2d': 'CELL_2D',
        'benchmark_resolutions': 'BENCHMARK_RESOLUTIONS',
        'benchmark_repeats': 'BENCHMARK_REPEATS',
        'benchmark_sweep_repeats': 'BENCHMARK_SWEEP_REPEATS',
    }, **overrides)
