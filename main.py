"""
main.py
=======
Entry point with two workflows (the GUI, gui.py, has the same two modes):

  python main.py solver [flags]      SOLVER COMPARISON: ONE cooler (parameters.py: GEOMETRY_SET P1/P2/P3 + overrides)
                                     in an ambient state, the selected solvers (SOLVERS: LMTD, NTU field / table,
                                     Cell 2D / 3D) compared, ambient and precooled scenario, plus the solver
                                     benchmarks. Default when no workflow is named.
  python main.py compare [analyses]  the standard coolers P1-P3 raced against each other over the Munich climate
                                     year with ONE solver (COMPARE_SOLVER): geometry, fan control, sensitivity ...

Solver-comparison run modes default to parameters.py (SOLVERS, INSIGHT_MODE, SOLVER_OUTPUTS: the catalogue of
src/outputs.py); each can be overridden per run:
    python main.py --no-insight --outputs                   # fast: summary only
    python main.py solver --cooler P3                       # the P3 cell instead of the preset
    python main.py solver --solvers lmtd ntu_table          # only these solvers (lmtd ntu_field ntu_table cell_2d cell_3d)
    python main.py solver --outputs cooler_results convergence   # these figures only
    python main.py solver --outputs caching resolution      # both solver-performance benchmarks
    python main.py solver --outputs annual                  # annual solver benchmark (8760 h, minutes; Cell on every 4th hour, extrapolated)
    python main.py solver --outputs annual --annual-cell-stride 1 --annual-time-budget 0   # the whole year for every solver (Cell 3D: ~2 h)

Cooler comparison (analyses: see comparison/analyses.py; none named = COMPARE_ANALYSES in parameters.py):
    python main.py compare geometry                         # plot all geometries
    python main.py compare fan_control --stride 24          # year with the variable-speed fan, every 24th hour
    python main.py compare all --solver lmtd                # everything with the fastest solver
"""

import argparse
import logging
import time

import matplotlib.pyplot as plt

import src.parameters as parameters
from src.cooler_sets import CELL_SETS
from src.outputs import OUTPUTS
from src.solvers import SOLVER_LABELS

WORKFLOWS = ('solver', 'compare')


def add_solver_arguments(parser):
    flag = lambda name, help_: parser.add_argument(f"--{name}", action=argparse.BooleanOptionalAction, default=None, help=help_)
    flag("insight", "per-iteration progress + detailed results in the terminal")
    parser.add_argument("--outputs", nargs="*", choices=list(OUTPUTS),
                        help="what to produce besides the summary (default SOLVER_OUTPUTS; none = only the summary): "
                             + ", ".join(f"{key} ({item.title})" for key, item in OUTPUTS.items()))
    parser.add_argument("--cooler", choices=list(CELL_SETS), help="cooler cell to use (default GEOMETRY_SET)")
    parser.add_argument("--solvers", nargs="+", choices=list(SOLVER_LABELS), help="solvers to run (default SOLVERS)")
    parser.add_argument("--annual-cell-stride", type=int, help="annual benchmark: Cell on every n-th hour (default ANNUAL_CELL_STRIDE; 1 = every hour)")
    parser.add_argument("--annual-time-budget", type=float, help="annual benchmark: seconds per solver before it is extrapolated (default ANNUAL_TIME_BUDGET_S; 0 = no limit)")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Dry-cooler model: solver comparison (one cooler) or cooler comparison (P1-P3 over the climate year).")
    sub = parser.add_subparsers(dest='workflow')
    add_solver_arguments(sub.add_parser('solver', help="solver comparison: one cooler, the selected solvers (default)"))
    from comparison.__main__ import add_arguments
    add_arguments(sub.add_parser('compare', help="the standard coolers raced against each other over the climate year"))
    argv = list(argv if argv is not None else __import__('sys').argv[1:])
    if not argv or argv[0] not in (*WORKFLOWS, '-h', '--help'):
        argv.insert(0, 'solver')                            # `python main.py --no-plots` = solver comparison
    return parser.parse_args(argv)


def setup_logging(insight_mode: bool = False):
    """Project messages only (logger "src"); third-party loggers stay at WARNING. Insight mode shows DEBUG
    (per-iteration progress, details), otherwise INFO (summary per scenario)."""
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    logging.getLogger("src").setLevel(logging.DEBUG if insight_mode else logging.INFO)
    logging.captureWarnings(True)    # route warnings.warn (e.g. correlation ranges) through logging


def run_solver(args):
    from src.run_modes import get_run_modes
    from src.scenario_pipeline import run_scenarios
    if args.cooler:
        parameters.GEOMETRY_SET = args.cooler
    modes = get_run_modes(insight_mode=args.insight, outputs=args.outputs, solvers=args.solvers,
                          annual_cell_stride=args.annual_cell_stride, annual_time_budget=args.annual_time_budget)
    setup_logging(modes.insight_mode)
    run_scenarios(modes)


def run_compare(args):
    from comparison.__main__ import apply_arguments
    from comparison.runner import run_analyses
    setup_logging()
    run_analyses(apply_arguments(args), plots_only=args.plots_only, keep_open=True)


def run():
    args = parse_args()
    t_start = time.perf_counter()
    (run_solver if args.workflow == 'solver' else run_compare)(args)
    logging.getLogger("src").info("Total run time: %.2f s", time.perf_counter() - t_start)
    if plt.get_fignums():
        plt.show()


if __name__ == "__main__":
    run()
