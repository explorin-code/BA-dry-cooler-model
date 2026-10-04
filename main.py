"""
main.py
=======
Entry point. Runs the dry-cooler scenario(s) via scenario_pipeline.run_scenarios()
and shows the resulting figures together. See scenario_pipeline.py for what
actually gets solved/plotted.

Run modes default to parameters.py (INSIGHT_MODE, PLOT_RESULTS,
PLOT_CONVERGENCE, BENCHMARK_MODE, RESOLUTION_MODE, CELL_2D); each can be overridden
per run, in any combination:
    python main.py --no-insight --no-plots --no-benchmark   # fast: summary only
    python main.py --convergence                            # + convergence figures
    python main.py --benchmark --resolution                 # + both performance figures
"""

import argparse
import logging
import time

import matplotlib.pyplot as plt

from src.run_modes import get_run_modes
from src.scenario_pipeline import run_scenarios


def parse_args():
    parser = argparse.ArgumentParser(description="Dry-cooler model: LMTD / NTU / Cell comparison.")
    parser.add_argument("--insight", action=argparse.BooleanOptionalAction, default=None,
                        help="per-iteration progress + detailed results in the terminal")
    parser.add_argument("--plots", action=argparse.BooleanOptionalAction, default=None,
                        help="cooler-results figures (inputs, results, spatial profiles)")
    parser.add_argument("--convergence", action=argparse.BooleanOptionalAction, default=None,
                        help="convergence + iteration-error figures")
    parser.add_argument("--benchmark", action=argparse.BooleanOptionalAction, default=None,
                        help="caching benchmark + solver-performance figure")
    parser.add_argument("--resolution", action=argparse.BooleanOptionalAction, default=None,
                        help="resolution sweep + its figure")
    parser.add_argument("--cell-2d", action=argparse.BooleanOptionalAction, default=None,
                        help="Cell solves one representative tube (2D) instead of all tubes (3D)")
    return parser.parse_args()


def setup_logging(insight_mode: bool):
    """Project messages only (logger "src"); third-party loggers stay at
    WARNING. Insight mode shows DEBUG (per-iteration progress, details),
    otherwise INFO (summary per scenario)."""
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    logging.getLogger("src").setLevel(logging.DEBUG if insight_mode else logging.INFO)
    logging.captureWarnings(True)    # route warnings.warn (e.g. correlation ranges) through logging


def run():
    args = parse_args()
    modes = get_run_modes(insight_mode=args.insight, plot_results=args.plots,
                          plot_convergence=args.convergence, benchmark_mode=args.benchmark,
                          resolution_mode=args.resolution, cell_2d=args.cell_2d)
    setup_logging(modes.insight_mode)

    t_start = time.perf_counter()
    run_scenarios(modes)
    logging.getLogger("src").info("Total run time: %.2f s", time.perf_counter() - t_start)

    if plt.get_fignums():
        plt.show()


if __name__ == "__main__":
    run()
