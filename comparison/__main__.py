"""python -m comparison [analysis ...] [--plots-only] [--stride N] [--solver S] [--sensitivity-set P2] [--sets P1 P2 ...]

Without names the analyses of COMPARE_ANALYSES in parameters.py run; `all` runs every analysis."""

import argparse
import logging

import src.parameters as parameters
from comparison.analyses import ANALYSES
from src.cooler_sets import CELL_SETS
from src.solvers import SOLVER_LABELS
from comparison.runner import run_analyses


def add_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('analyses', nargs='*', metavar='ANALYSIS', choices=[*ANALYSES, 'all', []],
                        help="analyses to run: " + ', '.join(ANALYSES) + ", all (default: COMPARE_ANALYSES in parameters.py)")
    parser.add_argument('--plots-only', action='store_true', help="re-plot from the stored CSVs without computing")
    parser.add_argument('--stride', type=int, help="year analyses: use every n-th hour (quick test)")
    parser.add_argument('--solver', choices=list(SOLVER_LABELS), help="solver of the analyses (default COMPARE_SOLVER)")
    parser.add_argument('--sensitivity-set', choices=list(CELL_SETS), help="cooler the sensitivity study varies (default COMPARE_SENSITIVITY_SET)")
    parser.add_argument('--sets', nargs='+', help="coolers to compare, e.g. P1 P3 (default COMPARE_SETS)")


def apply_arguments(args):
    """Command-line values replace the parameters.py values for this run; returns the analyses to run (or None)."""
    if args.stride:
        parameters.COMPARE_YEAR_STRIDE = args.stride
    if args.solver:
        parameters.COMPARE_SOLVER = args.solver
    if args.sets:
        parameters.COMPARE_SETS = args.sets
    if args.sensitivity_set:
        parameters.COMPARE_SENSITIVITY_SET = args.sensitivity_set
    if 'all' in args.analyses:
        return list(ANALYSES)
    return args.analyses or None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='python -m comparison', description=__doc__.splitlines()[0])
    add_arguments(parser)
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    logging.getLogger('src').setLevel(logging.INFO)
    run_analyses(apply_arguments(args), plots_only=args.plots_only)
