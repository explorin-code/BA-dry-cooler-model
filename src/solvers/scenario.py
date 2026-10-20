"""
scenario.py -- running the solvers: one by key (run_solver) or a selection on one operating point (solve_scenario).
"""

import dataclasses
import time

from src.solvers.cell import solve_cell
from src.solvers.common import fill_defaults
from src.solvers.lmtd import solve_lmtd
from src.solvers.ntu import solve_ntu
from src.solvers.results import SOLVER_LABELS, ScenarioResult, SolverResult


# key -> solver function (ops, geo, settings, profile) -> SolverResult. Cell 2D / 3D differ only in settings.cell_2d,
# NTU field / table in the mode.
SOLVERS = {
    'lmtd': lambda ops, geo, settings, profile: solve_lmtd(ops=ops, geo=geo, settings=settings),
    'ntu_field': lambda ops, geo, settings, profile: solve_ntu(ops=ops, geo=geo, settings=settings, mode='field'),
    'ntu_table': lambda ops, geo, settings, profile: solve_ntu(ops=ops, geo=geo, settings=settings, mode='table', profile=profile),
    'cell_2d': lambda ops, geo, settings, profile: solve_cell(ops=ops, geo=geo, settings=dataclasses.replace(settings, cell_2d=True)),
    'cell_3d': lambda ops, geo, settings, profile: solve_cell(ops=ops, geo=geo, settings=dataclasses.replace(settings, cell_2d=False)),
}
assert SOLVERS.keys() == SOLVER_LABELS.keys()


def run_solver(key, ops, geo, settings, profile=True) -> SolverResult:
    """Runs one of the five solvers (SOLVER_LABELS) on one operating point. profile: NTU table mode also solves
    the element field once at the end, for the grids (plots)."""
    try:
        solver = SOLVERS[key]
    except KeyError:
        raise ValueError(f"Unknown solver {key!r} -- one of {list(SOLVER_LABELS)}") from None
    return solver(ops, geo, settings, profile)


def scenario_keys(selected) -> list:
    """The solvers a scenario actually runs: ONE per algorithm family of the selection (LMTD, NTU, Cell). The variants
    of a family (NTU field / table, Cell 2D / 3D) give identical results -- they differ only in speed, which is the
    business of the benchmarks. Cell runs in 3D if selected, NTU as the field (the profile solver)."""
    keys = []
    if 'lmtd' in selected:
        keys.append('lmtd')
    if any(k in selected for k in ('ntu_field', 'ntu_table')):
        keys.append('ntu_field' if 'ntu_field' in selected else 'ntu_table')
    if any(k in selected for k in ('cell_2d', 'cell_3d')):
        keys.append('cell_3d' if 'cell_3d' in selected else 'cell_2d')
    return keys


def solve_scenario(ops, geo=None, solver_keys=None, settings=None) -> ScenarioResult:
    """Runs the selected solvers (default: all five) on one operating point and times each. The solvers share
    the property cache -- use benchmark_solvers.py for fair timings."""
    ops, geo, settings = fill_defaults(ops, geo, settings)
    keys = [k for k in (solver_keys or SOLVER_LABELS)]
    results = {}
    for key in keys:
        t_start = time.perf_counter()
        results[key] = run_solver(key, ops, geo, settings)
        results[key].solve_time = time.perf_counter() - t_start
    return ScenarioResult(results)
