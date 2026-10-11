"""
solvers -- the three solvers for the same dry-cooler problem, as a package:

  results.py     SolverResult, ScenarioResult, SOLVER_LABELS / SOLVER_FAMILY / DETAIL_ORDER (the five selectable solvers)
  common.py      defaults of a call, states at the mean temperatures
  field.py       the element field of NTU and Cell: coolant path, element relation, connections, the two traversals
  relaxation.py  the outer loop shared by LMTD and NTU (with the oscillation damping)
  lmtd.py        LMTD -- the whole exchanger as pure counterflow
  ntu_field.py   NTU -- the element-wise e-NTU field of one circuit (Cabezas-Gomez et al. 2007) and its outer-loop step
  ntu_table.py   NTU 'table' mode -- P(NTU, R) built once, stored in cache/, interpolated
  ntu.py         solve_it_NTU: field or table mode
  cell.py        Cell -- VDI cell method with local properties and k
  scenario.py    run_solver (one solver by key), scenario_keys / solve_scenario (a selection on one operating point)

Pure computation: progress is logged at DEBUG level (insight mode), plotting lives in plot_profiles.py / plot_solver.py.
The names below are the package's public interface; the NTU table internals are `src.solvers.ntu_table`.
"""

from src.solvers.field import (calc_element_effectiveness, calc_element_outlets, cell_path_order, get_coolant_inlet,
                               get_staggered_air_inlet, sweep_coolant_direction, sweep_two_way)
from src.solvers.cell import solve_it_cell
from src.solvers.lmtd import calc_LMTD, solve_it_LMTD
from src.solvers.ntu import solve_it_NTU
from src.solvers.ntu_field import solve_ntu_field
from src.solvers.results import DETAIL_ORDER, SOLVER_FAMILY, SOLVER_LABELS, ScenarioResult, SolverResult
from src.solvers.scenario import run_solver, scenario_keys, solve_scenario

__all__ = ['calc_element_effectiveness', 'calc_element_outlets', 'cell_path_order', 'get_coolant_inlet',
           'get_staggered_air_inlet', 'sweep_coolant_direction', 'sweep_two_way', 'solve_it_cell', 'calc_LMTD',
           'solve_it_LMTD', 'solve_it_NTU', 'solve_ntu_field', 'DETAIL_ORDER', 'SOLVER_FAMILY', 'SOLVER_LABELS',
           'ScenarioResult', 'SolverResult', 'run_solver', 'scenario_keys', 'solve_scenario']
