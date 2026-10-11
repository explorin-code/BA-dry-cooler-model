"""
ntu.py -- the NTU solver: 'field' mode solves the element field in every outer iteration (the internal profile comes
for free), 'table' mode looks up the tabulated P(NTU_a, R) (ntu_table.py) -- for many operating points at a constant row
count; no internal field unless profile=True (one field solve at the end).
"""

from functools import partial

import numpy as np

from src.solvers.common import fill_defaults, mean_states
from src.solvers.ntu_field import ntu_step
from src.solvers.ntu_table import get_ntu_table, ntu_table_step
from src.solvers.relaxation import relax_lmtd_ntu


def solve_it_NTU(omega: float = None, ops=None, geo=None, settings=None, n_elements: int = None,
                 mode: str = 'field', profile: bool = True):
    """omega: relaxation of the outer property/k loop (default
    settings.ntu_omega = 1.0, undamped); n_elements per tube and row
    (default settings.ntu_n_elements); mode 'field' or 'table' (run_solver passes it per
    solver variant). profile: in 'table' mode, solve the field
    once at the converged state for the grids (plots) -- pass False for
    multi-run use, the grids are then None."""
    ops, geo, settings = fill_defaults(ops, geo, settings)
    if omega is None:
        omega = settings.ntu_omega
    if n_elements is None:
        n_elements = settings.ntu_n_elements

    last_field = {}
    threshold = settings.convergence_threshold
    field_tol = threshold / 10                       # the inner field is always more accurate than the criterion it feeds
    field_step = partial(ntu_step, n_elements=n_elements, field_tol=field_tol,
                         field_max_iter=settings.ntu_field_max_iter, field_out=last_field)
    if mode == 'field':
        step_fn = field_step
    elif mode == 'table':
        step_fn = partial(ntu_table_step, table=get_ntu_table(geo.N_r, n_elements), field_step=field_step)
    else:
        raise ValueError(f"Unknown NTU mode {mode!r} -- 'field' or 'table'.")
    result = relax_lmtd_ntu(step_fn, "[NTU]", omega, threshold, ops, geo, settings)

    if mode == 'table' and profile:                            # one field solve for the grids
        field_step(result.k, *mean_states(ops, result.T_c_o, result.T_a_o), geo, ops)
    if 'T_c' in last_field:
        # Grids (row, tube, element) like Cell's: the one solved circuit, broadcast to all tubes.
        shape = (geo.N_r, geo.N_t, n_elements)
        widen = lambda grid: np.broadcast_to(grid[:, None, :], shape).copy()
        result.T_c_grid = widen(last_field['T_c'])
        result.T_a_grid = widen(last_field['T_a'])
        result.k_grid = np.full(shape, result.k)
        result.dQdL_grid = widen(last_field['Q_element'] / (geo.l / n_elements))
    return result
