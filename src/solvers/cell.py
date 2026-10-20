"""
cell.py -- the Cell solver: VDI cell method with local properties and k.

Method: VDI-Waermeatlas (2019), C1 3.1. Grid [row, tube, segment]; each cell is one small mixed(coolant)-unmixed(air)
cross-flow element with the Cabezas-Gomez element relation (as NTU) but local k and cp, traversed with the shared
sweep_two_way (field.py). Coolant / air connections: coolant_inlet, air_inlet (staggered neighbour averaging).
3D: all N_t tubes. 2D (settings.cell_2d): one representative tube -- the grid has a single tube column, so no neighbour
averaging of air; the result is broadcast to N_t tubes for output. The tube count of a run is always taken from the grid
shape, not from geo.
"""

from functools import partial
import logging

import numpy as np

from src.fluid_properties import get_air_properties, get_fluid_properties
from src.heat_transfer_core import calc_diagnostics, calc_overall_k
from src.solvers.common import fill_defaults, mean_states
from src.solvers.field import (air_inlet, calc_element_effectiveness, calc_element_outlets, cell_path_order,
                               coolant_inlet, sweep_two_way)
from src.solvers.results import SolverResult

logger = logging.getLogger(__name__)


def _cell_local_states(T_c_in, T_a_in, T_c_cur, T_a_cur, ops):
    """Local mean air temperature and fluid states of one cell; each fluid at
    its own local mean temperature (inlet + current outlet)/2 (thesis Eq.
    mean_temperatures, applied per cell)."""
    T_a_m = (T_a_in + T_a_cur) / 2.0
    props_c = get_fluid_properties(ops, (T_c_in + T_c_cur) / 2.0, ops.p_c)
    props_a = get_air_properties(ops, T_a_m, ops.p_a)
    return T_a_m, props_c, props_a


def _cell_element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur,
                  ops, geo, n_segments, m_dot_c_tube, m_dot_a_cell):
    """Outlets (T_c_out, T_a_out) of one cell: the same Cabezas-Gomez element
    relation as NTU (Eqs. 1, 3, 4), with Gamma and B from the cell's local k
    and cp instead of one value for the whole field."""
    T_a_m, props_c, props_a = _cell_local_states(T_c_in, T_a_in, T_c_cur, T_a_cur, ops)
    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_m)

    C_dot_c = m_dot_c_tube * props_c.cp
    C_dot_a = m_dot_a_cell * props_a.cp
    Gamma = calc_element_effectiveness(k_local * (geo.A / n_segments), C_dot_a)
    B = C_dot_a * Gamma / C_dot_c
    T_a_out, T_c_out = calc_element_outlets(T_c_in, T_a_in, Gamma, B)
    return T_c_out, T_a_out


def _relax_cell_grid(T_c, T_a, omega, Delta_T_hot, Delta_T_cold,
                     geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell, threshold, max_iter):
    """Iterates the cell field with sweep_two_way (shared with NTU) until both
    outlet temperatures change by less than threshold per iteration -- the
    same criterion as LMTD/NTU's outer loop."""
    history_hot, history_cold, history_T_c_o, history_T_a_o = [], [], [], []
    path = cell_path_order(geo.N_r, n_segments)
    element = partial(_cell_element, ops=ops, geo=geo, n_segments=n_segments,
                      m_dot_c_tube=m_dot_c_tube, m_dot_a_cell=m_dot_a_cell)

    for iteration in range(max_iter):
        change_c, change_a = sweep_two_way(T_c, T_a, element, ops.T_c_i, ops.T_a_i, path, omega)

        T_c_o = float(np.mean(T_c[0, :, n_segments - 1]))   # coolant leaves row 0
        T_a_o = float(np.mean(T_a[geo.N_r - 1, :, :]))       # air leaves the last row

        new_Delta_T_hot = ops.T_c_i - T_a_o
        new_Delta_T_cold = T_c_o - ops.T_a_i
        diff_hot = new_Delta_T_hot - Delta_T_hot
        diff_cold = new_Delta_T_cold - Delta_T_cold
        Delta_T_hot, Delta_T_cold = new_Delta_T_hot, new_Delta_T_cold
        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_T_c_o.append(T_c_o)
        history_T_a_o.append(T_a_o)

        logger.debug("[Cell] Air change: %.5f | Coolant change: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     change_a, change_c, T_c_o, T_a_o, omega)
        if abs(diff_hot) < threshold and abs(diff_cold) < threshold:
            break
    else:
        raise RuntimeError(f"[Cell] did not converge in {max_iter} iterations "
                           f"(last changes {diff_hot:.2e} / {diff_cold:.2e} K, threshold {threshold:.0e} K) -- "
                           f"cap too low, or threshold below the property-cache floor (~1.5e-5 K)?")

    return (T_c, T_a, T_c_o, T_a_o,
            history_hot, history_cold, history_T_c_o, history_T_a_o)


def _cell_local_grids(T_c, T_a, ops, geo, n_segments):
    """k, dQ/dL and Re/Nu per cell on the converged grids (same local states as during the iteration). dQ/dL is the
    coolant's energy balance over the cell (as in NTU's Q_element): the cell's outlet temperatures are not its mean
    temperatures, so k dA (T_c - T_a) would not add up to Q."""
    shape = T_c.shape
    k_grid, dQdL_grid = np.zeros(shape), np.zeros(shape)
    local = {key: np.zeros(shape) for key in ('Re_a', 'Nu_a', 'Re_c', 'Nu_c')}
    dL_cell = geo.l / n_segments
    m_dot_c_tube = ops.m_dot_c / geo.N_t

    for r in range(geo.N_r):
        for t in range(shape[1]):
            for s in range(n_segments):
                T_c_in = coolant_inlet(r, t, s, T_c, ops.T_c_i)
                T_a_m, props_c, props_a = _cell_local_states(T_c_in, air_inlet(r, t, s, T_a, ops.T_a_i),
                                                             T_c[r, t, s], T_a[r, t, s], ops)
                k_grid[r, t, s] = calc_overall_k(geo, ops, props_c, props_a, T_a_m)
                dQdL_grid[r, t, s] = m_dot_c_tube * props_c.cp * (T_c_in - T_c[r, t, s]) / dL_cell
                diag = calc_diagnostics(geo, ops, props_c, props_a, T_a_m)
                for key in local:
                    local[key][r, t, s] = diag[key]
    return k_grid, dQdL_grid, local


def solve_cell(n_segments: int = None, omega: float = None, ops=None, geo=None, settings=None):
    """omega: per-cell relaxation (default settings.cell_omega = 1.0; with
    omega < 1 the change-per-iteration criterion fires before convergence --
    see docs/model-notes.md). n_segments per tube and row (default settings.cell_n_segments)."""
    ops, geo, settings = fill_defaults(ops, geo, settings)
    if omega is None:
        omega = settings.cell_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments

    m_dot_c_tube = ops.m_dot_c / geo.N_t                       # per tube, also in 2D
    m_dot_a_cell = ops.m_dot_a / (geo.N_t * n_segments)
    n_tubes = 1 if settings.cell_2d else geo.N_t
    T_c = np.full((geo.N_r, n_tubes, n_segments), float(ops.T_c_i))
    T_a = np.full((geo.N_r, n_tubes, n_segments), float(ops.T_a_i))     # float: whole-number inputs must not make an integer grid

    (T_c, T_a, T_c_o, T_a_o,
     history_hot, history_cold, history_T_c_o, history_T_a_o) = _relax_cell_grid(
        T_c, T_a, omega, settings.Delta_T_hot_init, settings.Delta_T_cold_init,
        geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell,
        threshold=settings.convergence_threshold, max_iter=settings.cell_max_iter)

    # Whole-exchanger Q and diagnostics at the overall mean temperatures; the
    # reported k is the mean of the local k Cell actually used (equal cell areas).
    coolant_state, air_state = mean_states(ops, T_c_o, T_a_o)
    Q_dot = ops.m_dot_c * coolant_state.cp * (ops.T_c_i - T_c_o)
    diagnostics = calc_diagnostics(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                                   T_a_m=(ops.T_a_i + T_a_o) / 2.0)
    k_grid, dQdL_grid, local_diagnostics = _cell_local_grids(T_c, T_a, ops, geo, n_segments)

    if n_tubes == 1:                                           # 2D: broadcast to all tubes for output
        widen = lambda grid: np.broadcast_to(grid, (geo.N_r, geo.N_t, n_segments)).copy()
        T_c, T_a, k_grid, dQdL_grid = map(widen, (T_c, T_a, k_grid, dQdL_grid))
        local_diagnostics = {key: widen(grid) for key, grid in local_diagnostics.items()}

    return SolverResult(float(np.mean(k_grid)), Q_dot, T_c_o, T_a_o,
                        history_hot, history_cold, history_T_c_o, history_T_a_o, diagnostics,
                        T_c_grid=T_c, T_a_grid=T_a, k_grid=k_grid, dQdL_grid=dQdL_grid,
                        local_diagnostics=local_diagnostics)
