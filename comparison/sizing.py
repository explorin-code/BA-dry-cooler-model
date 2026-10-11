"""
sizing.py -- the sized module of a cooler cell: the geometry a comparison analyses.

Tube length and rows are the same fixed values for all coolers (COMPARE_TUBE_LENGTH_MM, COMPARE_N_ROWS) -- the length
changes the cooler's behaviour, not just its size. The tube count N_t is the cell's FIXED reference count
(`CellSet.N_t_ref`, or COMPARE_N_TUBES); `design_tubes` is how that count was derived: the SMALLEST number with
which the cooler reaches the target outlet temperature at the design ambient (dry, design face velocity).
`tubes_for_flow` is the other rule, a tube count from a tube Reynolds number (the tubes / wall studies). Duty and outlet
temperature at other conditions are RESULTS.
"""

import dataclasses
import logging
import math

from comparison import config, evaluate as ev
from src.cooler_sets import CELL_SETS
from src.fluid_properties import get_fluid_properties
from src.quiet import quiet_solvers

logger = logging.getLogger(__name__)

N_T_MAX = 200                  # upper limit of the design_tubes search


def tubes_for_flow(cfg, cell, pt, m_dot_c=None, Re=None) -> int:
    """N_t that carries m_dot_c with every tube at the tube Reynolds number Re (ideal equal split);
    pt: ambient state dict (only sets the property evaluation, the water enters at T_c_i)."""
    m_dot_c = cfg.m_dot_c if m_dot_c is None else m_dot_c
    Re = cfg.Re_c_design if Re is None else Re
    geo1 = ev.make_geometry(cell, 1, 1.0)
    ops1 = ev.make_ops(cfg, geo1, pt['T_a_i'], pt['phi'], pt['p_a'], m_dot_c=m_dot_c)
    state = get_fluid_properties(ops1, ops1.T_c_i, ops1.p_c)
    m_tube = Re * math.pi / 4 * (cell.d_i_mm / 1e3) * state.mu                  # m_dot = Re mu pi d_i / 4
    return max(1, round(m_dot_c / m_tube))


def design_tubes(cfg, key, m_dot_c=None) -> int:
    """Design-based tube count (how cooler_sets' N_t_ref was derived; `python -m comparison.sizing` prints it for the
    current settings next to N_t_ref -- refresh N_t_ref from it when the cells or the solver change): the smallest
    N_t (1 ... N_T_MAX) with which cooler `key` (cfg.n_rows rows, tube length cfg.l_tube) reaches cfg.T_c_target at the
    design ambient, dry, at the design face velocity. If the target is never reached, the N_t with the lowest outlet
    temperature is returned (with a warning)."""
    m_dot_c = cfg.m_dot_c if m_dot_c is None else m_dot_c
    pt, method = cfg.design_point(), cfg.solver(table_ok=False)
    cell = cfg.cell(key)
    best = (None, float('inf'))
    for n in range(1, N_T_MAX + 1):
        geo = ev.make_geometry(cell, n, cfg.l_tube)
        ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], m_dot_c=m_dot_c)
        with quiet_solvers():
            T_c_o = ev.run_point(cfg, geo, ops, method, 'dry')['T_c_o']
        if T_c_o <= cfg.T_c_target:
            logger.info("[Sizing] %s: %d tubes per row reach %.1f degC at the design point (T_c,o = %.2f)", key, n, cfg.T_c_target, T_c_o)
            return n
        if T_c_o < best[1]:
            best = (n, T_c_o)
    logger.warning("[Sizing] %s: target %.1f degC not reached with up to %d tubes -- using %d (T_c,o = %.2f)", key, cfg.T_c_target, N_T_MAX, *best)
    return best[0]


def get_sized(cfg, key, m_dot_c=None, N_r=None):
    """Geometry of cooler `key`: N_t from cfg.n_tubes, default the cell's fixed reference count N_t_ref (so a row
    sweep keeps the tube count), N_r: rows (default cfg.n_rows), tube length cfg.l_tube."""
    cell = cfg.cell(key)
    N_t = int(cfg.n_tubes) if cfg.n_tubes else CELL_SETS[key].N_t_ref
    if N_r is not None:
        cell = dataclasses.replace(cell, N_r=N_r)
    return ev.make_geometry(cell, N_t, cfg.l_tube)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    cfg = config.settings()
    for key in cfg.sets:
        print(f"{key}: design-based N_t = {design_tubes(cfg, key)} (fixed reference N_t_ref = {CELL_SETS[key].N_t_ref})")
