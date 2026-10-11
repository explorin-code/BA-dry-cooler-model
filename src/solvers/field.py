"""
field.py -- the element field shared by NTU and Cell: the coolant path, the element relation (Cabezas-Gomez, Navarro &
Saiz-Jabardo 2007, Eqs. 1, 3, 4), the connections between elements and the two traversals.
"""

import numpy as np
from math import exp


def cell_path_order(N_r, n_segments):
    """(row, segment) in coolant flow order: enters the last row, serpentine
    towards row 0, axial direction alternating per row. Used by NTU, Cell
    and analysis.py, so all three walk the same path."""
    path = []
    for r in range(N_r - 1, -1, -1):
        s_range = range(0, n_segments) if r % 2 == 0 else range(n_segments - 1, -1, -1)
        for s in s_range:
            path.append((r, s))
    return path


def calc_element_effectiveness(UA_e: float, C_c_e: float) -> float:
    """Element effectiveness Gamma = 1 - exp(-UA_e / C_c_e) of a small
    mixed(tube)-unmixed(air) cross-flow element, valid for C_c_e << C_h_e.
    Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), J. Heat Transfer
    129, p. 283, Eq. (1)."""
    return 1.0 - exp(-UA_e / C_c_e)


def calc_element_outlets(T_hot_in: float, T_cold_in: float, Gamma: float, B: float):
    """Element outlet temperatures (T_cold_out, T_hot_out), B = C_c Gamma / C_h.
    Source: Cabezas-Gomez et al. (2007), p. 283, Eqs. (3)-(4)."""
    T_cold_out = ((B + 2.0 * (1.0 - Gamma)) / (2.0 + B)) * T_cold_in + (2.0 * Gamma / (2.0 + B)) * T_hot_in
    T_hot_out = ((2.0 - B) / (2.0 + B)) * T_hot_in + (2.0 * B / (2.0 + B)) * T_cold_in
    return T_cold_out, T_hot_out


# --- Field connections and traversal (NTU + Cell) ----------------------------
# Grids are [row, tube, segment]; NTU solves one circuit, i.e. a single tube
# column. Both solvers pass an element function
#     element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur) -> (T_c_out, T_a_out)
# where T_c_cur/T_a_cur are the element's current outlet values (Cell needs
# them for its local mean temperatures, NTU ignores them). Which traversal a
# solver uses is just the function name called in its iteration loop:
# sweep_two_way (active) or sweep_coolant_direction (inactive, for comparison).

def coolant_inlet(r, t, s, T_c, T_c_i):
    """Coolant entering element (r, t, s): serpentine via return bends,
    entering at the last row, direction alternating per row (cell_path_order)."""
    N_r, _, n_segments = T_c.shape
    forward = r % 2 == 0                                       # even rows run s = 0 -> end
    s_first = 0 if forward else n_segments - 1
    if s == s_first:                                           # row entry: from the U-bend
        return T_c_i if r == N_r - 1 else T_c[r + 1, t, s]
    return T_c[r, t, s - 1 if forward else s + 1]


def _coolant_predecessors(path):
    """{(r, s): (r, s) of the element upstream on the coolant path, None at
    the circuit inlet} -- the same serpentine as coolant_inlet, looked up
    once per sweep instead of per element."""
    return dict(zip(path, [None] + path[:-1]))


def air_inlet(r, t, s, T_a, T_a_i):
    """Air entering element (r, t, s): mean of the two tubes of the row
    upstream it sits between (staggered bank); edge tube, a single tube
    column (NTU, Cell 2D) or row 0: straight through."""
    if r == 0:
        return T_a_i
    t_neighbour = t - 1 if r % 2 == 0 else t + 1
    if 0 <= t_neighbour < T_a.shape[1]:
        return (T_a[r - 1, t, s] + T_a[r - 1, t_neighbour, s]) / 2.0
    return T_a[r - 1, t, s]


def get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo):
    """Coolant entering cell (r, t, s) -- see coolant_inlet.
    Source: Kröger, Air-Cooled Heat Exchangers and Cooling Towers, Vol. 1, Sec. 3.5."""
    return coolant_inlet(r, t, s, T_c, ops.T_c_i)


def get_staggered_air_inlet(r, t, s, T_a, ops, geo):
    """Air entering cell (r, t, s) -- see air_inlet. Modelling assumption
    from the staggered bank (Kröger, Vol. 1, Sec. 5.1)."""
    return air_inlet(r, t, s, T_a, ops.T_a_i)


def _field_changes(T_c, T_a, T_c_old, T_a_old, omega):
    """Largest raw change per fluid over one iteration [K]. The stored value
    is old + omega (raw - old), so the raw change is the stored one / omega --
    computed once per sweep instead of per element. A diverged field (NaN/inf)
    raises: it would otherwise pass every 'change < tol' test. Divergence
    happens e.g. with omega != 1 per element: each pass is an exact
    sequential solve along its fluid, so per-element relaxation compounds
    along the whole path."""
    change_c = float(np.max(np.abs(T_c - T_c_old))) / omega
    change_a = float(np.max(np.abs(T_a - T_a_old))) / omega
    if not (np.isfinite(change_c) and np.isfinite(change_a)):
        raise RuntimeError("Field traversal diverged (non-finite temperatures) -- check omega.")
    return change_c, change_a


def sweep_two_way(T_c, T_a, element, T_c_i, T_a_i, path, omega=1.0):
    """ACTIVE traversal of NTU and Cell: one iteration over the field, each
    fluid swept in its own flow direction --
      pass 1: coolant along its path, air inlets from the previous iteration;
      pass 2: air row 0 -> last row, with the live air of the row upstream
              and the coolant from pass 1.
    Converges to the same field as sweep_coolant_direction (same element
    relation and connections) in fewer iterations: ~2x at 6 rows, ~3.5x at
    24 (see docs/review-needed.md), but evaluates every element twice per iteration.
    Updates T_c/T_a in place (relaxed with omega) and returns the largest
    raw change (coolant, air) [K]."""
    T_c_old, T_a_old = T_c.copy(), T_a.copy()
    N_r, n_tubes, n_segments = T_c.shape
    pred = _coolant_predecessors(path)
    relax = omega != 1.0

    for r, s in path:                                          # pass 1: coolant
        p = pred[r, s]
        for t in range(n_tubes):
            T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
            T_a_in = air_inlet(r, t, s, T_a_old, T_a_i)
            raw, _ = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a_old[r, t, s])
            T_c[r, t, s] = T_c_old[r, t, s] + omega * (raw - T_c_old[r, t, s]) if relax else raw

    for r in range(N_r):                                       # pass 2: air
        for t in range(n_tubes):
            for s in range(n_segments):
                p = pred[r, s]
                T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
                T_a_in = air_inlet(r, t, s, T_a, T_a_i)
                _, raw = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a[r, t, s])
                T_a[r, t, s] = T_a_old[r, t, s] + omega * (raw - T_a_old[r, t, s]) if relax else raw
    return _field_changes(T_c, T_a, T_c_old, T_a_old, omega)


def sweep_coolant_direction(T_c, T_a, element, T_c_i, T_a_i, path, omega=1.0):
    """INACTIVE -- not called anywhere; kept for comparison with
    sweep_two_way (swap the function name in solve_ntu_field and/or
    _relax_cell_grid). One pass along the coolant path only, as in
    Cabezas-Gomez et al. (2007), Fig. 2(a), item 1.7: each element updates
    both outlets, with air inlets from the previous iteration (the row
    upstream is only reached later on the coolant path). Same converged
    field, but more iterations, and at the same stopping threshold a larger
    remaining error (convergence is slower).
    Same in-place update and return value as sweep_two_way."""
    T_c_old, T_a_old = T_c.copy(), T_a.copy()
    n_tubes = T_c.shape[1]
    pred = _coolant_predecessors(path)
    relax = omega != 1.0

    for r, s in path:
        p = pred[r, s]
        for t in range(n_tubes):
            T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
            T_a_in = air_inlet(r, t, s, T_a_old, T_a_i)
            raw_c, raw_a = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a_old[r, t, s])
            T_c[r, t, s] = T_c_old[r, t, s] + omega * (raw_c - T_c_old[r, t, s]) if relax else raw_c
            T_a[r, t, s] = T_a_old[r, t, s] + omega * (raw_a - T_a_old[r, t, s]) if relax else raw_a
    return _field_changes(T_c, T_a, T_c_old, T_a_old, omega)
