"""
ntu_field.py -- NTU: the element-wise e-NTU field of ONE coolant circuit and the outer-loop step that uses it.

Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2006), "Thermal Performance of Multipass Parallel and Counter-Cross-Flow
Heat Exchangers", J. Heat Transfer 129, 282-290, DOI 10.1115/1.2430719. Arrangement G^c_(N_r,1). Paper notation inside
the element functions: h = coolant (tube fluid, mixed), c = air (external fluid, unmixed). Constant U and cp within one
field solve (as in the paper); the outer loop updates them between field solves -- revisit if results look off.
The field is traversed with the shared sweep_two_way (field.py), the same as Cell; the paper's own coolant-direction-only
march is the inactive sweep_coolant_direction (same converged field, more sweeps).
"""

from math import isclose
import warnings

import numpy as np

from src.solvers.field import (calc_element_effectiveness, calc_element_outlets, cell_path_order, coolant_inlet,
                               sweep_two_way)


def solve_ntu_field(N_r: int, n_elements: int, C_h_circuit: float, C_c_element: float,
                    UA_element: float, T_c_i: float, T_a_i: float,
                    tol: float, max_iter: int, T_start=None):
    """Converged element field of ONE coolant circuit (all circuits are
    identical). Plain numbers only (no geometry or operating-condition objects).
    Returns (T_c, T_a, Q_element, T_c_o, T_a_o, residuals);
    grids are [row, element] at the physical axial position.
    Assumption (not from the paper): air leaving element (r-1, e) enters
    (r, e), no averaging -- reproduces the paper's closed forms.
    Traversal: sweep_two_way (shared with Cell) instead of the paper's
    coolant-direction-only march (Fig. 2a, item 1.7, = sweep_coolant_direction):
    same converged field, fewer sweeps.
    T_start: optional (T_c, T_a) [row, element] grids to start from (warm
    start from the previous outer iteration); None = the paper's start from
    the inlet temperatures."""
    Gamma = calc_element_effectiveness(UA_element, C_c_element)
    B = C_c_element * Gamma / C_h_circuit                      # p. 283, below Eq. (4)
    if C_c_element / C_h_circuit > 0.1:                        # Eq. (1) needs C_c << C_h
        warnings.warn(f"NTU element field: C_c^e/C_h^e = {C_c_element / C_h_circuit:.3f} is not "
                      f"<< 1 -- increase n_elements (currently {n_elements}).")

    def element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur):   # constant Gamma, B
        T_a_out, T_c_out = calc_element_outlets(T_c_in, T_a_in, Gamma, B)
        return T_c_out, T_a_out

    path = cell_path_order(N_r, n_elements)
    if T_start is None:                                       # one circuit = one tube column
        T_c = np.full((N_r, 1, n_elements), float(T_c_i))
        T_a = np.full((N_r, 1, n_elements), float(T_a_i))     # assumed initial air field
    else:
        T_c = np.array(T_start[0], dtype=float)[:, None, :]
        T_a = np.array(T_start[1], dtype=float)[:, None, :]

    residuals = []
    for _ in range(max_iter):
        change_c, change_a = sweep_two_way(T_c, T_a, element, T_c_i, T_a_i, path)
        # Own criterion (max-norm); the paper compares mean air outlets (Fig. 2a).
        residuals.append(max(change_c, change_a))
        if residuals[-1] < tol:
            break
    else:
        raise RuntimeError(f"NTU element field did not converge in {max_iter} sweeps "
                           f"(last residual {residuals[-1]:.2e} K, tolerance {tol:.0e} K).")

    Q_element = np.array([[C_h_circuit * (coolant_inlet(r, 0, e, T_c, T_c_i) - T_c[r, 0, e])
                           for e in range(n_elements)] for r in range(N_r)])
    r_last, e_last = path[-1]
    T_c_o = float(T_c[r_last, 0, e_last])                     # after the last element
    T_a_o = float(np.mean(T_a[N_r - 1, 0, :]))                # mean air outlet (Fig. 2a)
    return T_c[:, 0, :], T_a[:, 0, :], Q_element, T_c_o, T_a_o, residuals


def ntu_step(k, coolant_state, air_state, geo, ops, n_elements, field_tol, field_max_iter, field_out=None):
    """One field solve at the current k and cp; returns outlets and Q.
    Element allocation per Fig. 2(a), items 1.5-1.6. field_out (dict)
    receives the field for the final grids; a field already in it (from the
    previous outer iteration) is the warm start of this solve -- k/cp change
    little between outer iterations, so it is nearly converged already
    (sweeps per field 9, 11, 11, 11, 11 -> 9, 9, 6, 3, 1; NTU ~13 -> ~7 ms)."""
    C_h = ops.m_dot_c * coolant_state.cp                       # paper h = coolant
    C_c = ops.m_dot_a * air_state.cp                           # paper c = air
    A_element = geo.A / n_elements                             # geo.A = one tube, one row
    if not isclose(geo.N_t * geo.N_r * n_elements * A_element, geo.A_tot, rel_tol=1e-9):
        raise ValueError("NTU element areas do not add up to geo.A_tot.")

    T_c, T_a, Q_element, T_c_o, T_a_o, residuals = solve_ntu_field(
        geo.N_r, n_elements, C_h / geo.N_t, C_c / (geo.N_t * n_elements), k * A_element,
        ops.T_c_i, ops.T_a_i, field_tol, field_max_iter,
        T_start=(field_out['T_c'], field_out['T_a']) if field_out and 'T_c' in field_out else None)

    Q_dot = C_h * (ops.T_c_i - T_c_o)                          # coolant side, as Cell

    if field_out is not None:
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        field_out.update(T_c=T_c, T_a=T_a, Q_element=Q_element, residuals=residuals,
                         # Paper's whole-exchanger groups (p. 286 Eqs. 7-8, p. 289)
                         C_star=C_min / C_max, NTU=k * geo.A_tot / C_min,
                         P=(T_a_o - ops.T_a_i) / (ops.T_c_i - ops.T_a_i), R=C_c / C_h)
    return T_c_o, T_a_o, Q_dot
