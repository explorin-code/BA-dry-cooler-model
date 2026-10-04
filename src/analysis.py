"""
analysis.py
============
Reconstructs a comparable spatial profile (temperature, k, dQ/dL along the
coolant's flow path) for each solver from its already-converged output.
Cell's and NTU's profiles are genuine reductions of the element/cell grids
they already computed. LMTD's is a post-hoc reconstruction from its own
converged scalars -- NOT something it tracks internally. See CLAUDE.md.
"""

from dataclasses import dataclass
import numpy as np

from src.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import cell_path_order


@dataclass
class ProfileData:
    x_frac: object      # 0 (coolant inlet) to 1 (coolant outlet)
    theta_c: object
    theta_a: object
    k: object
    dQdL: object


# AI-REVIEW: closed-form reconstruction of LMTD's implied exponential T(x)
# profile, derived from the log-mean assumption itself -- not something
# solve_it_LMTD computes. See CLAUDE.md.
def calc_lmtd_profile(lmtd_result, ops, geo, n_points: int = 200) -> ProfileData:
    Delta_theta_hot = ops.theta_c_i - lmtd_result.theta_a_o
    Delta_theta_cold = lmtd_result.theta_c_o - ops.theta_a_i

    theta_c_m = (ops.theta_c_i + lmtd_result.theta_c_o) / 2.0
    theta_a_m = (ops.theta_a_i + lmtd_result.theta_a_o) / 2.0
    coolant_state = get_fluid_properties(ops, theta_c_m, ops.p_c)
    air_state = get_air_properties(ops, theta_a_m, ops.p_a)
    C_dot_c = ops.m_dot_c * coolant_state.cp
    C_dot_a = ops.m_dot_a * air_state.cp

    x_frac = np.linspace(0.0, 1.0, n_points)
    ratio = Delta_theta_cold / Delta_theta_hot
    delta_T = Delta_theta_hot * ratio**x_frac
    Q = Delta_theta_hot * (1 - ratio**x_frac) / (1 / C_dot_c - 1 / C_dot_a)

    theta_c = ops.theta_c_i - Q / C_dot_c
    theta_a = lmtd_result.theta_a_o - Q / C_dot_a
    k_profile = np.full(n_points, lmtd_result.k)

    # Local dQ = k dA (T_B - T_A): Baehr, Thermodynamik, 16th ed. (2016), Sec. 3.1.
    A_tot = geo.A * geo.N_t * geo.N_r
    dQdL = lmtd_result.k * delta_T * (A_tot / geo.l)

    return ProfileData(x_frac=x_frac, theta_c=theta_c, theta_a=theta_a, k=k_profile, dQdL=dQdL)


# NTU's profile is a genuine reduction of its own converged element field
# (solve_it_NTU), along the same coolant path the field solver marches
# (cell_path_order). Same reduction as calc_cell_profile below; for NTU all
# tubes are identical, so mean/sum over tubes are exact copies/multiples.
def calc_ntu_profile(ntu_result, geo, n_elements: int) -> ProfileData:
    path = cell_path_order(geo.N_r, n_elements)
    n_points = len(path)

    theta_c = np.array([ntu_result.theta_c_grid[r, :, e].mean() for r, e in path])
    theta_a = np.array([ntu_result.theta_a_grid[r, :, e].mean() for r, e in path])
    k = np.array([ntu_result.k_grid[r, :, e].mean() for r, e in path])
    dQdL = np.array([ntu_result.dQdL_grid[r, :, e].sum() for r, e in path])   # extensive: sum over tubes

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, theta_c=theta_c, theta_a=theta_a, k=k, dQdL=dQdL)


# Cell's profile is a genuine reduction, not a reconstruction -- see
# solve_it_cell's own diagnostic pass for where this data actually comes from.
# path comes from solvers.cell_path_order -- the SAME function
# _relax_cell_grid's Pass 1 uses to drive its own loop, so this is provably
# the same path the solver itself walks, not a separate guess at it.
def calc_cell_profile(cell_result, geo, n_segments: int) -> ProfileData:
    path = cell_path_order(geo.N_r, n_segments)
    n_points = len(path)

    theta_c = np.array([cell_result.theta_c_grid[r, :, s].mean() for r, s in path])
    theta_a = np.array([cell_result.theta_a_grid[r, :, s].mean() for r, s in path])
    k = np.array([cell_result.k_grid[r, :, s].mean() for r, s in path])
    # dQdL is extensive (per-tube rates from parallel tubes add, not average) --
    # sum across tubes here, unlike the mean() used for the intensive quantities
    # above, so this matches LMTD/NTU's dQdL (both already scaled by N_t).
    dQdL = np.array([cell_result.dQdL_grid[r, :, s].sum() for r, s in path])

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, theta_c=theta_c, theta_a=theta_a, k=k, dQdL=dQdL)
