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
    T_coolant: object
    T_air: object
    k: object
    dQdL: object


# AI-REVIEW: closed-form reconstruction of LMTD's implied exponential T(x)
# profile, derived from the log-mean assumption itself -- not something
# solve_it_LMTD computes. See CLAUDE.md.
def calc_lmtd_profile(lmtd_result, ops, geo, n_points: int = 200) -> ProfileData:
    dT_hot = ops.T_coolant_in - lmtd_result.T_air_out
    dT_cold = lmtd_result.T_coolant_out - ops.T_air_in

    T_coolant_mean = (ops.T_coolant_in + lmtd_result.T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + lmtd_result.T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
    air_state = get_air_properties(ops, T_air_mean, ops.P_air)
    W_coolant = ops.m_dot_1 * coolant_state.cp
    W_air = ops.m_dot_2 * air_state.cp

    x_frac = np.linspace(0.0, 1.0, n_points)
    ratio = dT_cold / dT_hot
    delta_T = dT_hot * ratio**x_frac
    Q = dT_hot * (1 - ratio**x_frac) / (1 / W_coolant - 1 / W_air)

    T_coolant = ops.T_coolant_in - Q / W_coolant
    T_air = lmtd_result.T_air_out - Q / W_air
    k_profile = np.full(n_points, lmtd_result.k)

    A_total = geo.A * geo.n_tubes * geo.n_rows
    dQdL = lmtd_result.k * delta_T * (A_total / geo.l)

    return ProfileData(x_frac=x_frac, T_coolant=T_coolant, T_air=T_air, k=k_profile, dQdL=dQdL)


# NTU's profile is a genuine reduction of its own converged element field
# (solve_it_NTU), along the same coolant path the field solver marches
# (cell_path_order). Same reduction as calc_cell_profile below; for NTU all
# tubes are identical, so mean/sum over tubes are exact copies/multiples.
def calc_ntu_profile(ntu_result, geo, n_elements: int) -> ProfileData:
    path = cell_path_order(geo.n_rows, n_elements)
    n_points = len(path)

    T_coolant = np.array([ntu_result.T_c_grid[r, :, e].mean() for r, e in path])
    T_air = np.array([ntu_result.T_a_grid[r, :, e].mean() for r, e in path])
    k = np.array([ntu_result.k_grid[r, :, e].mean() for r, e in path])
    dQdL = np.array([ntu_result.dQdL_grid[r, :, e].sum() for r, e in path])   # extensive: sum over tubes

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, T_coolant=T_coolant, T_air=T_air, k=k, dQdL=dQdL)


# Cell's profile is a genuine reduction, not a reconstruction -- see
# solve_it_cell's own diagnostic pass for where this data actually comes from.
# path comes from solvers.cell_path_order -- the SAME function
# _relax_cell_grid's Pass 1 uses to drive its own loop, so this is provably
# the same path the solver itself walks, not a separate guess at it.
def calc_cell_profile(cell_result, geo, n_segments: int) -> ProfileData:
    path = cell_path_order(geo.n_rows, n_segments)
    n_points = len(path)

    T_coolant = np.array([cell_result.T_c_grid[r, :, s].mean() for r, s in path])
    T_air = np.array([cell_result.T_a_grid[r, :, s].mean() for r, s in path])
    k = np.array([cell_result.k_grid[r, :, s].mean() for r, s in path])
    # dQdL is extensive (per-tube rates from parallel tubes add, not average) --
    # sum across tubes here, unlike the mean() used for the intensive quantities
    # above, so this matches LMTD/NTU's dQdL (both already scaled by n_tubes).
    dQdL = np.array([cell_result.dQdL_grid[r, :, s].sum() for r, s in path])

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, T_coolant=T_coolant, T_air=T_air, k=k, dQdL=dQdL)
