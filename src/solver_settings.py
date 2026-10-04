"""
solver_settings.py
====================
Algorithm-tuning constants for the LMTD/NTU/Cell solvers (initial guesses,
convergence threshold, relaxation factor, iteration caps) -- distinct from
parameters.py, which holds physical/design inputs (geometry, operating
point) only. Raw values come from parameters.py; built the same way as
Geometry/OperatingConditions, via param_loader.from_parameters.
"""

from dataclasses import dataclass

from src.param_loader import from_parameters


@dataclass
class SolverSettings:
    # --- Shared by all three solvers --------------------------------------
    Delta_theta_hot_init: float          # initial Delta_theta_hot guess [K]
    Delta_theta_cold_init: float         # initial Delta_theta_cold guess [K]
    convergence_threshold: float   # convergence threshold on Delta_theta_hot/Delta_theta_cold change [K]
    central_omega: float           # under-relaxation factor of LMTD/NTU's shared outer loop
    cell_omega: float              # Cell's own relaxation factor (see parameters.CELL_OMEGA)
    cell_n_segments: int           # number of coolant-direction segments per tube pass (Cell only)
    ntu_n_elements: int            # number of elements per tube per row (NTU only)

    # --- Cell-only: relaxation (see solvers._relax_cell_grid) ---
    cell_threshold: float          # stop when no cell changes by more than this [K]
    cell_max_iter: int
    cell_min_iter: int

    # --- NTU-only: inner element-field iteration (see solvers.solve_ntu_field) ---
    ntu_field_threshold: float
    ntu_field_max_iter: int

    # --- Set from the run modes (CELL_2D / --cell-2d), not from parameters directly ---
    cell_2d: bool = False          # Cell: one representative tube instead of all


def get_solver_settings() -> SolverSettings:
    """Current solver tuning -- edit parameters.py to change the numbers."""
    return from_parameters(SolverSettings, {
        'Delta_theta_hot_init': 'DT_HOT_IT_INIT',
        'Delta_theta_cold_init': 'DT_COLD_IT_INIT',
        'convergence_threshold': 'CONVERGENCE_THRESHOLD',
        'central_omega': 'CENTRAL_OMEGA',
        'cell_omega': 'CELL_OMEGA',
        'cell_n_segments': 'CELL_N_SEGMENTS',
        'ntu_n_elements': 'NTU_N_ELEMENTS',
        'cell_threshold': 'CELL_THRESHOLD',
        'cell_max_iter': 'CELL_MAX_ITER',
        'cell_min_iter': 'CELL_MIN_ITER',
        'ntu_field_threshold': 'NTU_FIELD_THRESHOLD',
        'ntu_field_max_iter': 'NTU_FIELD_MAX_ITER',
    })
