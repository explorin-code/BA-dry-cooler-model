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
    Delta_T_hot_init: float          # initial Delta_T_hot guess [K]
    Delta_T_cold_init: float         # initial Delta_T_cold guess [K]
    central_omega: float           # under-relaxation factor of LMTD's outer loop
    outer_max_iter: int            # iteration cap of LMTD/NTU's outer loop
    ntu_omega: float               # NTU's own outer-loop relaxation factor (see parameters.NTU_OMEGA)
    cell_omega: float              # Cell's own relaxation factor (see parameters.CELL_OMEGA)
    cell_n_segments: int           # number of coolant-direction segments per tube pass (Cell only)
    ntu_n_elements: int            # number of elements per tube per row (NTU only)
    ntu_mode: str                  # 'field' or 'table' (see solvers.solve_it_NTU)

    # --- Convergence (all three solvers; per-solver overrides only as
    # hardcoded constants in solvers.py, for tests) ---
    convergence_threshold: float   # both outlet temperatures change < this per iteration [K]
    ntu_field_max_iter: int
    cell_max_iter: int

    # --- Set from the run modes (CELL_2D / --cell-2d), not from parameters directly ---
    cell_2d: bool = False          # Cell: one representative tube instead of all


def get_solver_settings() -> SolverSettings:
    """Current solver tuning -- edit parameters.py to change the numbers."""
    return from_parameters(SolverSettings, {
        'Delta_T_hot_init': 'DT_HOT_IT_INIT',
        'Delta_T_cold_init': 'DT_COLD_IT_INIT',
        'convergence_threshold': 'CONVERGENCE_THRESHOLD',
        'central_omega': 'CENTRAL_OMEGA',
        'outer_max_iter': 'OUTER_MAX_ITER',
        'ntu_omega': 'NTU_OMEGA',
        'cell_omega': 'CELL_OMEGA',
        'cell_n_segments': 'CELL_N_SEGMENTS',
        'ntu_n_elements': 'NTU_N_ELEMENTS',
        'ntu_mode': 'NTU_MODE',
        'cell_max_iter': 'CELL_MAX_ITER',
        'ntu_field_max_iter': 'NTU_FIELD_MAX_ITER',
    })
