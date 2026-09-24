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
    dT_hot_it_init: float          # initial dT_hot guess [K]
    dT_cold_it_init: float         # initial dT_cold guess [K]
    dT_kick: float                 # forced initial dT_hot/dT_cold to kickstart the while loop [K]
    convergence_threshold: float   # convergence threshold on dT_hot/dT_cold change [K]
    central_omega: float           # under-relaxation factor used for all three solvers
    cell_n_segments: int           # number of coolant-direction segments per tube pass (Cell only)

    # --- Cell-only: two-stage relaxation (see solvers._relax_cell_grid) ---
    cell_stage1_threshold: float
    cell_stage1_max_iter: int
    cell_stage1_min_iter: int
    cell_stage2_threshold: float
    cell_stage2_max_iter: int
    cell_stage2_min_iter: int


def get_solver_settings() -> SolverSettings:
    """Current solver tuning -- edit parameters.py to change the numbers."""
    return from_parameters(SolverSettings, {
        'dT_hot_it_init': 'DT_HOT_IT_INIT',
        'dT_cold_it_init': 'DT_COLD_IT_INIT',
        'dT_kick': 'DT_KICK',
        'convergence_threshold': 'CONVERGENCE_THRESHOLD',
        'central_omega': 'CENTRAL_OMEGA',
        'cell_n_segments': 'CELL_N_SEGMENTS',
        'cell_stage1_threshold': 'CELL_STAGE1_THRESHOLD',
        'cell_stage1_max_iter': 'CELL_STAGE1_MAX_ITER',
        'cell_stage1_min_iter': 'CELL_STAGE1_MIN_ITER',
        'cell_stage2_threshold': 'CELL_STAGE2_THRESHOLD',
        'cell_stage2_max_iter': 'CELL_STAGE2_MAX_ITER',
        'cell_stage2_min_iter': 'CELL_STAGE2_MIN_ITER',
    })
