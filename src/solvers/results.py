"""
results.py -- what the solvers return, and the registry of the five selectable solvers.
"""

from dataclasses import dataclass


@dataclass
class SolverResult:
    k: float
    Q_dot: float
    T_c_o: float
    T_a_o: float
    history_hot: list                # per-iteration change of Delta_T_hot [K]
    history_cold: list               # per-iteration change of Delta_T_cold [K]
    history_T_c_o: list
    history_T_a_o: list
    diagnostics: dict                # Pr/Re/Nu/alpha at the overall mean temperatures
    # (row, tube, segment) grids -- Cell and NTU only, LMTD leaves them None
    T_c_grid: object = None
    T_a_grid: object = None
    k_grid: object = None
    dQdL_grid: object = None
    local_diagnostics: dict = None   # Cell only: Re/Nu grids, not yet plotted
    solve_time: float = None         # wall time [s], set by solve_scenario


# The five selectable solvers: key -> (label, family). The family groups the variants that share an algorithm
# (the caching / resolution benchmarks compare families). Cell 2D = one representative tube, 3D = all tubes
# (identical results with uniform inlet air); NTU 'table' = P(NTU, R) lookup instead of the element field.
SOLVER_LABELS = {'lmtd': "LMTD", 'ntu_field': "NTU (field)", 'ntu_table': "NTU (table)",
                 'cell_2d': "Cell (2D)", 'cell_3d': "Cell (3D)"}
SOLVER_FAMILY = {'lmtd': "LMTD", 'ntu_field': "NTU", 'ntu_table': "NTU", 'cell_2d': "Cell", 'cell_3d': "Cell"}
DETAIL_ORDER = ('cell_3d', 'cell_2d', 'ntu_field', 'ntu_table', 'lmtd')    # most detailed first = the reference


@dataclass
class ScenarioResult:
    """The solvers' results on one operating point: {solver key: SolverResult}, in the order selected."""
    results: dict

    @property
    def reference_key(self) -> str:
        """The most detailed solver that ran (the deviations of the others are quoted against it, and it
        decides the precooling)."""
        return next(key for key in DETAIL_ORDER if key in self.results)

    @property
    def reference(self) -> SolverResult:
        return self.results[self.reference_key]
