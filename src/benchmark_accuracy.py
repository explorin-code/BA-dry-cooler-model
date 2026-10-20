"""
benchmark_accuracy.py
======================
Accuracy sweep (solver output 'accuracy'): how far the simpler solvers are from the most detailed one of the selection
(the reference, as everywhere: ScenarioResult.reference_key) when the water inlet temperature moves away from the run's
own (the property variation along the tube grows with it) -- at the run's geometry and ambient air point, the deviations
in Q [%] and in the pinch T_c,o - T_a,i [K].
"""

import logging

import matplotlib.pyplot as plt

from src import solvers
from src.figures import save_figure, SOLVER
from src.operating_conditions import ops_like
from src.plot_style import PRINT_DPI, PRINT_LW as LW, PRINT_WIDTH, TUM_GREY_50, print_context, solver_color
from src.quiet import quiet_solvers

logger = logging.getLogger(__name__)



def run_accuracy_sweep(ops, geo, settings, modes):
    """{'x': water inlet temperatures [°C] (modes.accuracy_t_coolant), 'x_run': the run's own, 'reference': solver key,
    'rows': {solver key: [(Q [W], T_c_o [°C], pinch [K]), ...]}}, one solver per algorithm family of the selection."""
    keys = solvers.scenario_keys(modes.solvers)
    if len(keys) < 2:
        raise ValueError("The accuracy sweep compares solvers: select at least two.")
    rows = {key: [] for key in keys}
    with quiet_solvers():
        for T_c_i in modes.accuracy_t_coolant:
            point = ops_like(geo, ops, T_c_i=T_c_i)
            for key in keys:
                r = solvers.run_solver(key, point, geo, settings, profile=False)
                rows[key].append((r.Q_dot, r.T_c_o, r.T_c_o - point.T_a_i))
    reference = next(key for key in solvers.DETAIL_ORDER if key in keys)
    return {'x': list(modes.accuracy_t_coolant), 'x_run': ops.T_c_i, 'reference': reference, 'rows': rows}


def deviations(sweep: dict, key: str):
    """(Q deviation [%], pinch deviation [K]) of solver `key` from the sweep's reference, per water temperature."""
    ref, rows = sweep['rows'][sweep['reference']], sweep['rows'][key]
    return ([(q - q0) / q0 * 100.0 for (q, _, _), (q0, _, _) in zip(rows, ref)],
            [p - p0 for (_, _, p), (_, _, p0) in zip(rows, ref)])


PANELS = {'q': (0, r"$\dot{Q}$ deviation [%]"), 'pinch': (1, r"$\Delta T_{pinch}$ deviation [K]")}


def plot_accuracy_sweep(sweep: dict, panels=('q', 'pinch'), name: str = 'accuracy_sweep'):
    """Graphs at print size (one per entry of `panels`: Q deviation [%] and / or pinch deviation [K]) of every non-reference solver from the reference
    against the water inlet temperature, one legend (in the first graph). The dotted line marks the run's own water temperature.
    panels=('pinch',): the pinch graph alone (e.g. for slides), saved under `name`."""
    reference = solvers.SOLVER_LABELS[sweep['reference']]
    with print_context(talk=True):
        fig, axes = plt.subplots(1, len(panels), figsize=(PRINT_WIDTH if len(panels) > 1 else PRINT_WIDTH * 0.62, 2.3 if len(panels) > 1 else 2.6), squeeze=False)
        axes = axes[0]
        for key in sweep['rows']:
            if key == sweep['reference']:
                continue
            deviation = deviations(sweep, key)
            for ax, panel in zip(axes, panels):
                ax.plot(sweep['x'], deviation[PANELS[panel][0]], marker='o', ms=2.5, color=solver_color(key), lw=LW, label=solvers.SOLVER_LABELS[key])
        for ax, panel in zip(axes, panels):
            ax.axhline(0.0, color=solver_color(sweep['reference']), lw=LW, label=f"{reference} (reference)")
            ax.axvline(sweep['x_run'], color=TUM_GREY_50, lw=0.9, ls=':', zorder=0, label=f"run: {sweep['x_run']:.0f} °C")
            ax.set_xlabel(r"$T_{c,i}$ [°C]")
            ax.set_ylabel(PANELS[panel][1])
        if panels[0] == 'q':
            low, high = axes[0].get_ylim()
            axes[0].set_ylim(low, high + 0.3 * (high - low))                  # room for the one legend in the upper left corner
            axes[0].legend(loc='upper left', ncol=2)
        else:
            axes[0].legend(loc='lower left')                                 # the pinch deviation falls away from the upper left: the lower left is free
        fig.tight_layout(w_pad=1.2)
        save_figure(fig, SOLVER, name, "Accuracy sweep", dpi=PRINT_DPI)
    return fig


def summary_lines(sweep: dict) -> list:
    """Largest deviation of every solver, for the log."""
    lines = []
    for key in sweep['rows']:
        if key != sweep['reference']:
            q, p = deviations(sweep, key)
            lines.append(f"{solvers.SOLVER_LABELS[key]} vs. {solvers.SOLVER_LABELS[sweep['reference']]}: "
                         f"|dQ| <= {max(map(abs, q)):.2f} %, |d pinch| <= {max(map(abs, p)):.2f} K")
    return lines
