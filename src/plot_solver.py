"""
plot_solver.py
===============
Solver figures -- how the solvers behave, as opposed to what the cooler
does (plot_profiles.plot_cooler_results):
  plot_convergence        -- outlet temperatures and iteration errors per
                             iteration (PLOT_CONVERGENCE), per scenario
  plot_solver_performance -- caching benchmark (BENCHMARK_MODE): all three
                             solvers per cache mode on one linear axis, their
                             own work without lookups and k (rest), and each
                             solver's time split per cache mode on its own
                             linear axis
The resolution benchmark has its own figure (benchmark_solvers.plot_resolution).
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.run_scenario import draw_convergence
from src.benchmark_solvers import (CACHE_MODES, draw_algorithm_comparison, draw_cache_panel,
                                   draw_mode_comparison, machine_info, resolution_labels)


def plot_convergence(result, ops, label: str, settings):
    """Returns the figure (does NOT call plt.show())."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 7.5))
    fig.canvas.manager.set_window_title(f"{label} -- Convergence")
    draw_convergence(ax1, ax2, result, ops, settings.central_omega, settings.cell_omega)
    resolutions = resolution_labels(settings)
    fig.suptitle(f"Solver Convergence — {label}\n"
                 f"(NTU: {resolutions['NTU']}, Cell: {resolutions['Cell']} per tube and row)",
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


def plot_solver_performance(label: str, settings, cache_results):
    """Returns the figure (does NOT call plt.show()). Top: the three solvers
    side by side per cache mode. Bottom: each solver's run time per cache
    mode, split into CoolProp / lookups / k / rest."""
    fig = plt.figure(figsize=(24, 13))
    fig.canvas.manager.set_window_title(f"{label} -- Solver performance (caching)")
    grid = GridSpec(2, 12, figure=fig, height_ratios=[1, 1])   # top: 4 panels x 3 cols, bottom: 3 x 4
    resolutions = resolution_labels(settings)

    for i, mode in enumerate(CACHE_MODES):
        draw_mode_comparison(fig.add_subplot(grid[0, 3 * i:3 * i + 3]), cache_results, mode, resolutions)
    draw_algorithm_comparison(fig.add_subplot(grid[0, 9:]), cache_results, resolutions)
    for i, name in enumerate(('LMTD', 'NTU', 'Cell')):
        draw_cache_panel(fig.add_subplot(grid[1, 4 * i:4 * i + 4]), name, cache_results,
                         show_legend=(name == 'Cell'), resolution=resolutions[name])

    fig.suptitle(f"Solver Performance: Property Caching — {label}\n"
                 f"(top: solvers compared per cache mode; bottom: each solver on its own linear time axis; "
                 f"{machine_info()})", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig
