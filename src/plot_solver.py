"""
plot_solver.py
===============
Solver figures -- how the solvers behave, as opposed to what the cooler
does (plot_profiles.plot_cooler_results):
  plot_convergence        -- outlet temperatures and iteration errors per
                             iteration (output 'convergence'), per scenario
  plot_solver_performance -- caching benchmark (output 'caching'): the selected
                             solvers per cache mode on one axis, and each
                             solver's time per cache mode on its own axis
                             (linear page: split into time parts); drawn on
                             linear and on logarithmic axes
The resolution benchmark has its own figure (benchmark_solvers.plot_resolution).
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from src.run_scenario import draw_convergence
from src.solvers import SOLVER_FAMILY
from src.benchmark_solvers import (CACHE_MODES, PART_COLORS, TIME_PARTS, draw_cache_panel, draw_mode_comparison,
                                   resolution_labels)
from src.figures import save_figure, slug
from src.formatting import machine_info


def plot_convergence(result, ops, label: str, settings):
    """Returns the figure (does NOT call plt.show())."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 7.5))
    draw_convergence(ax1, ax2, result, ops, settings.central_omega, settings.cell_omega)
    resolutions = resolution_labels(settings)
    used = [f"{SOLVER_FAMILY[key]}: {resolutions[key]}" for key in result.results if resolutions[key]]
    fig.suptitle(f"Convergence - {label}\n" + (f"({', '.join(used)} per tube and row)" if used else ""),
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_figure(fig, 'solver', slug(f"Convergence - {label}"), f"Convergence - {label}")
    return fig


def plot_solver_performance(settings, cache_results, log: bool = False):
    """Returns the figure (does NOT call plt.show()). log=True: logarithmic time axes (the solvers span ~4 orders of magnitude). Top: the solvers side by side per cache mode. Bottom: each
    solver's run time per cache mode, split into CoolProp / lookups / k / rest (and, for NTU table, the build of its
    table -- counted in the no-cache and cold-cache run, assumed done in the warm-cache run)."""
    keys = list(cache_results)
    title = f"Caching performance ({'log' if log else 'linear'})"
    fig = plt.figure(figsize=(24, 13.5))
    resolutions = resolution_labels(settings)
    margins = dict(left=0.06, right=0.985)
    grid_top = fig.add_gridspec(1, 3, top=0.86, bottom=0.57, wspace=0.3, **margins)
    grid_bottom = fig.add_gridspec(1, len(keys), top=0.45, bottom=0.10, wspace=0.14, **margins)

    for i, mode in enumerate(CACHE_MODES):
        draw_mode_comparison(fig.add_subplot(grid_top[0, i]), cache_results, mode, resolutions, log)
    for i, key in enumerate(keys):
        draw_cache_panel(fig.add_subplot(grid_bottom[0, i]), key, cache_results,
                         resolution=resolutions[key], show_modes=(i == 0), log=log)

    parts = [] if log else [(part, color) for part, color in zip(TIME_PARTS, PART_COLORS)
                            if any(per_mode[m]['split'].get(part, 0.0) > 0 for per_mode in cache_results.values()
                                   for m in CACHE_MODES)]
    if parts:                                                 # the time-part colours exist on the linear page only
        fig.legend([Patch(facecolor=color, edgecolor='black') for _, color in parts], [part for part, _ in parts],
                   loc='lower center', ncol=len(parts), fontsize=11, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"Caching performance ({'log' if log else 'linear'})\n"
                 f"(top: solvers compared per cache mode; bottom: each solver on its own {'logarithmic' if log else 'linear'} time axis; "
                 f"{machine_info()})", fontsize=15, fontweight="bold", y=0.985)
    save_figure(fig, 'solver', slug(title), title)
    return fig
