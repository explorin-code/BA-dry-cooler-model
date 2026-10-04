"""
plot_profiles.py
=================
Cooler-results figure (plot_cooler_results): what the cooler does, per
scenario -- input conditions and economics on top, then one row per solver
with its result box next to its spatial profile (temperature, k and local
dQ/dL on three y-axes, along the coolant's flow path, inlet to outlet). Each axis type (temperature / k /
dQ-dL) shares the same limits across all three panels, so panels are
directly comparable despite belonging to different solvers. Cell's and
NTU's data is genuinely computed; LMTD's is a post-hoc reconstruction --
see analysis.py and CLAUDE.md.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.solvers import ScenarioResult
from src.analysis import calc_lmtd_profile, calc_ntu_profile, calc_cell_profile
from src.plot_style import COLOR_LMTD, COLOR_NTU, COLOR_CELL, lighten_color, _nice_ticks, draw_text_grid
from src.run_scenario import input_condition_cells, economics_cells, output_condition_cells

COLOR_K = "#4d4d4d"       # dark gray -- neutral, not tied to hot/cold identity
COLOR_DQDL = "#1f77b4"    # blue -- distinct from the LMTD/NTU/Cell palette


def _data_range(*arrays):
    """Plain min/max across several arrays -- no padding; _nice_ticks
    (plot_style.py) adds its own margin by rounding outward to nice numbers."""
    values = np.concatenate([np.asarray(a).ravel() for a in arrays])
    return float(values.min()), float(values.max())


def _plot_solver_panel(ax_T, profile, color, T_ticks, k_ticks, dQdL_ticks, title):
    """One solver's panel: temperature on the primary axis, k and dQ/dL on
    two offset twin axes. Each axis's ylim is set to exactly span its own
    nice ticks, so the same number of nice ticks on every axis lines up
    across all three, even though the actual values differ per axis."""
    color_air = lighten_color(color, 0.55)

    ax_T.plot(profile.x_frac, profile.theta_c, color=color, linewidth=2.2, label=r"$T_c$ (coolant)")
    ax_T.plot(profile.x_frac, profile.theta_a, color=color_air, linewidth=2.2, linestyle='--', label=r"$T_a$ (air)")
    ax_T.set_ylabel("Temperature [°C]", color=color)
    ax_T.set_ylim(T_ticks[0], T_ticks[-1])
    ax_T.set_yticks(T_ticks)
    ax_T.tick_params(axis='y', labelcolor=color)

    ax_k = ax_T.twinx()
    ax_k.grid(False)   # only ax_T should draw gridlines -- twin axes' own grid was rendering on top of the legend
    ax_k.plot(profile.x_frac, profile.k, color=COLOR_K, linewidth=1.6, linestyle=':', label=r"$k$")
    ax_k.set_ylabel(r"$k$ [W/m²K]", color=COLOR_K)
    ax_k.set_ylim(k_ticks[0], k_ticks[-1])
    ax_k.set_yticks(k_ticks)
    ax_k.tick_params(axis='y', labelcolor=COLOR_K)

    ax_dQdL = ax_T.twinx()
    ax_dQdL.grid(False)
    ax_dQdL.spines['right'].set_position(('outward', 60))
    ax_dQdL.plot(profile.x_frac, profile.dQdL, color=COLOR_DQDL, linewidth=1.6, linestyle='-.', label=r"$dQ/dL$")
    ax_dQdL.set_ylabel(r"$dQ/dL$ [W/m]", color=COLOR_DQDL)
    ax_dQdL.set_ylim(dQdL_ticks[0], dQdL_ticks[-1])
    ax_dQdL.set_yticks(dQdL_ticks)
    ax_dQdL.tick_params(axis='y', labelcolor=COLOR_DQDL)

    ax_T.set_title(title, color=color, fontweight="bold")

    # Legend goes on ax_dQdL, not ax_T -- ax_dQdL was created last, so its
    # render pass happens last and draws on top of both other axes' lines
    # (a twin axes' whole render pass, not just its gridlines, paints over
    # earlier axes -- attaching the legend to ax_T left it behind ax_k's line).
    lines = ax_T.get_lines() + ax_k.get_lines() + ax_dQdL.get_lines()
    labels = [l.get_label() for l in lines]
    ax_dQdL.legend(lines, labels, fontsize=8, loc='best',
                    framealpha=1.0, facecolor='white', edgecolor='gray')


BOX_FACECOLORS = {'LMTD': "#eaf5ec", 'NTU': "#f3ecf5", 'Cell': "#fdf0e0"}


def plot_cooler_results(result: ScenarioResult, ops, geo, label: str, n_segments: int, n_elements: int,
                        W_pump=None, W_fan=None, m_dot_ev=None):
    """Returns the figure (does NOT call plt.show()). LMTD/NTU boxes also
    show their Q and pinch deviation from Cell. W_pump/W_fan/m_dot_ev: economics,
    None for whichever aren't computed."""
    profiles = {
        'LMTD': calc_lmtd_profile(result.lmtd, ops, geo),
        'NTU': calc_ntu_profile(result.ntu, geo, n_elements),
        'Cell': calc_cell_profile(result.cell, geo, n_segments),
    }
    solver_results = {'LMTD': result.lmtd, 'NTU': result.ntu, 'Cell': result.cell}
    colors = {'LMTD': COLOR_LMTD, 'NTU': COLOR_NTU, 'Cell': COLOR_CELL}

    T_ticks = _nice_ticks(*_data_range(*[a for p in profiles.values() for a in (p.theta_c, p.theta_a)]))
    k_ticks = _nice_ticks(*_data_range(*[p.k for p in profiles.values()]))
    dQdL_ticks = _nice_ticks(*_data_range(*[p.dQdL for p in profiles.values()]))

    fig = plt.figure(figsize=(21, 14))
    fig.canvas.manager.set_window_title(f"{label} -- Cooler results")
    grid = GridSpec(4, 2, figure=fig, height_ratios=[0.6, 1, 1, 1], width_ratios=[1.2, 1.45])

    # --- Header: inputs + geometry, economics below ----------------------
    ax_header = fig.add_subplot(grid[0, :])
    ax_header.axis("off")
    draw_text_grid(ax_header, 0.5, 0.66, input_condition_cells(ops, geo), "whitesmoke", "gray")
    draw_text_grid(ax_header, 0.5, 0.08, economics_cells(W_pump, W_fan, m_dot_ev), "#e8eef5", "#4a6fa5")

    # --- One row per solver: result box | profile -------------------------
    for row, name in enumerate(('LMTD', 'NTU', 'Cell'), start=1):
        ax_text = fig.add_subplot(grid[row, 0])
        ax_text.axis("off")
        reference = None if name == 'Cell' else result.cell
        draw_text_grid(ax_text, 0.5, 0.5, output_condition_cells(name, solver_results[name], ops.theta_a_i,
                                                                 reference=reference),
                       BOX_FACECOLORS[name], colors[name])
        ax_profile = fig.add_subplot(grid[row, 1])
        _plot_solver_panel(ax_profile, profiles[name], colors[name], T_ticks, k_ticks, dQdL_ticks, name)
        if name == 'Cell':
            ax_profile.set_xlabel("Coolant flow path (0 = inlet, 1 = outlet)")

    fig.suptitle(f"Cooler Results — {label}", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.95, 0.95])
    return fig
