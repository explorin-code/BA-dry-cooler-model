"""
plot_profiles.py
=================
Spatial-profile comparison: one panel per solver (LMTD/NTU/Cell), each
showing temperature, k, and local dQ/dL on three stacked y-axes, along the
coolant's flow path (inlet to outlet). Each axis type (temperature / k /
dQ-dL) shares the same limits across all three panels, so panels are
directly comparable despite belonging to different solvers. Cell's data is
genuinely computed; LMTD's/NTU's are post-hoc reconstructions -- see
analysis.py and CLAUDE.md.
"""

import numpy as np
import matplotlib.pyplot as plt

from src.solvers import ScenarioResult
from src.analysis import calc_lmtd_profile, calc_ntu_profile, calc_cell_profile
from src.run_scenario import COLOR_LMTD, COLOR_NTU, COLOR_CELL, lighten_color

COLOR_K = "#4d4d4d"       # dark gray -- neutral, not tied to hot/cold identity
COLOR_DQDL = "#1f77b4"    # blue -- distinct from the LMTD/NTU/Cell palette


def _data_range(*arrays):
    """Plain min/max across several arrays -- no padding; _nice_ticks below
    adds its own margin by rounding outward to nice numbers."""
    values = np.concatenate([np.asarray(a).ravel() for a in arrays])
    return float(values.min()), float(values.max())


_NICE_FRACTIONS = [1, 2, 2.5, 5, 10]


def _nice_step(raw_step):
    """Rounds raw_step up to the nearest 'nice' number (1/2/2.5/5/10 x
    10^n) -- so tick labels are whole/round numbers, not arbitrary fractions."""
    if raw_step <= 0:
        return 1.0
    exponent = np.floor(np.log10(raw_step))
    fraction = raw_step / 10**exponent
    for f in _NICE_FRACTIONS:
        if fraction <= f:
            return f * 10**exponent
    return 10 * 10**exponent


def _next_nice_step(step):
    """The next larger step in the same 1/2/2.5/5/10 x 10^n sequence."""
    exponent = np.floor(np.log10(step))
    fraction = round(step / 10**exponent, 6)
    idx = _NICE_FRACTIONS.index(fraction) if fraction in _NICE_FRACTIONS else 0
    if idx + 1 < len(_NICE_FRACTIONS):
        return _NICE_FRACTIONS[idx + 1] * 10**exponent
    return _NICE_FRACTIONS[0] * 10**(exponent + 1)


def _nice_ticks(data_lo, data_hi, n_intervals: int = 4):
    """Exactly n_intervals+1 evenly spaced NICE (round-number) ticks that
    fully cover [data_lo, data_hi]. Fixing n_intervals the same for every
    axis in a panel is what makes unrelated quantities (temperature, k,
    dQ/dL) land their own round-number ticks on the same horizontal lines
    -- the counts match even though the actual numbers don't."""
    span = data_hi - data_lo
    if span <= 0:
        span = abs(data_hi) if data_hi != 0 else 1.0
    step = _nice_step(span / n_intervals)
    for _ in range(4):   # a couple of nice-step bumps is always enough in practice
        nice_lo = np.floor(data_lo / step) * step
        ticks = nice_lo + step * np.arange(n_intervals + 1)
        if ticks[-1] >= data_hi - 1e-9:
            return ticks
        step = _next_nice_step(step)
    return ticks


def _plot_solver_panel(ax_T, profile, color, T_ticks, k_ticks, dQdL_ticks, title):
    """One solver's panel: temperature on the primary axis, k and dQ/dL on
    two offset twin axes. Each axis's ylim is set to exactly span its own
    nice ticks, so the same number of nice ticks on every axis lines up
    across all three, even though the actual values differ per axis."""
    color_air = lighten_color(color, 0.55)

    ax_T.plot(profile.x_frac, profile.T_coolant, color=color, linewidth=2.2, label="Coolant T")
    ax_T.plot(profile.x_frac, profile.T_air, color=color_air, linewidth=2.2, linestyle='--', label="Air T")
    ax_T.set_ylabel("Temperature [°C]", color=color)
    ax_T.set_ylim(T_ticks[0], T_ticks[-1])
    ax_T.set_yticks(T_ticks)
    ax_T.tick_params(axis='y', labelcolor=color)

    ax_k = ax_T.twinx()
    ax_k.grid(False)   # only ax_T should draw gridlines -- twin axes' own grid was rendering on top of the legend
    ax_k.plot(profile.x_frac, profile.k, color=COLOR_K, linewidth=1.6, linestyle=':', label="k")
    ax_k.set_ylabel("k [W/m²K]", color=COLOR_K)
    ax_k.set_ylim(k_ticks[0], k_ticks[-1])
    ax_k.set_yticks(k_ticks)
    ax_k.tick_params(axis='y', labelcolor=COLOR_K)

    ax_dQdL = ax_T.twinx()
    ax_dQdL.grid(False)
    ax_dQdL.spines['right'].set_position(('outward', 60))
    ax_dQdL.plot(profile.x_frac, profile.dQdL, color=COLOR_DQDL, linewidth=1.6, linestyle='-.', label="dQ/dL")
    ax_dQdL.set_ylabel("dQ/dL [W/m]", color=COLOR_DQDL)
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


def plot_profiles(result: ScenarioResult, ops, geo, label: str, n_segments: int):
    lmtd_p = calc_lmtd_profile(result.lmtd, ops, geo)
    ntu_p = calc_ntu_profile(result.ntu, ops, geo)
    cell_p = calc_cell_profile(result.cell, geo, n_segments)

    T_ticks = _nice_ticks(*_data_range(lmtd_p.T_coolant, lmtd_p.T_air, ntu_p.T_coolant, ntu_p.T_air,
                                        cell_p.T_coolant, cell_p.T_air))
    k_ticks = _nice_ticks(*_data_range(lmtd_p.k, ntu_p.k, cell_p.k))
    dQdL_ticks = _nice_ticks(*_data_range(lmtd_p.dQdL, ntu_p.dQdL, cell_p.dQdL))

    fig, (ax_lmtd, ax_ntu, ax_cell) = plt.subplots(3, 1, figsize=(11, 13))
    fig.canvas.manager.set_window_title(f"{label} -- Profiles")

    _plot_solver_panel(ax_lmtd, lmtd_p, COLOR_LMTD, T_ticks, k_ticks, dQdL_ticks, "LMTD")
    _plot_solver_panel(ax_ntu, ntu_p, COLOR_NTU, T_ticks, k_ticks, dQdL_ticks, "NTU")
    _plot_solver_panel(ax_cell, cell_p, COLOR_CELL, T_ticks, k_ticks, dQdL_ticks, "Cell")

    ax_cell.set_xlabel("Coolant flow path (0 = inlet, 1 = outlet)")

    fig.suptitle(f"Spatial Profile Comparison — {label}\n"
                 f"(LMTD/NTU profiles are reconstructions, not directly solved — see CLAUDE.md)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.88, 0.93])
    return fig
