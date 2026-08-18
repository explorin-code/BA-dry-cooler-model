"""
run_scenario.py
================
Builds and returns the convergence figure for an already-solved
ScenarioResult (see solvers.solve_scenario()). Pure display -- no solving.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import seaborn as sns

from src.solvers import ScenarioResult
from src.operating_conditions import OperatingConditions
from src.economics import calc_total_power

sns.set_theme(style="whitegrid", context="talk")

COLOR_LMTD = "#1b7837"  # green
COLOR_NTU = "#762a83"   # purple
COLOR_CELL = "#e08214"  # orange


def lighten_color(color, factor=0.5):
    """Blend `color` toward white. factor=0 -> unchanged, factor=1 -> white."""
    r, g, b = mc.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


def plot_with_tail(ax, series, max_len, color, label, linewidth=2.2, marker='o', markersize=3):
    """
    Plot `series` as a solid line for the iterations it actually ran.
    If it converged before `max_len`, extend it as a flat dashed line
    (starting from its last value) so all solvers span the same x-range.
    """
    n = len(series)
    ax.plot(range(n), series, color=color, linewidth=linewidth,
             marker=marker, markersize=markersize, label=label, zorder=3)
    if n < max_len:
        ax.plot(range(n - 1, max_len), [series[-1]] * (max_len - n + 1),
                 color=color, linewidth=linewidth, linestyle='--', alpha=0.55, zorder=2)


def format_input_conditions(ops: OperatingConditions) -> str:
    """
    Build the 'Input Conditions' string: inlet temps for both media plus
    all three flow-rate representations (velocity / mass / volume flow)
    for both media, regardless of which one was actually specified.
    """
    coolant_line = (
        f"Coolant ({ops.coolant_type}):  "
        f"T_in = {ops.T_coolant_in:5.1f} °C   "
        f"w = {ops.w_coolant:6.3f} m/s   "
        f"ṁ = {ops.m_coolant:6.3f} kg/s   "
        f"V̇ = {ops.V_coolant:7.5f} m³/s"
    )
    air_line = (
        f"Air (phi = {ops.phi_air:5.1f}):        "
        f"T_in = {ops.T_air_in:5.1f} °C   "
        f"w = {ops.w_o:6.3f} m/s   "
        f"ṁ = {ops.m_o:6.3f} kg/s   "
        f"V̇ = {ops.V_o:7.5f} m³/s"
    )
    return coolant_line + "\n" + air_line


def format_geometry_info(geo) -> str:
    """
    Build the geometric-info row for the input box: number of tube rows,
    number of tubes per row, and the frontal (inflow) area expressed as
    height x width = area.
    """
    return (
        f"Geometry:      "
        f"n_rows = {geo.n_rows:3d}   "
        f"n_tubes = {geo.n_tubes:3d}   "
        f"{geo.height:.3f} m x {geo.width:.3f} m = {geo.inflow_cross_section:.4f} m²"
    )


def format_output_conditions(label: str, T_coolant_out: float, T_air_out: float,
                              dQ: float, diagnostics: dict) -> str:
    """
    Build a per-solver 'Output' string: outlet temperatures, heat transfer
    rate, and the final-iteration dimensionless groups (Pr/Re/Nu) for both
    the coolant and air side.
    """
    header = f"{label} — Results"
    temps_line = (
        f"T_coolant_out = {T_coolant_out:5.2f} °C   "
        f"T_air_out = {T_air_out:5.2f} °C   "
        f"Q = {dQ/1000:6.2f} kW"
    )
    coolant_line = (
        f"Coolant:  Pr = {diagnostics['Pr_coolant']:6.3f}   "
        f"Re = {diagnostics['Re_coolant']:8.1f}   "
        f"Nu = {diagnostics['Nu_coolant']:7.2f}   "
        f"α_i = {diagnostics['alpha_i']:7.1f} W/m²K"
    )
    air_line = (
        f"Air:      Pr = {diagnostics['Pr_air']:6.3f}   "
        f"Re = {diagnostics['Re_air']:8.1f}   "
        f"Nu = {diagnostics['Nu_air']:7.2f}   "
        f"α_R = {diagnostics['alpha_R']:7.1f} W/m²K"
    )
    return header + "\n" + temps_line + "\n" + coolant_line + "\n" + air_line


def format_economics(W_pump, W_fan, m_water) -> str:
    """Build the 'Economics' box: pump/fan/total power and water usage.
    Any value left as None (not yet available) prints as 'n/a'."""
    pump_str = f"{W_pump:7.2f} W" if W_pump is not None else "    n/a"
    fan_str = f"{W_fan:7.2f} W" if W_fan is not None else "    n/a"
    W_total = calc_total_power(W_pump, W_fan)
    total_str = f"{W_total:7.2f} W" if W_total is not None else "    n/a"
    water_str = f"{m_water * 1000:6.3f} g/s" if m_water is not None else "   n/a"

    return (
        f"Economics — Pump: {pump_str}   Fan: {fan_str}   "
        f"Total: {total_str}   Water: {water_str}"
    )


def plot_scenario(result: ScenarioResult, ops: OperatingConditions, geo, label: str, omega: float,
                   W_pump=None, W_fan=None, m_water=None):
    """Builds the convergence figure for an already-solved ScenarioResult
    and returns it (does NOT call plt.show() -- the caller decides when to
    display, so multiple scenarios' windows can be shown together).
    W_pump/W_fan/m_water are optional economics figures to display -- pass
    None for whichever aren't computed yet."""
    lmtd, ntu, cell = result.lmtd, result.ntu, result.cell

    print(f"[{label}] [LMTD] converged! k = {lmtd.k:.2f} W/m2K, Q = {lmtd.dQ/1000:.2f} kW, "
          f"iterations = {len(lmtd.history_hot)}")
    print(f"[{label}] [NTU]  converged! k = {ntu.k:.2f} W/m2K, Q = {ntu.dQ/1000:.2f} kW, "
          f"iterations = {len(ntu.history_hot)}")
    print(f"[{label}] [Cell] converged! k = {cell.k:.2f} W/m2K, Q = {cell.dQ/1000:.2f} kW, "
          f"iterations = {len(cell.history_hot)}")

    # --- Input conditions actually used by the solvers -------------------
    input_conditions_text = format_input_conditions(ops) + "\n" + format_geometry_info(geo)

    # --- Per-solver output fields ------------------------------------------
    output_text_lmtd = format_output_conditions("LMTD", lmtd.T_coolant_out, lmtd.T_air_out, lmtd.dQ, lmtd.diagnostics)
    output_text_ntu = format_output_conditions("NTU", ntu.T_coolant_out, ntu.T_air_out, ntu.dQ, ntu.diagnostics)
    output_text_cell = format_output_conditions("Cell", cell.T_coolant_out, cell.T_air_out, cell.dQ, cell.diagnostics)
    economics_text = format_economics(W_pump, W_fan, m_water)

    # --- Colors: one hue per solver, air = lighter tone of the coolant hue ---
    color_lmtd_air = lighten_color(COLOR_LMTD, 0.55)
    color_ntu_air = lighten_color(COLOR_NTU, 0.55)
    color_cell_air = lighten_color(COLOR_CELL, 0.55)

    color_lmtd_hoterr = lighten_color(COLOR_LMTD, 0.45)
    color_ntu_hoterr = lighten_color(COLOR_NTU, 0.45)
    color_cell_hoterr = lighten_color(COLOR_CELL, 0.45)

    max_len_temp = max(len(lmtd.history_T_coolant), len(ntu.history_T_coolant), len(cell.history_T_coolant))
    max_len_err = max(len(lmtd.history_hot), len(ntu.history_hot), len(cell.history_hot))

    # ====================================================================
    # All three plots share ONE figure/window, side by side
    # ====================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 8.5))
    fig.canvas.manager.set_window_title(label)

    # --- Left: outlet temperatures -------------------------------------
    plot_with_tail(ax1, lmtd.history_T_coolant, max_len_temp, COLOR_LMTD, "LMTD - Coolant out")
    plot_with_tail(ax1, lmtd.history_T_air, max_len_temp, color_lmtd_air, "LMTD - Air out")
    plot_with_tail(ax1, ntu.history_T_coolant, max_len_temp, COLOR_NTU, "NTU - Coolant out")
    plot_with_tail(ax1, ntu.history_T_air, max_len_temp, color_ntu_air, "NTU - Air out")
    plot_with_tail(ax1, cell.history_T_coolant, max_len_temp, COLOR_CELL, "Cell - Coolant out")
    plot_with_tail(ax1, cell.history_T_air, max_len_temp, color_cell_air, "Cell - Air out")

    ax1.axhline(ops.T_air_in, color='blue', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label="T_air_in")
    ax1.axhline(ops.T_coolant_in, color='red', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label="T_coolant_in")

    ax1.set_title(f"Outlet Temperatures ($\\omega$ = {omega})")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Temperature [°C]")
    ax1.legend(fontsize=9)

    # --- Right: errors (dT_hot / dT_cold) -------------------------------
    plot_with_tail(ax2, lmtd.history_hot, max_len_err, color_lmtd_hoterr, "LMTD - Error dT_hot")
    plot_with_tail(ax2, lmtd.history_cold, max_len_err, COLOR_LMTD, "LMTD - Error dT_cold")
    plot_with_tail(ax2, ntu.history_hot, max_len_err, color_ntu_hoterr, "NTU - Error dT_hot")
    plot_with_tail(ax2, ntu.history_cold, max_len_err, COLOR_NTU, "NTU - Error dT_cold")
    plot_with_tail(ax2, cell.history_hot, max_len_err, color_cell_hoterr, "Cell - Error dT_hot")
    plot_with_tail(ax2, cell.history_cold, max_len_err, COLOR_CELL, "Cell - Error dT_cold")

    ax2.axhline(0, color='black', linewidth=0.8, linestyle=':')
    ax2.set_title(f"Iteration Errors ($\\omega$ = {omega})")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Difference (Current - Previous) [K]")
    ax2.legend(fontsize=9)

    fig.suptitle(f"LMTD vs. NTU vs. Cell Solver Convergence — {label}", fontsize=18, fontweight="bold", y=0.99)

    fig.text(
        0.5, 0.925,
        input_conditions_text,
        ha="center", va="top",
        fontsize=10.5, family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="whitesmoke", edgecolor="gray", alpha=0.9),
    )
    fig.text(
        0.5, 0.83,
        output_text_cell,
        ha="center", va="top",
        fontsize=8.2, family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#fdf0e0", edgecolor=COLOR_CELL, alpha=0.9),
    )
    fig.text(
        0.32, 0.745,
        output_text_lmtd,
        ha="center", va="top",
        fontsize=8.2, family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#eaf5ec", edgecolor=COLOR_LMTD, alpha=0.9),
    )
    fig.text(
        0.68, 0.745,
        output_text_ntu,
        ha="center", va="top",
        fontsize=8.2, family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f3ecf5", edgecolor=COLOR_NTU, alpha=0.9),
    )
    fig.text(
        0.5, 0.665,
        economics_text,
        ha="center", va="top",
        fontsize=8.2, family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#e8eef5", edgecolor="#4a6fa5", alpha=0.9),
    )

    fig.tight_layout(rect=[0, 0, 1, 0.62])
    return fig
