"""
run_scenario.py
================
Shared display pieces for an already-solved ScenarioResult (see
solvers.solve_scenario()): the text formatters used by the cooler-results
figure (plot_profiles.plot_cooler_results) and the terminal summary, and
the convergence panels used by the convergence figure
(plot_solver.plot_convergence). Pure display -- no solving.
"""

import seaborn as sns

from src.solvers import ScenarioResult
from src.operating_conditions import OperatingConditions
from src.economics import calc_total_power
from src.plot_style import COLOR_LMTD, COLOR_NTU, COLOR_CELL, lighten_color

sns.set_theme(style="whitegrid", context="talk")


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
    """Inlet temps plus all three flow-rate representations (w/ṁ/V̇) for
    both media, regardless of which one was actually specified."""
    coolant_line = (
        f"Coolant ({ops.coolant_type}):  "
        f"T_in = {ops.T_c_i:5.1f} °C   "
        f"w = {ops.w_c:6.3f} m/s   "
        f"ṁ = {ops.m_dot_c:6.3f} kg/s   "
        f"V̇ = {ops.V_dot_c:7.5f} m³/s"
    )
    air_line = (
        f"Air (phi = {ops.phi * 100:3.0f} %):     "
        f"T_in = {ops.T_a_i:5.1f} °C   "
        f"w = {ops.w_fr:6.3f} m/s   "
        f"ṁ = {ops.m_dot_a:6.3f} kg/s   "
        f"V̇ = {ops.V_dot_a:7.5f} m³/s"
    )
    return coolant_line + "\n" + air_line


def format_geometry_info(geo) -> str:
    """Tube rows/count and frontal (inflow) area as height x width = area."""
    return (
        f"Geometry:      "
        f"n_rows = {geo.N_r:3d}   "
        f"n_tubes = {geo.N_t:3d}   "
        f"{geo.l:.3f} m x {geo.b:.3f} m = {geo.A_fr:.4f} m²"
    )


def calc_pinch(T_c_o: float, T_a_i: float) -> float:
    """Cold-end temperature approach (pinch point) of the counterflow
    arrangement: coolant outlet vs. air inlet [K]."""
    return T_c_o - T_a_i


def calc_deviation(value: float, reference: float) -> float:
    """Relative deviation of value from reference [%]."""
    return (value - reference) / reference * 100.0


def format_output_conditions(label: str, result, T_a_i: float, reference=None) -> str:
    """Terminal text: outlet temperatures, pinch and Q; for LMTD/NTU
    (reference = the Cell SolverResult) their deviation from Cell; then
    final-iteration Pr/Re/Nu for both sides. The figure uses
    output_condition_cells() instead (subscripts, aligned columns)."""
    diagnostics = result.diagnostics
    pinch = calc_pinch(result.T_c_o, T_a_i)
    T_co, T_ao, dT_pinch = "T_c,o", "T_a,o", "ΔT_pinch"
    alpha_c, alpha_a = "α_c", "α_a"

    lines = [
        f"{label} — Results",
        f"{T_co} = {result.T_c_o:5.2f} °C   {T_ao} = {result.T_a_o:5.2f} °C   "
        f"{dT_pinch} = {pinch:5.2f} K   Q = {result.Q_dot/1000:6.2f} kW",
    ]
    if reference is not None:
        pinch_ref = calc_pinch(reference.T_c_o, T_a_i)
        lines.append(f"vs. Cell:   Q {calc_deviation(result.Q_dot, reference.Q_dot):+5.2f} %   "
                     f"{dT_pinch} {calc_deviation(pinch, pinch_ref):+5.2f} %")
    lines += [
        "",
        f"Coolant:  Pr = {diagnostics['Pr_c']:6.3f}   "
        f"Re = {diagnostics['Re_c']:8.1f}   "
        f"Nu = {diagnostics['Nu_c']:7.2f}   "
        f"{alpha_c} = {diagnostics['alpha_c']:7.1f} W/m²K",
        f"Air:      Pr = {diagnostics['Pr_a']:6.3f}   "
        f"Re = {diagnostics['Re_a']:8.1f}   "
        f"Nu = {diagnostics['Nu_a']:7.2f}   "
        f"{alpha_a} = {diagnostics['alpha_a']:7.1f} W/m²K",
    ]
    return "\n".join(lines)


# -----------------------------------------------------------------------------
# Figure cells for plot_style.draw_text_grid: subscripts as mathtext, one
# value per cell so columns align. Values keep fixed-width formats, so equal
# quantities line up digit by digit within a column.
# -----------------------------------------------------------------------------

def input_condition_cells(ops: OperatingConditions, geo) -> list:
    """Blocks for the input box: both media (same columns), then geometry."""
    media = [
        [f"Coolant ({ops.coolant_type}):", r"$T_{c,i}$ =", f"{ops.T_c_i:5.1f} °C",
         r"$w$ =", f"{ops.w_c:6.3f} m/s", r"$\dot{m}$ =", f"{ops.m_dot_c:6.3f} kg/s",
         r"$\dot{V}$ =", f"{ops.V_dot_c:8.5f} m³/s"],
        [f"Air (φ = {ops.phi * 100:.0f} %):", r"$T_{a,i}$ =", f"{ops.T_a_i:5.1f} °C",
         r"$w$ =", f"{ops.w_fr:6.3f} m/s", r"$\dot{m}$ =", f"{ops.m_dot_a:6.3f} kg/s",
         r"$\dot{V}$ =", f"{ops.V_dot_a:8.5f} m³/s"],
    ]
    geometry = [
        ["Geometry:", r"$n_{rows}$ =", f"{geo.N_r:d}", r"$n_{tubes}$ =", f"{geo.N_t:d}",
         r"$H \times W$ =", f"{geo.l:.3f} m × {geo.b:.3f} m",
         r"$A_{front}$ =", f"{geo.A_fr:.4f} m²"],
    ]
    return [media, geometry]


def economics_cells(P_pump, P_fan, m_dot_ev) -> list:
    """Blocks for the economics box."""
    from src.economics import calc_total_power
    P_total = calc_total_power(P_pump, P_fan)
    fmt_w = lambda v: f"{v:7.2f} W" if v is not None else "n/a"
    return [[["Economics:", r"$P_{pump}$ =", fmt_w(P_pump), r"$P_{fan}$ =", fmt_w(P_fan),
              r"$P_{total}$ =", fmt_w(P_total), r"$\dot{m}_{water}$ =",
              f"{m_dot_ev * 1000:6.3f} g/s" if m_dot_ev is not None else "n/a"]]]


def output_condition_cells(label: str, result, T_a_i: float, reference=None) -> list:
    """Blocks for one solver's result box: outlets/pinch/Q, for LMTD/NTU the
    deviation from Cell directly below the pinch and Q values, an empty
    line, then Pr/Re/Nu/alpha for both sides (own aligned columns)."""
    d = result.diagnostics
    pinch = calc_pinch(result.T_c_o, T_a_i)
    outputs = [
        ["Outputs:", r"$T_{c,o}$ =", f"{result.T_c_o:5.2f} °C", r"$T_{a,o}$ =",
         f"{result.T_a_o:5.2f} °C", r"$\Delta T_{pinch}$ =", f"{pinch:5.2f} K", r"$Q$ =",
         f"{result.Q_dot / 1000:6.2f} kW"],
    ]
    if reference is not None:
        pinch_ref = calc_pinch(reference.T_c_o, T_a_i)
        outputs.append(["vs. Cell:", "", "", "", "", "", f"{calc_deviation(pinch, pinch_ref):+5.2f} %",
                        "", f"{calc_deviation(result.Q_dot, reference.Q_dot):+6.2f} %"])
    numbers = [
        ["Coolant:", r"$Pr$ =", f"{d['Pr_c']:6.3f}", r"$Re$ =", f"{d['Re_c']:7.1f}",
         r"$Nu$ =", f"{d['Nu_c']:6.2f}", r"$\alpha_c$ =", f"{d['alpha_c']:7.1f} W/m²K"],
        ["Air:", r"$Pr$ =", f"{d['Pr_a']:6.3f}", r"$Re$ =", f"{d['Re_a']:7.1f}",
         r"$Nu$ =", f"{d['Nu_a']:6.2f}", r"$\alpha_a$ =", f"{d['alpha_a']:7.1f} W/m²K"],
    ]
    return [[f"{label} — Results"], outputs, [""], numbers]


def format_economics(P_pump, P_fan, m_dot_ev) -> str:
    """Build the 'Economics' box: pump/fan/total power and water usage.
    Any value left as None (not yet available) prints as 'n/a'."""
    pump_str = f"{P_pump:7.2f} W" if P_pump is not None else "    n/a"
    fan_str = f"{P_fan:7.2f} W" if P_fan is not None else "    n/a"
    P_total = calc_total_power(P_pump, P_fan)
    total_str = f"{P_total:7.2f} W" if P_total is not None else "    n/a"
    water_str = f"{m_dot_ev * 1000:6.3f} g/s" if m_dot_ev is not None else "   n/a"

    return (
        f"Economics — Pump: {pump_str}   Fan: {fan_str}   "
        f"Total: {total_str}   Water: {water_str}"
    )


def draw_convergence(ax1, ax2, result: ScenarioResult, ops: OperatingConditions, omega: float,
                     cell_omega: float):
    """Convergence panels for an already-solved ScenarioResult: outlet
    temperatures per iteration (ax1) and per-iteration errors (ax2).
    omega: LMTD/NTU's relaxation factor, cell_omega: Cell's own."""
    omega_text = f"$\\omega$: LMTD/NTU = {omega}, Cell = {cell_omega}"
    lmtd, ntu, cell = result.lmtd, result.ntu, result.cell

    # --- Colors: one hue per solver, air = lighter tone of the coolant hue ---
    color_lmtd_air = lighten_color(COLOR_LMTD, 0.55)
    color_ntu_air = lighten_color(COLOR_NTU, 0.55)
    color_cell_air = lighten_color(COLOR_CELL, 0.55)

    color_lmtd_hoterr = lighten_color(COLOR_LMTD, 0.45)
    color_ntu_hoterr = lighten_color(COLOR_NTU, 0.45)
    color_cell_hoterr = lighten_color(COLOR_CELL, 0.45)

    max_len_temp = max(len(lmtd.history_T_c_o), len(ntu.history_T_c_o), len(cell.history_T_c_o))
    max_len_err = max(len(lmtd.history_hot), len(ntu.history_hot), len(cell.history_hot))

    # --- Left: outlet temperatures -------------------------------------
    plot_with_tail(ax1, lmtd.history_T_c_o, max_len_temp, COLOR_LMTD, r"LMTD – $T_{c,o}$")
    plot_with_tail(ax1, lmtd.history_T_a_o, max_len_temp, color_lmtd_air, r"LMTD – $T_{a,o}$")
    plot_with_tail(ax1, ntu.history_T_c_o, max_len_temp, COLOR_NTU, r"NTU – $T_{c,o}$")
    plot_with_tail(ax1, ntu.history_T_a_o, max_len_temp, color_ntu_air, r"NTU – $T_{a,o}$")
    plot_with_tail(ax1, cell.history_T_c_o, max_len_temp, COLOR_CELL, r"Cell – $T_{c,o}$")
    plot_with_tail(ax1, cell.history_T_a_o, max_len_temp, color_cell_air, r"Cell – $T_{a,o}$")

    ax1.axhline(ops.T_a_i, color='blue', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label=r"$T_{a,i}$")
    ax1.axhline(ops.T_c_i, color='red', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label=r"$T_{c,i}$")

    ax1.set_title(f"Outlet Temperatures ({omega_text})")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Temperature [°C]")
    ax1.legend(fontsize=9)

    # --- Right: errors (Delta_T_hot / Delta_T_cold) -------------------------------
    plot_with_tail(ax2, lmtd.history_hot, max_len_err, color_lmtd_hoterr, r"LMTD – error $\Delta T_{hot}$")
    plot_with_tail(ax2, lmtd.history_cold, max_len_err, COLOR_LMTD, r"LMTD – error $\Delta T_{cold}$")
    plot_with_tail(ax2, ntu.history_hot, max_len_err, color_ntu_hoterr, r"NTU – error $\Delta T_{hot}$")
    plot_with_tail(ax2, ntu.history_cold, max_len_err, COLOR_NTU, r"NTU – error $\Delta T_{cold}$")
    plot_with_tail(ax2, cell.history_hot, max_len_err, color_cell_hoterr, r"Cell – error $\Delta T_{hot}$")
    plot_with_tail(ax2, cell.history_cold, max_len_err, COLOR_CELL, r"Cell – error $\Delta T_{cold}$")

    ax2.axhline(0, color='black', linewidth=0.8, linestyle=':')
    ax2.set_title(f"Iteration Errors ({omega_text})")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Difference (Current - Previous) [K]")
    ax2.legend(fontsize=9)
