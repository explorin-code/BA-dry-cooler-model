"""
run_scenario.py
================
Shared display pieces for an already-solved ScenarioResult (see
solvers.solve_scenario()): the text formatters used by the cooler-results
figure (plot_profiles.plot_cooler_results) and the terminal summary, and
the convergence panels used by the convergence figure
(plot_solver.plot_convergence). Pure display -- no solving.
"""

from src.solvers import SOLVER_FAMILY, ScenarioResult
from src.inputs.operating_conditions import OperatingConditions
from src.auxiliaries.economics import calc_total_power
from src.plots.plot_style import family_color, lighten_color


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


def calc_pinch(T_c_o: float, T_amb: float) -> float:
    """Pinch point of the WHOLE cooler: coolant outlet minus the AMBIENT air that enters it [K]. With the precooling pad the
    air entering the coil is colder than the ambient; the pad is part of the cooler, so the pinch is not taken against that
    colder air (it can become negative: the pad cools the water below the dry-bulb temperature of the ambient)."""
    return T_c_o - T_amb


def calc_deviation(value: float, reference: float) -> float:
    """Relative deviation of value from reference [%]."""
    return (value - reference) / reference * 100.0


def format_output_conditions(label: str, result, T_a_i: float, reference=None, reference_label: str = "Cell") -> str:
    """Terminal text: outlet temperatures, pinch and Q; for LMTD/NTU
    (reference = the reference SolverResult) their deviation from Cell; then
    final-iteration Pr/Re/Nu for both sides. The figure uses output_condition_cells() instead (subscripts, aligned columns)."""
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
        lines.append(f"vs. {reference_label}:   Q {calc_deviation(result.Q_dot, reference.Q_dot):+5.2f} %   "
                     f"{dT_pinch} {result.T_c_o - reference.T_c_o:+5.2f} K")
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
        ["Geometry:", r"$N_{r}$ =", f"{geo.N_r:d}", r"$N_{t}$ =", f"{geo.N_t:d}",
         r"$l \times b$ =", f"{geo.l:.3f} m × {geo.b:.3f} m",
         r"$A_{fr}$ =", f"{geo.A_fr:.4f} m²"],
    ]
    return [media, geometry]


def economics_cells(P_pump, P_fan, m_dot_ev) -> list:
    """Blocks for the economics box."""
    from src.auxiliaries.economics import calc_total_power
    P_total = calc_total_power(P_pump, P_fan)
    fmt_w = lambda v: f"{v:7.2f} W" if v is not None else "n/a"
    return [[["Economics:", r"$P_{pump}$ =", fmt_w(P_pump), r"$P_{fan}$ =", fmt_w(P_fan),
              r"$P_{tot}$ =", fmt_w(P_total), r"$\dot{m}_{ev}$ =",
              f"{m_dot_ev * 1000:6.3f} g/s" if m_dot_ev is not None else "n/a"]]]


def output_condition_cells(label: str, result, T_a_i: float, reference=None, reference_label: str = "Cell",
                           with_numbers: bool = True) -> list:
    """Blocks for one solver's result box: outlets/pinch/Q, below them (when a reference result is given) the deviation
    from it -- pinch in K, Q in % -- and, with_numbers, an empty line and Pr/Re/Nu/alpha for both sides."""
    d = result.diagnostics
    pinch = calc_pinch(result.T_c_o, T_a_i)
    outputs = [
        ["Outputs:", r"$T_{c,o}$ =", f"{result.T_c_o:5.2f} °C", r"$T_{a,o}$ =",
         f"{result.T_a_o:5.2f} °C", r"$\Delta T_{pinch}$ =", f"{pinch:5.2f} K", r"$Q$ =",
         f"{result.Q_dot / 1000:6.2f} kW"],
    ]
    if reference is not None:
        outputs.append([f"vs. {reference_label}:", "", "", "", "", "", f"{result.T_c_o - reference.T_c_o:+5.2f} K",
                        "", f"{calc_deviation(result.Q_dot, reference.Q_dot):+6.2f} %"])
    blocks = [[f"{label} — Results"], outputs]
    if with_numbers:
        blocks += [[""], [
            ["Coolant:", r"$Pr$ =", f"{d['Pr_c']:6.3f}", r"$Re$ =", f"{d['Re_c']:7.1f}",
             r"$Nu$ =", f"{d['Nu_c']:6.2f}", r"$\alpha_c$ =", f"{d['alpha_c']:7.1f} W/m²K"],
            ["Air:", r"$Pr$ =", f"{d['Pr_a']:6.3f}", r"$Re$ =", f"{d['Re_a']:7.1f}",
             r"$Nu$ =", f"{d['Nu_a']:6.2f}", r"$\alpha_a$ =", f"{d['alpha_a']:7.1f} W/m²K"],
        ]]
    return blocks


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
    """Convergence panels for an already-solved ScenarioResult: outlet temperatures per iteration (ax1) and
    per-iteration errors (ax2), one colour per solver that ran. omega: LMTD/NTU's relaxation factor, cell_omega:
    Cell's own."""
    omega_text = f"$\\omega$: LMTD/NTU = {omega}, Cell = {cell_omega}"
    results = result.results
    max_len_temp = max(len(r.history_T_c_o) for r in results.values())
    max_len_err = max(len(r.history_hot) for r in results.values())

    for key, r in results.items():
        label, color = SOLVER_FAMILY[key], family_color(key)
        plot_with_tail(ax1, r.history_T_c_o, max_len_temp, color, rf"{label} – $T_{{c,o}}$")
        plot_with_tail(ax1, r.history_T_a_o, max_len_temp, lighten_color(color, 0.55), rf"{label} – $T_{{a,o}}$")
        plot_with_tail(ax2, r.history_hot, max_len_err, lighten_color(color, 0.55), rf"{label} – error $\Delta T_{{hot}}$")
        plot_with_tail(ax2, r.history_cold, max_len_err, color, rf"{label} – error $\Delta T_{{cold}}$")

    ax1.axhline(ops.T_a_i, color='blue', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label=r"$T_{a,i}$")
    ax1.axhline(ops.T_c_i, color='red', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label=r"$T_{c,i}$")
    ax1.set_title(f"Outlet Temperatures ({omega_text})")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Temperature [°C]")
    ax1.legend(fontsize=9)

    ax2.axhline(0, color='black', linewidth=0.8, linestyle=':')
    ax2.set_title(f"Iteration Errors ({omega_text})")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Difference (Current - Previous) [K]")
    ax2.legend(fontsize=9)
