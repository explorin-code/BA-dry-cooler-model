"""
main.py
=======
Runs both solvers (LMTD and NTU) back to back with the same
under-relaxation factor and plots their convergence side by side:
outlet temperatures on the left, iteration errors on the right. The
grey box up top echoes the exact input conditions (and geometry) that
were fed into both runs, and the green/purple boxes below it summarize
each solver's converged result.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import seaborn as sns

from src.solvers import solve_it_LMTD, solve_it_NTU
from src.parameters import OperatingConditions
from src.dry_cooler_physics import get_geometry

sns.set_theme(style="whitegrid", context="talk")

# ------------------------------------------------------------------
# Centralized under-relaxation factor.
# Both solvers get the SAME omega so the step size (and therefore the
# iteration count / x-axis) is directly comparable between them.
# ------------------------------------------------------------------
CENTRAL_OMEGA = 0.2

# ------------------------------------------------------------------
# Solver base colors: green (LMTD) / purple (NTU).
# Deliberately avoiding red/blue since that pairing tends to read as
# "cold/hot" -- these colors carry no thermal connotation.
# (Taken from the ColorBrewer PRGn diverging pair -- colorblind-safe
# and visually distinct.)
# ------------------------------------------------------------------
COLOR_LMTD = "#1b7837"  # green
COLOR_NTU = "#762a83"   # purple


def lighten_color(color, factor=0.5):
    """Blend `color` toward white. factor=0 -> unchanged, factor=1 -> white."""
    r, g, b = mc.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


def plot_with_tail(ax, series, max_len, color, label, linewidth=2.2, marker='o', markersize=3):
    """
    Plot `series` as a solid line for the iterations it actually ran.
    If it converged before `max_len`, extend it as a flat dashed line
    (starting from its last value) so both solvers span the same x-range.
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
        f"Air ({ops.air_type}):        "
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


# ------------------------------------------------------------------
# Run both solvers with the same omega
# ------------------------------------------------------------------
(k_lmtd, Q_lmtd, Tc_lmtd, Ta_lmtd,
 hist_hot_lmtd, hist_cold_lmtd,
 hist_Tc_lmtd, hist_Ta_lmtd,
 diag_lmtd) = solve_it_LMTD(omega=CENTRAL_OMEGA)

print(f"[LMTD] converged! k = {k_lmtd:.2f} W/m2K, Q = {Q_lmtd/1000:.2f} kW, "
      f"iterations = {len(hist_hot_lmtd)}")

(k_ntu, Q_ntu, Tc_ntu, Ta_ntu,
 hist_hot_ntu, hist_cold_ntu,
 hist_Tc_ntu, hist_Ta_ntu,
 diag_ntu) = solve_it_NTU(omega=CENTRAL_OMEGA)

print(f"[NTU]  converged! k = {k_ntu:.2f} W/m2K, Q = {Q_ntu/1000:.2f} kW, "
      f"iterations = {len(hist_hot_ntu)}")

# ------------------------------------------------------------------
# Input conditions actually used by the solvers (both solvers construct
# OperatingConditions() with the same defaults, so a fresh instance here
# reflects exactly what was fed into the runs above).
# ------------------------------------------------------------------
ops = OperatingConditions()
geo = get_geometry()
input_conditions_text = format_input_conditions(ops) + "\n" + format_geometry_info(geo)

# ------------------------------------------------------------------
# Per-solver output fields: outlet temps, Q, and final-iteration Pr/Re/Nu
# ------------------------------------------------------------------
output_text_lmtd = format_output_conditions("LMTD", Tc_lmtd, Ta_lmtd, Q_lmtd, diag_lmtd)
output_text_ntu = format_output_conditions("NTU", Tc_ntu, Ta_ntu, Q_ntu, diag_ntu)

# ------------------------------------------------------------------
# Colors: one hue per solver, air = lighter tone of the coolant hue
# ------------------------------------------------------------------
color_lmtd_air = lighten_color(COLOR_LMTD, 0.55)
color_ntu_air = lighten_color(COLOR_NTU, 0.55)

color_lmtd_hoterr = lighten_color(COLOR_LMTD, 0.45)
color_ntu_hoterr = lighten_color(COLOR_NTU, 0.45)

max_len_temp = max(len(hist_Tc_lmtd), len(hist_Tc_ntu))
max_len_err = max(len(hist_hot_lmtd), len(hist_hot_ntu))

# ====================================================================
# Both plots share ONE figure/window, side by side
# ====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 8.5))

# --- Left: outlet temperatures -------------------------------------
plot_with_tail(ax1, hist_Tc_lmtd, max_len_temp, COLOR_LMTD, "LMTD - Coolant out")
plot_with_tail(ax1, hist_Ta_lmtd, max_len_temp, color_lmtd_air, "LMTD - Air out")
plot_with_tail(ax1, hist_Tc_ntu, max_len_temp, COLOR_NTU, "NTU - Coolant out")
plot_with_tail(ax1, hist_Ta_ntu, max_len_temp, color_ntu_air, "NTU - Air out")

ax1.axhline(ops.T_air_in, color='blue', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label="T_air_in")
ax1.axhline(ops.T_coolant_in, color='red', linewidth=1.5, linestyle='--', alpha=0.8, zorder=1, label="T_coolant_in")

ax1.set_title(f"Outlet Temperatures ($\\omega$ = {CENTRAL_OMEGA})")
ax1.set_xlabel("Iteration")
ax1.set_ylabel("Temperature [°C]")
ax1.legend()

# --- Right: errors (dT_hot / dT_cold) -------------------------------
# cold error = full (dark) color, hot error = lighter tone
plot_with_tail(ax2, hist_hot_lmtd, max_len_err, color_lmtd_hoterr, "LMTD - Error dT_hot")
plot_with_tail(ax2, hist_cold_lmtd, max_len_err, COLOR_LMTD, "LMTD - Error dT_cold")
plot_with_tail(ax2, hist_hot_ntu, max_len_err, color_ntu_hoterr, "NTU - Error dT_hot")
plot_with_tail(ax2, hist_cold_ntu, max_len_err, COLOR_NTU, "NTU - Error dT_cold")

ax2.axhline(0, color='black', linewidth=0.8, linestyle=':')
ax2.set_title(f"Iteration Errors ($\\omega$ = {CENTRAL_OMEGA})")
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Difference (Current - Previous) [K]")
ax2.legend()

fig.suptitle("LMTD vs. NTU Solver Convergence", fontsize=18, fontweight="bold", y=0.99)

# --- Input Conditions field: inlet temps + all 3 flow-rate reps for both media ---
fig.text(
    0.5, 0.925,
    input_conditions_text,
    ha="center", va="top",
    fontsize=10.5, family="monospace",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="whitesmoke", edgecolor="gray", alpha=0.9),
)

# --- Output fields: one per solver, placed side-by-side below the input box ---
fig.text(
    0.27, 0.83,
    output_text_lmtd,
    ha="center", va="top",
    fontsize=8.8, family="monospace",
    bbox=dict(boxstyle="round,pad=0.45", facecolor="#eaf5ec", edgecolor=COLOR_LMTD, alpha=0.9),
)
fig.text(
    0.73, 0.83,
    output_text_ntu,
    ha="center", va="top",
    fontsize=8.8, family="monospace",
    bbox=dict(boxstyle="round,pad=0.45", facecolor="#f3ecf5", edgecolor=COLOR_NTU, alpha=0.9),
)

fig.tight_layout(rect=[0, 0, 1, 0.72])

plt.show()