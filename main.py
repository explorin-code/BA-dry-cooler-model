import matplotlib.pyplot as plt
import matplotlib.colors as mc
import seaborn as sns

from src.solvers import solve_it_LMTD, solve_it_NTU

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


# ------------------------------------------------------------------
# Run both solvers with the same omega
# ------------------------------------------------------------------
(k_lmtd, Q_lmtd, Tc_lmtd, Ta_lmtd,
 hist_hot_lmtd, hist_cold_lmtd,
 hist_Tc_lmtd, hist_Ta_lmtd) = solve_it_LMTD(omega=CENTRAL_OMEGA)

print(f"[LMTD] converged! k = {k_lmtd:.2f} W/m2K, Q = {Q_lmtd/1000:.2f} kW, "
      f"iterations = {len(hist_hot_lmtd)}")

(k_ntu, Q_ntu, Tc_ntu, Ta_ntu,
 hist_hot_ntu, hist_cold_ntu,
 hist_Tc_ntu, hist_Ta_ntu) = solve_it_NTU(omega=CENTRAL_OMEGA)

print(f"[NTU]  converged! k = {k_ntu:.2f} W/m2K, Q = {Q_ntu/1000:.2f} kW, "
      f"iterations = {len(hist_hot_ntu)}")

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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 6.5))

# --- Left: outlet temperatures -------------------------------------
plot_with_tail(ax1, hist_Tc_lmtd, max_len_temp, COLOR_LMTD, "LMTD - Coolant out")
plot_with_tail(ax1, hist_Ta_lmtd, max_len_temp, color_lmtd_air, "LMTD - Air out")
plot_with_tail(ax1, hist_Tc_ntu, max_len_temp, COLOR_NTU, "NTU - Coolant out")
plot_with_tail(ax1, hist_Ta_ntu, max_len_temp, color_ntu_air, "NTU - Air out")

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

fig.suptitle("LMTD vs. NTU Solver Convergence", fontsize=18, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])

plt.show()
