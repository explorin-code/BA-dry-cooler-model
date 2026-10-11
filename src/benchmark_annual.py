"""
benchmark_annual.py
====================
Annual benchmark (output 'annual' / main.py --outputs annual): how long does a year-round simulation take? The selected
solvers (SOLVERS) run over the 8760 hourly operating points of the Munich climate year at the current geometry;
only the inlet air changes from hour to hour, coolant inlet and flows stay at their parameters.py values.

- Solvers: every selected variant -- LMTD, NTU (field), NTU (table) over every hour, Cell (2D / 3D) on every
  ANNUAL_CELL_STRIDE-th hour only (a quarter by default) as the accuracy / cost comparison; a solver's year total is
  extrapolated when it covers only part of the year. The year panel of the figure draws LMTD, NTU (table) and one Cell
  variant only (NTU field is a fallback if neither LMTD nor NTU table ran).
- NTU table: shown twice -- with the table already built and with the build included (the P(NTU, R) table is
  rebuilt for this benchmark, never read from disk; the build is the price of the first year, nothing afterwards).
- Weather: the Munich climate year of src/weather.py (the same data the cooler comparison uses), or
  ANNUAL_WEATHER_CSV (columns T_air_C, phi_percent[, p_Pa], 8760 rows).
- Caching: the property cache is on (the best case for repeated runs) and is emptied before each solver's
  year, so every solver starts equally cold and warms up over its own year.
- Time budget: a solver that exceeds ANNUAL_TIME_BUDGET_S (in practice Cell) is stopped and its total extrapolated
  linearly from the hours it covered. All solvers process their hours in the same seeded random order.

The precooling decision is not part of this benchmark (dry operation only).
"""


import logging
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

import src.solvers as solvers
from src import weather
from src import fluid_properties
from src.benchmark_solvers import build_fresh_table, preserved_ntu_tables
from src.figures import save_figure
from src.formatting import fmt_duration, machine_info
from src.quiet import silence_all
from src.operating_conditions import get_operating_conditions
from src.plot_style import lighten_color, solver_color

logger = logging.getLogger(__name__)

N_HOURS = 8760
COLORS = {label: solver_color(key) for key, label in solvers.SOLVER_LABELS.items()}
COLORS['NTU (table, incl. build)'] = lighten_color(COLORS['NTU (table)'], 0.5)
COLORS['NTU (table, prebuilt)'] = COLORS['NTU (table)']


# =============================================================================
# Benchmark
# =============================================================================

def _annual_solvers(selected) -> list:
    """The solver keys of the annual run: the selected ones, in the order of solvers.SOLVER_LABELS."""
    return [key for key in solvers.SOLVER_LABELS if key in selected]


def _run_year(solve, ops_year, hours, budget):
    """Solves the hours in the given order until done or the budget is used up.
    Returns per-hour wall times and Q (NaN where not covered)."""
    times, Q = np.full(N_HOURS, np.nan), np.full(N_HOURS, np.nan)
    fluid_properties.clear_cache()
    t_run = time.perf_counter()
    with silence_all():
        for h in hours:
            t_start = time.perf_counter()
            Q[h] = solve(ops_year[h]).Q_dot
            times[h] = time.perf_counter() - t_start
            if time.perf_counter() - t_run > budget:
                break
    return times, Q


def run_annual_benchmark(geo, settings, modes):
    """results: weather, per solver the hourly times / Q, the number of hours measured, the year total (measured or
    extrapolated), the NTU table build time."""
    climate = weather.climate_year(modes.annual_weather_csv)
    T_air, phi = climate['T_air_C'].values, climate['phi_percent'].values
    source = modes.annual_weather_csv or "Munich climate year 2025-10-09 .. 2026-10-08 (DWD 03379)"
    T_c_i = get_operating_conditions(geo).T_c_i
    if T_air.max() >= T_c_i:
        raise ValueError(f"Weather reaches {T_air.max():.1f} °C, not below the coolant inlet {T_c_i:.1f} °C.")

    t_start = time.perf_counter()
    ops_year = [get_operating_conditions(geo, **state) for _, state in weather.hourly_states(climate)]
    t_ops = time.perf_counter() - t_start
    order = np.random.default_rng(1).permutation(N_HOURS)
    stride = max(int(modes.annual_cell_stride), 1)
    logger.info("[Annual] %s: %d operating points (T_a,i %.1f .. %.1f °C) built in %.1f s; budget %.0f s per solver",
                source, N_HOURS, T_air.min(), T_air.max(), t_ops, modes.annual_time_budget)

    results = {'T_air': T_air, 'phi': phi, 'source': source, 'order': order, 'solvers': {}, 't_build': 0.0}
    with preserved_ntu_tables():
        try:
            for key in _annual_solvers(modes.solvers):
                name = solvers.SOLVER_LABELS[key]
                every = stride if key.startswith('cell') else 1
                hours = order[order % every == 0]
                logger.info("[Annual] %s on %d hours ...", name, len(hours))
                if key == 'ntu_table':
                    with silence_all():
                        results['t_build'] = build_fresh_table(geo.N_r, settings.ntu_n_elements)
                solve = lambda o, key=key: solvers.run_solver(key, o, geo, settings, profile=False)
                times, Q = _run_year(solve, ops_year, hours, modes.annual_time_budget)
                n_done, t_done = int((~np.isnan(times)).sum()), float(np.nansum(times))
                results['solvers'][name] = {
                    'times': times, 'Q': Q, 'hours': hours, 'n_done': n_done, 'n_target': len(hours), 'every': every,
                    'measured': t_done, 'total': t_done * N_HOURS / n_done,           # = measured when complete
                    'complete': n_done == N_HOURS}
        finally:
            fluid_properties.clear_cache()

    _log_annual(results)
    return results


def _entries(results):
    """The rows of the figures: name, solver result, time offset [s] (the NTU table appears twice: prebuilt, and with
    its build time added)."""
    entries = []
    for name, r in results['solvers'].items():
        if name == 'NTU (table)':
            entries.append(('NTU (table, prebuilt)', r, 0.0))
            entries.append(('NTU (table, incl. build)', r, results['t_build']))
        else:
            entries.append((name, r, 0.0))
    return entries


def _log_annual(results):
    s, t_build = results['solvers'], results['t_build']
    logger.info("[Annual benchmark] %s, %s", results['source'], machine_info())
    for name, r, offset in _entries(results):
        logger.info("  %-24s %1s%9s for 8760 h  = %.2f ms/h", name, "" if r['complete'] else "≈",
                    fmt_duration(r['total'] + offset), (r['total'] + offset) / N_HOURS * 1e3)
    if 'NTU (table)' in s:
        logger.info("  (NTU table build: %s, once)", fmt_duration(t_build))
    if 'Cell (2D)' in s or 'Cell (3D)' in s:
        cell = s.get('Cell (3D)', s.get('Cell (2D)'))
        for name in ('LMTD', 'NTU (table)'):
            if name in s:
                both = ~np.isnan(s[name]['Q']) & ~np.isnan(cell['Q'])
                if both.any():
                    dev = (s[name]['Q'][both] / cell['Q'][both] - 1) * 100
                    logger.info("  %s vs. Cell: Q mean %+.2f %%, max |dev| %.2f %% over %d h",
                                name, dev.mean(), np.abs(dev).max(), both.sum())


# =============================================================================
# Figure
# =============================================================================

def _axis_scale(totals) -> float:
    """The time axis of panel (a) spans the fastest rows; one that is more than 20x slower than the one before it
    (Cell next to LMTD / NTU) would flatten all the others and is marked 'off the scale' instead."""
    ordered = sorted(totals)
    scale = ordered[0]
    for t in ordered[1:]:
        if t > 20 * scale:
            break
        scale = t
    return scale


def _cumulative(r, order):
    """Cumulative wall time over the operating points processed (covered part)."""
    t = r['times'][r['hours']][:r['n_done']]
    return np.arange(1, r['n_done'] + 1), np.cumsum(t)


def plot_annual_benchmark(results):
    """(a) cumulative wall time over the operating points solved; (b) total time for the year per solver, log
    scale (NTU table with and without its build, Cell's extrapolated part hatched); (c) the year's inlet air
    temperature and the resulting heat flows."""
    s, order = results['solvers'], results['order']
    entries = _entries(results)
    fig = plt.figure(figsize=(19, 13))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.4, 1], width_ratios=[1.35, 1])
    ax_cum, ax_bar, ax_year = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

    # --- (a) cumulative time ---------------------------------------------------
    totals = [r['total'] + offset for _, r, offset in entries]
    y_top = 1.15 * _axis_scale(totals)
    off_scale = []
    for name, r, offset in entries:
        x, y = _cumulative(r, order)
        style = dict(color=COLORS[name], linewidth=2.2)
        if offset:
            ax_cum.plot([0, 0], [0, offset], **style)                  # the build, before the first point
            style['linestyle'] = '--'
        ax_cum.plot(np.concatenate([[0], x]), offset + np.concatenate([[0], y]), label=name, **style)
        if r['n_done'] < N_HOURS:
            ax_cum.plot([r['n_done'], N_HOURS], [offset + y[-1], offset + r['total']], color=COLORS[name],
                        linestyle=':', linewidth=2.2)
            if offset + r['total'] > y_top:
                off_scale.append((name, r))
    if off_scale:
        ax_cum.plot([], [], color='0.4', linestyle=':', linewidth=2.2, label="extrapolated to 8760 h")
    for i, (name, r) in enumerate(off_scale):
        ax_cum.annotate(f"{name} off the scale: ≈ {fmt_duration(r['total'])} for 8760 h",
                        (0.98, 0.97 - 0.1 * i), xycoords='axes fraction', ha='right', va='top', fontsize=11,
                        color=COLORS[name], fontweight='bold')
    ax_cum.set_xlim(0, N_HOURS)
    ax_cum.set_ylim(0, y_top)
    ax_cum.set_xlabel("Operating points (hours) solved")
    ax_cum.set_ylabel("Cumulative wall time [s]")
    ax_cum.set_title("(a) Cumulative wall time over the year", fontsize=13)
    ax_cum.legend(fontsize=10, loc='upper left', bbox_to_anchor=(0.0, 0.97 - 0.1 * len(off_scale)) if off_scale else None)
    ax_cum.grid(alpha=0.3)

    # --- (b) totals, log scale ---------------------------------------------------
    x_min = min(totals) / 5
    extrapolated_labelled = False
    for row, (name, r, offset) in enumerate(entries):
        bar = dict(height=0.6, edgecolor='black', linewidth=0.6)
        total = r['total'] + offset
        if not r['complete']:                                          # measured part + extrapolated rest (hatched)
            ax_bar.barh(row, r['measured'] + offset - x_min, left=x_min, color=COLORS[name], **bar)
            ax_bar.barh(row, r['total'] - r['measured'], left=r['measured'] + offset,
                        color=lighten_color(COLORS[name], 0.5), hatch='xx', **bar,
                        label=None if extrapolated_labelled else "extrapolated to 8760 h")
            extrapolated_labelled = True
            text = f" ≈ {fmt_duration(total)}"
        else:
            ax_bar.barh(row, total - x_min, left=x_min, color=COLORS[name], **bar)
            text = f" {fmt_duration(total)}"
        ax_bar.annotate(text, (total, row), va='center', fontsize=11)
    ax_bar.set_xscale('log')
    ax_bar.set_xlim(x_min, max(totals) * 30)
    ax_bar.set_yticks(range(len(entries)), [e[0] for e in entries])
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Wall time for 8760 operating points [s] (log)")
    ax_bar.set_title("(b) Total for one year", fontsize=13)
    handles, labels = ax_bar.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax_bar.legend(unique.values(), unique.keys(), fontsize=10, loc='upper right')
    ax_bar.grid(axis='x', which='both', alpha=0.3)
    ax_bar.grid(axis='y', visible=False)

    # --- (c) the year ------------------------------------------------------------
    hours = np.arange(N_HOURS)
    ax_year.plot(hours / 24, results['T_air'], color='0.55', linewidth=0.5, label=r"$T_\mathrm{a,i}$")
    ax_year.set_xlabel("Day of the year")
    ax_year.set_ylabel(r"$T_\mathrm{a,i}$ [°C]")
    ax_q = ax_year.twinx()
    ax_q.grid(False)                                                   # one set of grid lines (the temperature axis')
    # the heat flows of LMTD and NTU (table) over the year, plus the samples of ONE Cell variant (2D, else 3D; none
    # if Cell did not run); NTU (field) is not drawn -- only if neither LMTD nor NTU (table) ran, as a fallback
    curves = [n for n in ('LMTD', 'NTU (table)') if n in s] or [n for n in ('NTU (field)',) if n in s]
    for name in curves:
        ok = ~np.isnan(s[name]['Q'])
        ax_q.plot(hours[ok] / 24, s[name]['Q'][ok] / 1e3, color=COLORS[name], linewidth=0.5, label=rf"$\dot Q$ {name}")
    cell = next((n for n in ('Cell (2D)', 'Cell (3D)') if n in s), None)
    if cell:
        ok = ~np.isnan(s[cell]['Q'])
        ax_q.plot(hours[ok] / 24, s[cell]['Q'][ok] / 1e3, '.', color=COLORS[cell], markersize=4,
                  label=f"$\\dot Q$ {cell} ({ok.sum()} h covered)")
    ax_q.set_ylabel(r"$\dot Q$ [kW]")
    ax_year.set_xlim(0, 365)
    ax_year.set_title(f"(c) Operating points: {results['source']}", fontsize=13, loc='left')
    lines = ax_year.get_legend_handles_labels()
    lines_q = ax_q.get_legend_handles_labels()
    ax_year.legend(lines[0] + lines_q[0], lines[1] + lines_q[1], fontsize=10, loc='lower right',
                   bbox_to_anchor=(1.0, 1.0), ncol=4, frameon=False, markerscale=3)

    fig.suptitle(f"Annual performance — {', '.join(s)} over 8760 operating points\n"
                 f"(property cache on, emptied per solver{'; NTU table rebuilt' if 'NTU (table)' in s else ''}; {machine_info()})",
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_figure(fig, 'solver', 'annual_performance', "Annual performance")
    return fig
