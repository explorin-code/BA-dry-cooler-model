"""
benchmark_annual.py
====================
Annual benchmark (ANNUAL_MODE / main.py --annual): runs five solver variants --
LMTD, NTU in 'field' and 'table' mode, Cell in 3D (all tubes) and 2D (one
representative tube; identical results with uniform inlet air) -- over one year of
hourly operating points (8760) at the current geometry, the use case of a
year-round simulation. Only the inlet air (temperature, humidity) changes
from hour to hour; coolant inlet and flows stay at their parameters.py values.

- Weather: ANNUAL_WEATHER_CSV (columns T_air_C, phi_percent, 8760 rows) or,
  if None, a synthetic Munich-like year (reproducible, see synthetic_year) --
  for a performance comparison only the spread of operating points matters,
  not real weather.
- Caching: the property cache is on (the best case for repeated runs) and is
  emptied before each solver's year, so every solver starts equally cold and
  warms up over its own year.
- NTU table: the P(NTU, R) table is rebuilt for this benchmark (never read
  from disk), and its build time is reported separately from the hourly
  calculation -- that split is the point of the comparison.
- Time budget: a solver that exceeds ANNUAL_TIME_BUDGET_S (in practice Cell)
  is stopped and its total extrapolated linearly from the hours it covered.
  All solvers process the hours in the same seeded random order, so a
  stopped run covers the whole year evenly rather than only January. Its
  cache was still warming up, so the extrapolation is slightly pessimistic.

The precooling decision is not part of this benchmark (dry operation only).
"""

import dataclasses
import logging
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

import src.solvers as solvers
from src.benchmark_solvers import _quiet, _clear_property_cache, machine_info
from src.operating_conditions import get_operating_conditions
from src.plot_style import COLOR_LMTD, COLOR_NTU, COLOR_CELL, lighten_color

logger = logging.getLogger(__name__)

N_HOURS = 8760
SOLVER_NAMES = ('LMTD', 'NTU (field)', 'NTU (table)', 'Cell (3D)', 'Cell (2D)')
COLORS = {'LMTD': COLOR_LMTD, 'NTU (field)': COLOR_NTU, 'NTU (table)': lighten_color(COLOR_NTU, 0.45),
          'Cell (3D)': COLOR_CELL, 'Cell (2D)': lighten_color(COLOR_CELL, 0.45)}


# =============================================================================
# Weather
# =============================================================================

def synthetic_year(seed: int = 0):
    """Hourly (T_air [°C], phi [%]) for one synthetic Munich-like year.
    Assumption, not measured data: seasonal cosine (mean 9 °C, amplitude
    9.5 K, coldest around 20 January -- roughly Munich's monthly means),
    daily cycle peaking at 15:00 with a larger amplitude in summer, plus
    multi-day weather swings (AR(1) noise, ~2 days correlation, 3 K).
    Relative humidity is higher in winter and at night (anti-correlated
    with the daily temperature swing), clipped to 25-98 %."""
    rng = np.random.default_rng(seed)
    h = np.arange(N_HOURS)
    day, hour = h / 24.0, h % 24
    season = -np.cos(2 * np.pi * (day - 20) / 365)                      # -1 in January, +1 in July
    daily = np.cos(2 * np.pi * (hour - 15) / 24) * (3.0 + 1.5 * season)

    a = np.exp(-1 / 48)                                                 # AR(1), ~48 h memory
    noise_T, noise_phi = np.zeros(N_HOURS), np.zeros(N_HOURS)
    for i in range(1, N_HOURS):
        noise_T[i] = a * noise_T[i - 1] + 3.0 * np.sqrt(1 - a**2) * rng.standard_normal()
        noise_phi[i] = a * noise_phi[i - 1] + 8.0 * np.sqrt(1 - a**2) * rng.standard_normal()

    T_air = 9.0 + 9.5 * season + daily + noise_T
    phi = np.clip(78.0 - 8.0 * season - 2.5 * daily + noise_phi, 25.0, 98.0)
    return T_air, phi


def load_weather(csv_path):
    """(T_air [°C], phi [%]) from a CSV with header columns T_air_C and phi_percent."""
    data = np.genfromtxt(csv_path, delimiter=',', names=True)
    T_air, phi = np.asarray(data['T_air_C'], float), np.asarray(data['phi_percent'], float)
    if len(T_air) != N_HOURS:
        raise ValueError(f"{csv_path}: expected {N_HOURS} hourly rows, found {len(T_air)}.")
    return T_air, phi


# =============================================================================
# Benchmark
# =============================================================================

def _build_ntu_table(geo, settings):
    """Rebuilds the NTU table (never from disk) and installs it for this run;
    returns the build time [s]."""
    from scipy.interpolate import RectBivariateSpline
    t_start = time.perf_counter()
    P = solvers.build_ntu_table(geo.N_r, settings.ntu_n_elements)
    solvers._ntu_tables[(geo.N_r, settings.ntu_n_elements)] = RectBivariateSpline(
        solvers.NTU_TABLE_NTU, solvers.NTU_TABLE_R, P, kx=3, ky=3)
    return time.perf_counter() - t_start


def _run_year(name, solve, ops_year, order, budget):
    """Solves the hours in `order` until done or the budget is used up.
    Returns per-hour wall times and Q (NaN where not covered)."""
    times, Q = np.full(N_HOURS, np.nan), np.full(N_HOURS, np.nan)
    _clear_property_cache()
    t_run = time.perf_counter()
    with _quiet():
        for h in order:
            t_start = time.perf_counter()
            Q[h] = solve(ops_year[h]).Q_dot
            times[h] = time.perf_counter() - t_start
            if time.perf_counter() - t_run > budget:
                break
    return times, Q


def run_annual_benchmark(geo, settings, modes):
    """results: weather, per-solver hourly times/Q in processing order, NTU
    table build time, totals (measured or extrapolated)."""
    T_air, phi = (load_weather(modes.annual_weather_csv) if modes.annual_weather_csv
                  else synthetic_year())
    source = modes.annual_weather_csv or "synthetic Munich-like year"
    T_c_i = get_operating_conditions(geo).T_c_i
    if T_air.max() >= T_c_i:
        raise ValueError(f"Weather reaches {T_air.max():.1f} °C, not below the coolant inlet {T_c_i:.1f} °C.")

    t_start = time.perf_counter()
    ops_year = [get_operating_conditions(geo, T_a_i=float(T), phi=float(p) / 100.0) for T, p in zip(T_air, phi)]
    t_ops = time.perf_counter() - t_start
    order = np.random.default_rng(1).permutation(N_HOURS)
    logger.info("[Annual] %s: %d operating points (T_a,i %.1f .. %.1f °C) built in %.1f s; budget %.0f s per solver",
                source, N_HOURS, T_air.min(), T_air.max(), t_ops, modes.annual_time_budget)

    t_build = None
    settings_3d = dataclasses.replace(settings, cell_2d=False)
    settings_2d = dataclasses.replace(settings, cell_2d=True)
    runs = {
        'LMTD': lambda o: solvers.solve_it_LMTD(ops=o, geo=geo, settings=settings),
        'NTU (field)': lambda o: solvers.solve_it_NTU(ops=o, geo=geo, settings=settings, mode='field'),
        'NTU (table)': lambda o: solvers.solve_it_NTU(ops=o, geo=geo, settings=settings, mode='table', profile=False),
        'Cell (3D)': lambda o: solvers.solve_it_cell(ops=o, geo=geo, settings=settings_3d),
        'Cell (2D)': lambda o: solvers.solve_it_cell(ops=o, geo=geo, settings=settings_2d),
    }
    results = {'T_air': T_air, 'phi': phi, 'source': source, 'order': order, 'solvers': {}}
    saved_tables = dict(solvers._ntu_tables)
    try:
        for name in SOLVER_NAMES:
            logger.info("[Annual] %s ...", name)
            if name == 'NTU (table)':
                with _quiet():
                    t_build = _build_ntu_table(geo, settings)
            times, Q = _run_year(name, runs[name], ops_year, order, modes.annual_time_budget)
            done = ~np.isnan(times)
            n_done, t_done = int(done.sum()), float(np.nansum(times))
            results['solvers'][name] = {
                'times': times, 'Q': Q, 'n_done': n_done, 'measured': t_done,
                'total': t_done * N_HOURS / n_done,                    # = measured when complete
                'complete': n_done == N_HOURS}
    finally:
        solvers._ntu_tables.clear()
        solvers._ntu_tables.update(saved_tables)
        _clear_property_cache()
    results['t_build'] = t_build

    _log_annual(results)
    return results


def fmt_duration(seconds: float) -> str:
    if seconds < 1:
        return f"{seconds * 1e3:.0f} ms"
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def _log_annual(results):
    s, t_build = results['solvers'], results['t_build']
    logger.info("[Annual benchmark] %s, %s", results['source'], machine_info())
    for name in SOLVER_NAMES:
        r = s[name]
        extra = (f" + table build {fmt_duration(t_build)}" if name == 'NTU (table)' else "")
        status = "" if r['complete'] else f"  (stopped after {r['n_done']} h, extrapolated)"
        logger.info("  %-12s %10s for 8760 h%s  = %.2f ms/h%s", name, fmt_duration(r['total']), extra,
                    r['total'] / N_HOURS * 1e3, status)
    field, table = s['NTU (field)'], s['NTU (table)']
    both = ~np.isnan(field['Q']) & ~np.isnan(table['Q'])
    if both.any():
        logger.info("  NTU table vs. field: Q max deviation %.1e %% over %d h",
                    np.max(np.abs(table['Q'][both] / field['Q'][both] - 1)) * 100, both.sum())
    tip = tipping_point(results)
    if tip is not None:
        logger.info("  NTU table (incl. build) overtakes NTU field after %d operating points", tip)


def tipping_point(results):
    """Operating points after which NTU table incl. its build is cheaper
    than NTU field (both in processing order); None if never."""
    s = results['solvers']
    order = results['order']
    if not (s['NTU (field)']['complete'] and s['NTU (table)']['complete']):
        return None
    field = np.cumsum(np.nan_to_num(s['NTU (field)']['times'][order]))
    table = results['t_build'] + np.cumsum(np.nan_to_num(s['NTU (table)']['times'][order]))
    cheaper = np.nonzero(table < field)[0]
    return int(cheaper[0]) + 1 if len(cheaper) else None


# =============================================================================
# Figure
# =============================================================================

def _cumulative(r, order):
    """Cumulative wall time over the operating points processed (covered part)."""
    t = r['times'][order][:r['n_done']]
    return np.arange(1, r['n_done'] + 1), np.cumsum(t)


def plot_annual_benchmark(results):
    """(a) cumulative wall time over the operating points, linear, zoomed on
    LMTD/NTU so the table build and the tipping point are visible; (b) total
    time for the year per solver, log scale, with NTU table split into build
    and calculation and Cell's extrapolated part hatched; (c) the year's
    inlet air temperature and the resulting heat flow."""
    s, order, t_build = results['solvers'], results['order'], results['t_build']
    fig = plt.figure(figsize=(19, 13))
    fig.canvas.manager.set_window_title("Annual benchmark -- 8760 operating points")
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.4, 1], width_ratios=[1.35, 1])
    ax_cum, ax_bar, ax_year = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

    # --- (a) cumulative time ---------------------------------------------------
    off_scale = []
    y_top = 1.15 * max(s['LMTD']['total'], s['NTU (field)']['total'], t_build + s['NTU (table)']['total'])
    ax_cum.axhspan(0, t_build, color=COLORS['NTU (table)'], alpha=0.25,
                   label=f"NTU table: building P(NTU, R) once ({fmt_duration(t_build)})")
    for name in SOLVER_NAMES:
        r = s[name]
        x, y = _cumulative(r, order)
        offset = t_build if name == 'NTU (table)' else 0.0
        style = dict(color=COLORS[name], linewidth=2.2)
        if name == 'NTU (table)':
            ax_cum.plot([0, 0], [0, t_build], **style)                # the build, before the first point
            style['linestyle'] = '--'
        ax_cum.plot(np.concatenate([[0], x]), offset + np.concatenate([[0], y]), label=name, **style)
        if not r['complete']:
            ax_cum.plot([r['n_done'], N_HOURS], [offset + y[-1], offset + r['total']], color=COLORS[name],
                        linestyle=':', linewidth=2.2)
            off_scale.append(name) if offset + r['total'] > y_top else None
    if any(not s[n]['complete'] for n in SOLVER_NAMES):
        ax_cum.plot([], [], color='0.4', linestyle=':', linewidth=2.2, label="extrapolated (stopped at the budget)")
    for i, name in enumerate(off_scale):
        r = s[name]
        ax_cum.annotate(f"{name} off the scale: ≈ {fmt_duration(r['total'])} for 8760 h\n"
                        f"(measured {r['n_done']} h in {fmt_duration(r['measured'])}, extrapolated)",
                        (0.98, 0.97 - 0.1 * i), xycoords='axes fraction', ha='right', va='top', fontsize=11,
                        color=COLORS[name], fontweight='bold')
    tip = tipping_point(results)
    if tip is not None:
        y_tip = np.cumsum(np.nan_to_num(s['NTU (field)']['times'][order]))[tip - 1]
        ax_cum.plot(tip, y_tip, 'o', color='black', markersize=8, zorder=5)
        ax_cum.annotate(f"tipping point: {tip} operating points", (tip, y_tip), xytext=(15, -25),
                        textcoords='offset points', fontsize=11, arrowprops=dict(arrowstyle='->'))
    ax_cum.set_xlim(0, N_HOURS)
    ax_cum.set_ylim(0, y_top)
    ax_cum.set_xlabel("Operating points (hours) solved")
    ax_cum.set_ylabel("Cumulative wall time [s]")
    ax_cum.set_title("(a) Cumulative wall time over the year", fontsize=13)
    ax_cum.legend(fontsize=10, loc='upper left',
                  bbox_to_anchor=(0.0, 0.97 - 0.1 * len(off_scale)) if off_scale else None)
    ax_cum.grid(alpha=0.3)

    # --- (b) totals, log scale ---------------------------------------------------
    totals = [s[n]['total'] + (t_build if n == 'NTU (table)' else 0.0) for n in SOLVER_NAMES]
    x_min = min(min(totals), t_build) / 5
    extrapolated_labelled = False
    for row, name in enumerate(SOLVER_NAMES):
        r = s[name]
        bar = dict(height=0.6, edgecolor='black', linewidth=0.6)
        if name == 'NTU (table)':
            ax_bar.barh(row, t_build - x_min, left=x_min, color=COLORS[name], hatch='//', **bar,
                        label="table build (once)")
            ax_bar.barh(row, r['total'], left=t_build, color=COLORS[name], **bar)
            text = f" {fmt_duration(t_build + r['total'])}\n (build {fmt_duration(t_build)} + {fmt_duration(r['total'])})"
        elif not r['complete']:
            ax_bar.barh(row, r['measured'] - x_min, left=x_min, color=COLORS[name], **bar)
            ax_bar.barh(row, r['total'] - r['measured'], left=r['measured'], color=lighten_color(COLORS[name], 0.5),
                        hatch='xx', **bar, label=None if extrapolated_labelled else "extrapolated to 8760 h")
            extrapolated_labelled = True
            text = f" ≈ {fmt_duration(r['total'])}\n (measured {r['n_done']} h)"
        else:
            ax_bar.barh(row, r['total'] - x_min, left=x_min, color=COLORS[name], **bar)
            text = f" {fmt_duration(r['total'])}"
        ax_bar.annotate(text, (totals[row], row), va='center', fontsize=11)
    ax_bar.set_xscale('log')
    ax_bar.set_xlim(x_min, max(totals) * 30)
    ax_bar.set_yticks(range(len(SOLVER_NAMES)), SOLVER_NAMES)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Wall time for 8760 operating points [s] (log)")
    ax_bar.set_title("(b) Total for one year", fontsize=13)
    ax_bar.legend(fontsize=10, loc='upper right')
    ax_bar.grid(axis='x', which='both', alpha=0.3)

    # --- (c) the year ------------------------------------------------------------
    hours = np.arange(N_HOURS)
    ax_year.plot(hours / 24, results['T_air'], color='0.55', linewidth=0.5, label=r"$T_\mathrm{a,i}$")
    ax_year.set_xlabel("Day of the year")
    ax_year.set_ylabel(r"$T_\mathrm{a,i}$ [°C]")
    ax_q = ax_year.twinx()
    for name in ('NTU (table)', 'Cell (2D)', 'Cell (3D)'):
        Q = s[name]['Q']
        ok = ~np.isnan(Q)
        if name.startswith('Cell'):
            ax_q.plot(hours[ok] / 24, Q[ok] / 1e3, '.', color=COLORS[name], markersize=4,
                      label=f"$\\dot Q$ {name} ({ok.sum()} h covered)")
        else:
            ax_q.plot(hours[ok] / 24, Q[ok] / 1e3, color=COLORS['NTU (field)'], linewidth=0.5,
                      label=r"$\dot Q$ NTU (table)")
    ax_q.set_ylabel(r"$\dot Q$ [kW]")
    ax_year.set_xlim(0, 365)
    ax_year.set_title(f"(c) Operating points: {results['source']}", fontsize=13, pad=28)
    lines = ax_year.get_legend_handles_labels()
    lines_q = ax_q.get_legend_handles_labels()
    ax_year.legend(lines[0] + lines_q[0], lines[1] + lines_q[1], fontsize=10, loc='lower center',
                   bbox_to_anchor=(0.5, 1.06), ncol=4, frameon=False)

    fig.suptitle("Annual Benchmark — LMTD, NTU (field / table) and Cell (3D / 2D) over 8760 operating points\n"
                 f"(property cache on, emptied per solver; NTU table rebuilt; {machine_info()})",
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig
