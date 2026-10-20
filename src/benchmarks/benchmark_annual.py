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
- Weather: the Munich climate year of src/inputs/weather.py (the same data the cooler comparison uses), or
  ANNUAL_WEATHER_CSV (columns T_air_C, phi_percent[, p_Pa], 8760 rows).
- Caching: the property cache is on (the best case for repeated runs) and is emptied before each solver's
  year, so every solver starts equally cold and warms up over its own year.
- Time budget: a solver that exceeds ANNUAL_TIME_BUDGET_S (in practice Cell) is stopped and its total extrapolated
  linearly from the hours it covered; a budget of 0 means no limit. All solvers process their hours in the same seeded
  random order. THE WHOLE YEAR with every hour solved by every solver (some minutes for Cell 2D, ~2 h for Cell 3D):
      MPLBACKEND=Agg venv/bin/python main.py solver --solvers lmtd ntu_field ntu_table cell_2d cell_3d --outputs annual \
          --annual-cell-stride 1 --annual-time-budget 0
- The raw results are saved before the figures are drawn (annual_results.pkl next to them, `figures.load_data(figures.SOLVER, RESULTS_FILE)`), so a long run is
  never lost to a plotting problem and the figures can be redrawn without it. The figures: plot_solver.plot_annual.

The precooling decision is not part of this benchmark (dry operation only).
"""


import logging
import time

import numpy as np

import src.solvers as solvers
from src.inputs import weather
from src.thermo import fluid_properties
from src.plots import figures
from src.benchmarks.benchmark_solvers import build_fresh_table, preserved_ntu_tables
from src.formatting import fmt_duration, machine_info
from src.quiet import silence_all
from src.inputs.operating_conditions import get_operating_conditions

logger = logging.getLogger(__name__)

N_HOURS = 8760
RESULTS_FILE = 'annual_results.pkl'            # the raw results, next to the figures (figures.save_data / load_data)


# =============================================================================
# Benchmark
# =============================================================================

def _annual_solvers(selected) -> list:
    """The solver keys of the annual run: the selected ones, in the order of solvers.SOLVER_LABELS."""
    return [key for key in solvers.SOLVER_LABELS if key in selected]


def _run_year(solve, ops_year, hours, budget):
    """Solves the hours in the given order until done or the budget is used up.
    Returns per-hour wall times, Q and coolant outlet temperature (NaN where not covered)."""
    times, Q, T_c_o = (np.full(N_HOURS, np.nan) for _ in range(3))
    fluid_properties.clear_cache()
    t_run = time.perf_counter()
    with silence_all():
        for h in hours:
            t_start = time.perf_counter()
            result = solve(ops_year[h])
            times[h] = time.perf_counter() - t_start
            Q[h], T_c_o[h] = result.Q_dot, result.T_c_o
            if budget and time.perf_counter() - t_run > budget:                 # budget 0: no limit
                break
    return times, Q, T_c_o


def run_annual_benchmark(geo, settings, modes):
    """results: weather, per solver the hourly times / Q / T_c,o, the number of hours measured, the year total (measured or
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
                times, Q, T_c_o = _run_year(solve, ops_year, hours, modes.annual_time_budget)
                n_done, t_done = int((~np.isnan(times)).sum()), float(np.nansum(times))
                results['solvers'][name] = {
                    'times': times, 'Q': Q, 'T_c_o': T_c_o, 'hours': hours, 'n_done': n_done, 'n_target': len(hours), 'every': every,
                    'measured': t_done, 'total': t_done * N_HOURS / n_done,           # = measured when complete
                    'complete': n_done == N_HOURS}
        finally:
            fluid_properties.clear_cache()

    results['machine'] = machine_info()
    logger.info("[Annual] raw results saved: %s", figures.save_data(figures.SOLVER, RESULTS_FILE, results))
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
