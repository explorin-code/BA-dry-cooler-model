"""
benchmark_solvers.py
=====================
Solver benchmarks (outputs 'caching' / 'resolution'; main.py --outputs caching resolution): the computational cost of
the selected solvers (SOLVERS: LMTD, NTU field / table, Cell 2D / 3D) on one operating point. Runs its own fresh,
silent solves -- never the normal run's, whose property cache is already warm and whose terminal output itself
costs time. Two parts:

1. Caching (run_cache_benchmark): each solver at its configured resolution in three property-cache modes -- no
   cache (every lookup calls CoolProp), cold cache (cache on, starting empty: a normal single run) and warm
   cache (pure algorithm cost) -- with the time split into lookup (CoolProp), lookup (cached), k evaluation and the
   rest of the solver, plus work counts. NTU (table) additionally pays the table build in the no-cache and
   cold-cache runs (a first run needs the table); in the warm-cache runs the table is assumed to be built already.
   Drawn by plot_solver.plot_caching.
2. Resolution (run_resolution_benchmark): NTU (n_elements) and Cell (n_segments) over BENCHMARK_RESOLUTIONS, cache on
   and emptied before every run (= a normal run): wall time, outer iterations and the pinch point; plus the NTU table build time
   up to BENCHMARK_TABLE_MAX_N. NTU runs as the field solver only (the table gives the same results). LMTD has no resolution.
   Drawn by plot_solver.plot_resolution.

The measurement uses the hooks of src/profiling.py (`@timed` on the property lookups, CoolProp calls and the k
evaluation, active only inside `profiling.record()`) and fluid_properties.caching() / clear_cache() for the cache modes -- no
function is patched. The hooks' own overhead (~2-4 % of Cell's time, also outside benchmarks) is included in the totals.
"""


from contextlib import contextmanager
import gc
import logging
import statistics
import time

import dataclasses


import src.thermo.fluid_properties as fluid_properties
import src.solvers as solvers
from src.solvers import ntu_table
from src.solvers.common import mean_states
from src import profiling
from src.formatting import machine_info
from src.quiet import silence_all
from src.plots.plot_style import (TUM_ACCENT_BLUE, TUM_BLUE, TUM_BLUE_540, TUM_GREY_20, TUM_IVORY,
                            solver_color)

logger = logging.getLogger(__name__)

RESULTS_FILE = 'resolution_results.pkl'            # the raw results of the resolution sweep, next to the figure (figures.save_data / load_data)
TABLE_SOLVE_REPEATS = 15                           # a table solve takes ~1 ms: more repeats (right after a long build the first runs are disturbed; N = 30 once showed 2 ms)
CACHE_MODES = ('no cache', 'cold cache', 'warm cache')
TIME_PARTS = ('lookup (CoolProp)', 'lookup (cached)', 'k evaluation', 'rest of solver', 'table build')
PART_COLORS = (TUM_BLUE_540, TUM_BLUE, TUM_ACCENT_BLUE, TUM_GREY_20, TUM_IVORY)      # TUM blues, grey, ivory


# =============================================================================
# Measurement helpers
# =============================================================================

@contextmanager
def _cache_mode(mode: str, fn):
    """Puts the property cache into `mode` before a measured run of fn(): 'no cache' switches the cache off,
    'cold cache' starts it empty, 'warm cache' fills it with one unmeasured run of fn()."""
    if mode == 'no cache':
        with fluid_properties.caching(False):
            yield
        return
    fluid_properties.clear_cache()
    if mode == 'warm cache':
        fn()
    yield


def _measured_run(fn):
    """One counted run of fn(): total wall time, time split and counts (src/profiling.py hooks)."""
    with profiling.record() as rec:
        t_start = time.perf_counter()
        output = fn()
        total = time.perf_counter() - t_start
    t_lookups = rec.time['lookup']
    return {
        'total': total,
        'split': {
            'lookup (CoolProp)': rec.time['coolprop'],
            'lookup (cached)': max(t_lookups - rec.time['coolprop'], 0.0),
            'k evaluation': rec.time['k'],
            'rest of solver': max(total - t_lookups - rec.time['k'], 0.0),
        },
        'iterations': len(output.history_hot),
        'k evaluations': rec.calls['k'],
        'property lookups': rec.calls['lookup'],
        'CoolProp calls': rec.calls['coolprop'],
        'Q [kW]': output.Q_dot / 1000,
    }


@contextmanager
def preserved_ntu_tables():
    """The in-memory NTU tables as they were, restored afterwards (a benchmark installs freshly built ones)."""
    saved = dict(ntu_table.tables)
    try:
        yield
    finally:
        ntu_table.tables.clear()
        ntu_table.tables.update(saved)


def build_fresh_table(N_r: int, n_elements: int) -> float:
    """Builds the NTU table P(NTU, R) from scratch (never from disk), installs it for this run and returns the
    build time [s]. The benchmarks' table covers NTU_a <= 6 (ntu_table.NTU_TABLE_NTU_BENCHMARK, the range of the reported build times); the
    tables the studies use (get_ntu_table) reach NTU_a = 60 and take ~3x longer to build."""
    t_start = time.perf_counter()
    axis = ntu_table.NTU_TABLE_NTU_BENCHMARK
    ntu_table.install_table(N_r, n_elements, ntu_table.build_ntu_table(N_r, n_elements, axis), axis)
    return time.perf_counter() - t_start


def _solver_runs(ops, geo, settings, resolution=None):
    """One zero-argument callable per solver key; NTU elements / Cell segments at `resolution` (None -> the
    settings' own)."""
    if resolution is not None:
        settings = dataclasses.replace(settings, ntu_n_elements=resolution, cell_n_segments=resolution)
    return {key: (lambda key=key: solvers.run_solver(key, ops, geo, settings, profile=False))
            for key in solvers.SOLVER_LABELS}


# =============================================================================
# 1. Caching
# =============================================================================

def run_cache_benchmark(ops, geo, settings, modes):
    """results[solver key][cache mode] = _measured_run dict (median of benchmark_repeats runs by total time; 'no
    cache' runs once -- it is slow and its spread is small relative to its size). NTU (table): the build time is
    added to 'no cache' and 'cold cache' (split part 'table build'), not to 'warm cache'."""
    t_benchmark = time.perf_counter()
    results = {}
    with silence_all(), preserved_ntu_tables():
        for key, fn in _solver_runs(ops, geo, settings).items():
            if key not in modes.solvers:
                continue
            if key == 'ntu_table':
                t_build = build_fresh_table(geo.N_r, settings.ntu_n_elements)       # installed: runs below find it
            results[key] = {}
            for mode in CACHE_MODES:
                runs = []
                for _ in range(1 if mode == 'no cache' else modes.benchmark_repeats):
                    with _cache_mode(mode, fn):
                        runs.append(_measured_run(fn))
                runs.sort(key=lambda r: r['total'])
                run = runs[len(runs) // 2]
                if key == 'ntu_table':
                    build = 0.0 if mode == 'warm cache' else t_build
                    run['total'] += build
                    run['split']['table build'] = build
                results[key][mode] = run
    fluid_properties.clear_cache()

    logger.info("[Benchmark: caching] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    logger.info("  %-12s %-11s %10s %6s %9s %10s %9s %8s",
                "", "", "time", "iter", "k evals", "lookups", "CoolProp", "Q [kW]")
    for key, per_mode in results.items():
        for mode, r in per_mode.items():
            logger.info("  %-12s %-11s %7.1f ms %6d %9d %10d %9d %8.2f", solvers.SOLVER_LABELS[key], mode,
                        r['total'] * 1e3, r['iterations'], r['k evaluations'], r['property lookups'],
                        r['CoolProp calls'], r['Q [kW]'])
    return results


def counts_table_latex(results) -> str:
    """The work counts of the caching benchmark as a LaTeX table (booktabs), one row per solver: outer iterations, k evaluations,
    property lookups and the CoolProp calls without a cache / with a cold / with a warm cache (a cache hit is no call). The
    counts of one solve at the configured resolution; NTU (table) = one lookup solve, its table build is not counted."""
    group = lambda n: f"{n:,}".replace(",", r"\,")
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"Solver & Iterations & $k$ evaluations & Property lookups & \multicolumn{3}{c}{CoolProp calls} \\",
             r"\cmidrule(l){5-7}", r" & & & & no cache & cold cache & warm cache \\", r"\midrule"]
    for key, per_mode in results.items():
        cold = per_mode['cold cache']
        lines.append(" & ".join([solvers.SOLVER_LABELS[key], group(cold['iterations']), group(cold['k evaluations']),
                                 group(cold['property lookups'])]
                                + [group(per_mode[mode]['CoolProp calls']) for mode in CACHE_MODES]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def resolution_labels(settings) -> dict:
    """Resolution each solver runs at (LMTD has none), by solver key."""
    elements, segments = f"{settings.ntu_n_elements} elements", f"{settings.cell_n_segments} segments"
    return {'lmtd': "", 'ntu_field': elements, 'ntu_table': elements, 'cell_2d': segments, 'cell_3d': segments}


# =============================================================================
# 2. Resolution
# =============================================================================

# series of the resolution benchmark: key -> (label, colour, line style)
RESOLUTION_SERIES = {
    'cell_3d': ("Cell (3D)", solver_color('cell_3d'), '-'),
    'cell_2d': ("Cell (2D)", solver_color('cell_2d'), '-'),
    'ntu_field': ("NTU (field)", solver_color('ntu_field'), '-'),      # also stands for NTU (table): same results
}


def _sweep_point(run, ops, n, repeats):
    """One point of the resolution sweep: (n, median wall time [s] of `repeats` runs, iterations, pinch [K]); the property
    cache is emptied before every run (= a normal single run). Returns (row, the solver's result)."""
    times = []
    for _ in range(repeats):
        fluid_properties.clear_cache()
        gc.collect()                                              # garbage of earlier work (a long table build) must not land in the timing
        t_start = time.perf_counter()
        output = run()
        times.append(time.perf_counter() - t_start)
    return (n, statistics.median(times), len(output.history_hot), output.T_c_o - ops.T_a_i), output


def run_resolution_benchmark(ops, geo, settings, modes):
    """results[series key] = list of (n, time [s], iterations, pinch [K]) over modes.benchmark_resolutions for the
    selected NTU / Cell variants (series: RESOLUTION_SERIES), results['lmtd'] = (time [s], iterations, pinch [K]).
    NTU (field) runs once for both NTU variants' RESULTS (the table gives the same results: its P(NTU, R) lookup is the field solved
    once and interpolated); results['ntu_variants'] names which of the two were selected. With the table selected, its tables are
    built fresh for every swept n up to modes.benchmark_table_max_n (the build time grows steeply with n): results['table_build'] =
    [(n, seconds)] and results['ntu_table'] = [(n, time [s] of a solve with the table prebuilt, iterations, pinch [K])];
    results['table_covers'] lists the n whose table covers the operating point (R = C_a/C_c <= 0.1 n) -- for the others the
    solve falls back to the field solve.
    Cache on and emptied before every run (= a normal run); times are the median of modes.benchmark_sweep_repeats runs."""
    t_benchmark = time.perf_counter()
    swept = [k for k in ('ntu_field', 'cell_2d', 'cell_3d')
             if k in modes.solvers or (k == 'ntu_field' and 'ntu_table' in modes.solvers)]
    results = {key: [] for key in swept}
    results['ntu_variants'] = '/'.join(v for v in ('field', 'table') if f'ntu_{v}' in modes.solvers)
    repeats = modes.benchmark_sweep_repeats
    with silence_all():
        for n in modes.benchmark_resolutions:
            runs = _solver_runs(ops, geo, settings, n)
            for key in swept:
                results[key].append(_sweep_point(runs[key], ops, n, repeats)[0])
        if 'lmtd' in modes.solvers:
            fluid_properties.clear_cache()
            t_start = time.perf_counter()
            output = _solver_runs(ops, geo, settings)['lmtd']()
            results['lmtd'] = (time.perf_counter() - t_start, len(output.history_hot), output.T_c_o - ops.T_a_i)
        if 'ntu_table' in modes.solvers:
            results['table_build'], results['ntu_table'], results['table_covers'] = [], [], []
            with preserved_ntu_tables():
                for n in (n for n in modes.benchmark_resolutions if n <= modes.benchmark_table_max_n):
                    results['table_build'].append((n, build_fresh_table(geo.N_r, n)))
                    row, output = _sweep_point(_solver_runs(ops, geo, settings, n)['ntu_table'], ops, n, max(repeats, TABLE_SOLVE_REPEATS))
                    coolant_state, air_state = mean_states(ops, output.T_c_o, output.T_a_o)
                    R = ops.m_dot_a * air_state.cp / (ops.m_dot_c * coolant_state.cp)
                    results['ntu_table'].append(row)
                    if R <= ntu_table.ntu_table_R(n)[-1]:           # else the 'table' solve was a field solve (outside the table)
                        results['table_covers'].append(n)
    fluid_properties.clear_cache()

    logger.info("[Benchmark: resolution] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    for key, (name, _, _) in RESOLUTION_SERIES.items():
        if results.get(key):
            logger.info("  %-28s %s", name, "   ".join(f"n={n}: {t * 1e3:.0f} ms/{it} it" for n, t, it, _ in results[key]))
    if results.get('table_build'):
        logger.info("  %-28s %s", "NTU table build", "   ".join(f"n={n}: {t:.1f} s" for n, t in results['table_build']))
        logger.info("  %-28s %s", "NTU (table, prebuilt)", "   ".join(f"n={n}: {t * 1e3:.2f} ms" for n, t, *_ in results['ntu_table']))
    return results
