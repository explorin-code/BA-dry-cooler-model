"""Smoke tests: every workflow and benchmark runs end to end on a reduced problem, and the numbers that must not
move without a decision (the fixed tube counts, the independence of the geometry from the operating inputs, the
air-side power) stay put. All of it runs in well under a minute (`venv/bin/python -m pytest`)."""

import numpy as np
import pytest

import src.parameters as parameters
from src.cooler_sets import CELL_SETS
from src.dry_cooler_physics import get_geometry
from src.run_modes import get_run_modes
from src.scenario_pipeline import run_scenarios


FAST_SOLVERS = ['lmtd', 'ntu_field']          # Cell 3D costs seconds per scenario: only the tests that need it ask for it


def modes(**kwargs):
    base = dict(insight_mode=False, outputs=[], solvers=FAST_SOLVERS)
    return get_run_modes(**{**base, **kwargs})


# --- geometry -----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("key, tubes", [('P1', 34), ('P2', 32), ('P3', 30)])
def test_default_tube_count_is_the_reference(key, tubes):
    parameters.GEOMETRY_SET = key
    assert get_geometry().N_t == CELL_SETS[key].N_t_ref == tubes


def test_geometry_does_not_follow_operating_inputs():
    reference = get_geometry().N_t
    for name, value in (('T_COOLANT_IN', 80.0), ('T_COOLANT_TARGET_OUT', 30.0), ('M_COOLANT', 0.28),
                        ('T_AIR_IN', 35.0), ('W_O', 1.5), ('PHI_AIR_PERCENT', 20.0)):
        setattr(parameters, name, value)
        assert get_geometry().N_t == reference, name


def test_reference_tube_counts_match_the_design_sizing():
    """N_t_ref is a copy of comparison.sizing.design_tubes -- refresh it when this fails (cells or solver changed)."""
    from comparison import config, sizing
    cfg = config.settings()
    for key in ('P1', 'P2', 'P3'):
        assert sizing.design_tubes(cfg, key) == CELL_SETS[key].N_t_ref, key


# --- solver comparison -------------------------------------------------------------------------------------

def test_solver_workflow_dry_and_precooled():
    parameters.T_AIR_IN = 27.0                   # P2 at 2.5 m/s misses 25 degC -> precooling engaged
    results = run_scenarios(modes(solvers=['lmtd', 'ntu_field', 'cell_3d']))        # the pinned fan numbers use Cell's outlet air
    assert results['precooled'] is not None
    dry, wet = results['ambient'], results['precooled']
    assert set(dry['result'].results) == {'lmtd', 'ntu_field', 'cell_3d'}
    assert dry['m_dot_ev'] is None and wet['m_dot_ev'] > 0
    assert wet['result'].reference.T_c_o < dry['result'].reference.T_c_o
    # the fan power moves with the air state only -- these two numbers are the fan basis of the model (P2, 27 degC)
    assert dry['P_fan'] == pytest.approx(386.66, abs=0.05)
    assert wet['P_fan'] == pytest.approx(734.58, abs=0.05)
    assert wet['m_dot_ev'] * 1000 == pytest.approx(10.8, abs=0.05)


def test_water_and_fan_do_not_depend_on_the_coolant_temperature():
    """The pad only sees the air; the geometry is fixed (the bug of 2026-10: it used to be re-sized)."""
    parameters.T_AIR_IN = 27.0
    cold = run_scenarios(modes())['precooled']
    parameters.T_COOLANT_IN = 80.0
    hot = run_scenarios(modes())['precooled']
    assert hot['m_dot_ev'] == pytest.approx(cold['m_dot_ev'], rel=1e-9)


def test_figures_are_built():
    parameters.T_AIR_IN = 27.0
    run_scenarios(modes(outputs=['cooler_results', 'convergence']))
    import matplotlib.pyplot as plt
    titles = {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}
    assert {'Cooler results - Dry', 'Cooler results - Precooled', 'Convergence - Dry', 'Convergence - Precooled'} <= titles


def test_outputs_follow_the_selection_in_catalogue_order():
    results = run_scenarios(modes(outputs=['resolution', 'cooler_results'], benchmark_resolutions=[2, 3]))
    assert list(results['outputs']) == ['cooler_results', 'resolution']          # catalogue order, only the selected


def test_unknown_output_is_rejected():
    with pytest.raises(ValueError, match="Unknown item"):
        run_scenarios(modes(outputs=['nonsense']))


# --- solver benchmarks ---------------------------------------------------------------------------------------

def test_caching_and_resolution_benchmarks():
    results = run_scenarios(modes(outputs=['caching', 'resolution'], benchmark_repeats=1, benchmark_resolutions=[2, 3],
                                  solvers=['lmtd', 'ntu_field', 'cell_2d']))['outputs']
    assert set(results['caching']) == {'lmtd', 'ntu_field', 'cell_2d'}
    assert all(set(per_mode) == {'no cache', 'cold cache', 'warm cache'} for per_mode in results['caching'].values())
    from src.benchmark_solvers import counts_table_latex
    table = counts_table_latex(results['caching'])
    assert table.count("\\\\\n") == 5 and 'Cell (2D)' in table and r'\begin{tabular}' in table
    resolution = results['resolution']
    assert 'table_build' not in resolution                       # NTU (table) not selected
    assert [row[0] for row in resolution['ntu_field']] == [2, 3] and [row[0] for row in resolution['cell_2d']] == [2, 3]
    assert 'lmtd' in resolution
    import matplotlib.pyplot as plt
    titles = {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}
    assert {'Caching times', 'Caching time split', 'Caching speed-up', 'Caching presentation', 'Resolution performance'} <= titles


def test_annual_benchmark_runs_and_extrapolates():
    results = run_scenarios(modes(outputs=['annual'], solvers=['lmtd'], annual_time_budget=1.0))['outputs']['annual']
    lmtd = results['solvers']['LMTD']
    assert lmtd['n_done'] > 0 and lmtd['total'] > 0
    assert lmtd['complete'] == (lmtd['n_done'] == 8760)
    assert not np.isnan(lmtd['T_c_o'][lmtd['hours'][0]])


def test_annual_benchmark_without_budget_solves_everything_and_draws_the_deviation(tmp_path):
    from src.benchmark_annual import RESULTS_FILE
    from src.figures import SOLVER, load_data
    results = run_scenarios(modes(outputs=['annual'], solvers=['lmtd', 'cell_2d'], annual_time_budget=0.0,
                                  annual_cell_stride=730))['outputs']['annual']
    assert results['solvers']['LMTD']['complete'] and results['solvers']['LMTD']['n_done'] == 8760
    cell = results['solvers']['Cell (2D)']
    assert cell['n_done'] == cell['n_target'] == 12 and not cell['complete']             # every 730th hour, extrapolated
    assert load_data(SOLVER, RESULTS_FILE)['solvers']['Cell (2D)']['n_done'] == 12                        # the raw results were saved
    import matplotlib.pyplot as plt
    titles = {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}
    assert {'Annual times', 'Annual deviation'} <= titles


def test_accuracy_sweep_compares_to_the_reference():
    from src.benchmark_accuracy import deviations
    sweep = run_scenarios(modes(outputs=['accuracy'], solvers=['lmtd', 'ntu_field', 'cell_2d'],
                                accuracy_t_coolant=[37.0, 80.0]))['outputs']['accuracy']
    assert sweep['reference'] == 'cell_2d' and sweep['x'] == [37.0, 80.0] and sweep['x_run'] == 37.0
    assert deviations(sweep, 'cell_2d') == ([0.0, 0.0], [0.0, 0.0])
    dQ, dpinch = deviations(sweep, 'lmtd')
    assert all(abs(d) < 5.0 for d in dQ) and dQ != [0.0, 0.0]
    import matplotlib.pyplot as plt
    assert 'Accuracy sweep' in {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}


def test_profile_grids_solve_a_second_water_temperature():
    cases = run_scenarios(modes(outputs=['profile_grids'], solvers=['lmtd', 'ntu_field'],
                                profile_t_coolant=90.0))['outputs']['profile_grids']
    assert [ops.T_c_i for ops, _ in cases] == [37.0, 90.0]
    assert cases[1][1].reference.Q_dot > 3 * cases[0][1].reference.Q_dot
    import matplotlib.pyplot as plt
    titles = {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}
    assert {'Profile grid (absolute)', 'Profile grid (normalised)', 'Profile graphs (normalised, solvers combined)',
            'Profile results table'} <= titles


def test_accuracy_sweep_needs_two_solvers():
    with pytest.raises(ValueError, match="at least two"):
        run_scenarios(modes(outputs=['accuracy'], solvers=['lmtd']))


# --- cooler comparison ---------------------------------------------------------------------------------------

FAST = ['geometry', 'ambient_cases', 'airflow_sweep', 'sensitivity', 'rows', 'tubes', 'wall']


@pytest.mark.parametrize("analysis", FAST)
def test_comparison_analysis(analysis, outputs, monkeypatch):
    from comparison import sensitivity
    monkeypatch.setattr(sensitivity, 'N_POINTS', 3)                      # the shape of the study, not its resolution
    from comparison.runner import run_analyses
    results = run_analyses([analysis])
    assert results[analysis] and all(len(df) for df in results[analysis].values())
    assert list((outputs / "figures" / "ch4_results").glob("*/*.png"))   # at least one figure was saved


def test_comparison_year_with_fan_control(outputs):
    parameters.COMPARE_YEAR_STRIDE = 48
    from comparison.runner import run_analyses
    data = run_analyses(['fan_control'])['fan_control']
    control = data['control']
    assert set(control['set']) == {'P1', 'P2', 'P3'}
    assert (control['w_fr'] <= parameters.COMPARE_W_FR_DESIGN + 1e-9).all()
    assert (control['w_fr'] >= parameters.COMPARE_W_FR_MIN - 1e-9).all()
    assert data['fixed']['w_fr'].eq(parameters.COMPARE_W_FR_DESIGN).all()
    # the controlled hours that are not limited land on the target within TARGET_TOL, and every hour is either met or limited
    assert control.loc[~control['limited'], 'target_met'].all()
    assert not control.loc[control['limited'], 'target_met'].any()
    from comparison import config, operation
    summary = operation.annual_summary(config.settings(), control, hours_per_row=48)
    assert (summary.h_target_met + summary.h_limited == len(control) / 3 * 48).all()


def test_pinch_is_that_of_the_whole_cooler():
    """The pinch is taken against the AMBIENT air also when precooled; the coil's own approach is pinch_coil."""
    from comparison import config, evaluate as ev, sizing
    cfg = config.settings()
    geo = sizing.get_sized(cfg, 'P2')
    ops = ev.make_ops(cfg, geo, 30.0, 0.3, 95500.0)
    dry, wet = ev.run_point(cfg, geo, ops, 'ntu_field', 'dry'), ev.run_point(cfg, geo, ops, 'ntu_field', 'precooled')
    assert dry['pinch'] == dry['pinch_coil'] == pytest.approx(dry['T_c_o'] - 30.0)
    assert wet['pinch'] == pytest.approx(wet['T_c_o'] - 30.0)
    assert wet['pinch_coil'] == pytest.approx(wet['T_c_o'] - wet['T_a_i']) and wet['T_a_i'] < 30.0 and wet['pinch'] < wet['pinch_coil']


def test_geometry_has_the_print_size_drawing(outputs):
    from comparison.runner import run_analyses
    run_analyses(['geometry'])
    assert (outputs / "figures" / "ch4_results" / "geometry" / "geometry_print.png").exists()


def test_comparison_pinch_fan_figures_and_data(outputs):
    parameters.COMPARE_YEAR_STRIDE = 48
    from comparison import operation
    from comparison.runner import run_analyses
    data = run_analyses(['pinch_fan'])['pinch_fan']
    year = data['year']
    assert not data['control_dry'][data['control_dry']['mode'] == 'precooled'].shape[0]                 # no pad in the dry adaptive year
    assert data['control_dry'].limited.sum() >= data['control'].limited.sum()                          # dry misses at least the hours the pad misses
    assert set(year['mode']) == {'dry', 'precooled'} and len(year) == 2 * len(year[year['mode'] == 'dry'])    # both modes, every hour
    wet = year[year['mode'] == 'precooled']
    assert np.allclose(wet.pinch, wet.T_c_o - wet.T_a_amb)                                            # pinch of the whole cooler
    assert (wet.pinch_coil >= wet.pinch - 1e-9).all()                                                 # the pad only cools the air
    assert 15 < operation.target_crossing(year[year.set == 'P2'], 'dry', 25.0) < 25                   # the dry outlet passes 25 degC near 20 degC ambient
    import matplotlib.pyplot as plt
    assert {'Pinch year', 'Pinch speed', 'Controlled fan'} <= {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()} or \
        {p.stem for p in (outputs / "figures" / "ch4_results" / "cooler_performance").glob("*.png")} >= {'pinch_year', 'pinch_year_clouds', 'missed_hours', 'pinch_speed', 'controlled_fan', 'fan_cop', 'annual_costs', 'monthly_cost', 'water_vs_hours'}


def test_variation_table_and_figures(outputs):
    """The geometry variations: a reduced range set through the real code path (worker processes, tables, plot); the cooler's own value is in every curve,
    the dotted / open-marker classes of the figure exist, and the costs add up."""
    from comparison import config, variation
    from comparison.plots_print import variation_figures
    parameters.COMPARE_YEAR_STRIDE = 96
    cfg = config.settings()
    ranges = {'N_r': (6, 6, 1), 'N_t': (20, 40, 2), 'l': (1.0, 3.0, 2), 's_l_mm': (20.0, 30.0, 2), 's_f_mm': (2.0, 4.0, 2), 'delta_f_mm': (0.1, 0.2, 2)}
    df = variation.run(cfg, ranges=ranges)
    assert set(df.set) == set(cfg.sets) and set(df.param) == set(ranges)
    assert (df[df.standard].groupby(['set', 'param']).size() == 1).all()                              # one standard point per curve
    std = df[df.standard & (df.param == 'N_r')].set_index('set')
    assert (std.pinch > 0).all() and (std.E_fan_kWh > 0).all() and (std.m_Cu > 0).all()
    costs = variation.add_costs(cfg, df)
    assert np.allclose(costs.cost_total, costs.cost_fan + costs.cost_pump + costs.cost_water + costs.cost_material)
    from comparison import operation
    assert np.isclose(costs[costs.standard & (costs.set == 'P2')].cost_material.iloc[0], operation.material_writeoff(cfg, 'P2'))     # same write-off as the annual summary
    variation_figures(df, cfg)
    assert {'variation_costs.png', 'variation_totals.png'} <= {p.name for p in (outputs / "figures" / "ch4_results" / "cooler_performance").glob("*.png")}


def test_fixed_year_rule_matches_the_solved_fixed_year():
    """operation.fixed_year picks dry / precooled per hour from the both-modes year exactly as year_fixed solves it."""
    from comparison import config, operation
    cfg = config.settings()
    parameters.COMPARE_YEAR_STRIDE = 96
    both = operation.year_modes(cfg, stride=96)
    derived = operation.fixed_year(both)
    solved = operation.year_fixed(cfg, stride=96)
    key = ['set', 'hour']
    merged = derived.merge(solved, on=key, suffixes=('_derived', '_solved'))
    assert len(merged) == len(solved) == len(derived)
    assert np.allclose(merged.P_fan_derived, merged.P_fan_solved) and (merged.mode_derived == merged.mode_solved).all()


def test_comparison_sizing_is_the_reference_count():
    from comparison import config, sizing
    cfg = config.settings()
    assert [sizing.get_sized(cfg, k).N_t for k in cfg.sets] == [34, 32, 30]


def test_comparison_settings_are_read_when_requested_not_at_import():
    from comparison import config
    before = config.settings()
    parameters.COMPARE_W_FR_DESIGN = 2.0
    parameters.NTU_OMEGA = 0.5
    parameters.COMPARE_SETS = ['P3']
    after = config.settings()
    assert before.w_fr_design == 3.0 and after.w_fr_design == 2.0          # a Settings object is a snapshot
    assert after.solver_settings.ntu_omega == 0.5 and after.solver_settings.cell_2d    # the solver tuning follows parameters.py
    assert list(after.sets) == ['P3']


# --- instrumentation -----------------------------------------------------------------------------------------

def test_profiling_counts_and_cache_switch_restore_everything():
    from src import fluid_properties, profiling, solvers
    from src.dry_cooler_physics import get_geometry
    from src.operating_conditions import get_operating_conditions
    from src.solver_settings import get_solver_settings
    geo = get_geometry()
    ops, settings = get_operating_conditions(geo=geo), get_solver_settings()
    lookups_before = (fluid_properties._fluid_lookup, fluid_properties._air_lookup)

    with profiling.record() as outer:
        with fluid_properties.caching(False):
            solvers.run_solver('lmtd', ops, geo, settings, profile=False)
    assert outer.calls['k'] == 19 and outer.calls['lookup'] == 38           # LMTD: 19 iterations, 2 lookups each
    assert outer.calls['coolprop'] == 5 * 19 + 5 * 19                      # no cache: every lookup calls CoolProp (5 coolant, 5 air)
    assert outer.time['lookup'] >= outer.time['coolprop'] > 0

    assert (fluid_properties._fluid_lookup, fluid_properties._air_lookup) == lookups_before    # cache switch restored
    fluid_properties.clear_cache()
    solvers.run_solver('lmtd', ops, geo, settings, profile=False)           # outside a recorder: runs, counts nothing
    with profiling.record() as inner:
        solvers.run_solver('lmtd', ops, geo, settings, profile=False)       # warm cache: no CoolProp call at all
    assert inner.calls['coolprop'] == 0 and inner.calls['lookup'] == 38


def test_solver_figures_are_saved_with_their_titles(tmp_path):
    parameters.T_AIR_IN = 27.0
    run_scenarios(modes(outputs=['cooler_results', 'convergence']))
    saved = {p.name for p in (tmp_path / "figures" / "ch4_results" / "solver").glob("*.png")}
    assert saved == {'cooler_results_dry.png', 'cooler_results_precooled.png', 'convergence_dry.png', 'convergence_precooled.png'}


def test_shared_helpers():
    from src.formatting import fmt_duration, fmt_ms
    from src import weather
    assert fmt_ms(0.04) == "0.04 ms" and fmt_ms(22632.14) == "22632.1 ms"
    assert fmt_duration(0.5) == "500 ms" and fmt_duration(90) == "90.0 s" and fmt_duration(600) == "10.0 min"
    states = weather.hourly_states(weather.load_year(), stride=24)
    assert len(states) == 365 and set(states[0][1]) == {'T_a_i', 'phi', 'p_a'} and 0 < states[0][1]['phi'] <= 1


# --- GUI -------------------------------------------------------------------------------------------------------

def _app():
    from pathlib import Path
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(str(Path(__file__).resolve().parent.parent / "gui.py"), default_timeout=250).run()


def test_gui_loads_in_both_modes_with_the_catalogue_selections():
    app = _app()
    assert not app.exception
    assert [e.label for e in app.sidebar.expander][:3] == ['Run', 'Operating conditions', 'Coolant']
    assert {m.label: m.value for m in app.multiselect}['figures and benchmarks'] == ['cooler_results', 'caching', 'resolution']
    app.sidebar.radio[0].set_value("Cooler comparison").run()
    assert not app.exception
    assert {m.label: m.value for m in app.multiselect}['analyses'] == ['geometry']


def test_gui_flow_given_as_volume_flow_gives_the_same_result():
    """The 'exactly one nonzero' triple is ONE choice + one number; 1 kg/s as 0.001007 m3/s is the same cooler run."""
    def run(basis):
        app = _app()
        if basis:                                                  # outside the form: reruns at once -- before the form's widgets are set
            app.selectbox(key="flow_basis_water flow").set_value(basis).run()
        next(m for m in app.multiselect if m.label == 'figures and benchmarks').set_value([])
        next(m for m in app.multiselect if m.label == 'solvers that run').set_value(['lmtd'])
        next(b for b in app.button if 'Run' in b.label).click().run()
        assert not app.exception
        return app.markdown[0].value.split('\n')[2]                       # the first result row of the table
    row_mass = run(None)
    row_volume = run('V_COOLANT')
    assert row_mass.split('|')[3:6] == row_volume.split('|')[3:6]         # Q, T_c,o, T_a,o identical to the printed digits
    assert 'LMTD' in row_mass


def test_resolution_benchmark_times_the_table_build_up_to_the_maximum():
    resolution = run_scenarios(modes(outputs=['resolution'], benchmark_resolutions=[3, 15, 20], benchmark_table_max_n=15,
                                     solvers=['ntu_table']))['outputs']['resolution']
    assert [n for n, _ in resolution['table_build']] == [3, 15] and all(t > 0 for _, t in resolution['table_build'])
    assert resolution['ntu_variants'] == 'table' and [row[0] for row in resolution['ntu_field']] == [3, 15, 20]
    # all swept n up to the maximum are timed; the design point has R = C_a/C_c ~ 1.4: the table of 3 elements (R up to 0.3)
    # does not cover it (a field solve ran), the one of 15 does
    table = resolution['ntu_table']
    assert [row[0] for row in table] == [3, 15] and resolution['table_covers'] == [15]
    assert abs(table[1][3] - resolution['ntu_field'][1][3]) < 1e-3
