"""Smoke tests: every workflow and benchmark runs end to end on a reduced problem, and the numbers that must not
move without a decision (the fixed tube counts, the independence of the geometry from the operating inputs, the
air-side power) stay put. All of it runs in well under a minute (`venv/bin/python -m pytest`)."""

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
    resolution = results['resolution']
    assert [row[0] for row in resolution['ntu_field']] == [2, 3] and [row[0] for row in resolution['cell_2d']] == [2, 3]
    assert 'lmtd' in resolution
    import matplotlib.pyplot as plt
    titles = {plt.figure(n).canvas.manager.get_window_title() for n in plt.get_fignums()}
    assert {'Caching performance (linear)', 'Caching performance (log)', 'Resolution performance'} <= titles


def test_annual_benchmark_runs_and_extrapolates():
    results = run_scenarios(modes(outputs=['annual'], solvers=['lmtd'], annual_time_budget=1.0))['outputs']['annual']
    lmtd = results['solvers']['LMTD']
    assert lmtd['n_done'] > 0 and lmtd['total'] > 0
    assert lmtd['complete'] == (lmtd['n_done'] == 8760)


# --- cooler comparison ---------------------------------------------------------------------------------------

FAST = ['geometry', 'ambient_cases', 'airflow_sweep', 'sensitivity', 'rows', 'tubes', 'wall']


@pytest.mark.parametrize("analysis", FAST)
def test_comparison_analysis(analysis, outputs, monkeypatch):
    from comparison import sensitivity
    monkeypatch.setattr(sensitivity, 'N_POINTS', 3)                      # the shape of the study, not its resolution
    from comparison.runner import run_analyses
    results = run_analyses([analysis])
    assert results[analysis] and all(len(df) for df in results[analysis].values())
    assert list((outputs / "figures" / "ch4").glob("*.png"))             # at least one figure was saved


def test_comparison_year_with_fan_control(outputs):
    parameters.COMPARE_YEAR_STRIDE = 48
    from comparison.runner import run_analyses
    data = run_analyses(['fan_control'])['fan_control']
    control = data['control']
    assert set(control['set']) == {'P1', 'P2', 'P3'}
    assert (control['w_fr'] <= parameters.COMPARE_W_FR_DESIGN + 1e-9).all()
    assert (control['w_fr'] >= parameters.COMPARE_W_FR_MIN - 1e-9).all()
    assert data['fixed']['w_fr'].eq(parameters.COMPARE_W_FR_DESIGN).all()


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
    saved = {p.name for p in (tmp_path / "figures" / "solver").glob("*.png")}
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
