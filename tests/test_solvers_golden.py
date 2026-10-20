"""The solvers' numbers, pinned: Q, both outlet temperatures, k and the iteration count of each of the five solvers at
a design, a hot and the coldest operating point of the P2 cooler, plus values of the NTU table (tests/golden_solvers.json,
recorded from the code before solvers.py became a package -- the split reproduced them bit for bit), and the literature
validation of the NTU element field (Cabezas-Gomez et al. 2007, Tables 1 and 5). A change here is a change of the model."""

import json
from pathlib import Path

import pytest

from src import solvers
from src.dry_cooler_physics import get_geometry
from src.operating_conditions import get_operating_conditions
from src.solver_settings import get_solver_settings
from src.solvers import ntu_table

GOLDEN = json.loads((Path(__file__).parent / "golden_solvers.json").read_text())
POINTS = {'design': dict(T_a_i=20.0, phi=0.6, p_a=95500.0), 'hot': dict(T_a_i=32.0, phi=0.3, p_a=95500.0),
          'coldest': dict(T_a_i=-10.7, phi=0.88, p_a=95230.0)}


@pytest.mark.parametrize("point", POINTS)
@pytest.mark.parametrize("key", list(solvers.SOLVER_LABELS))
def test_solver_results_are_unchanged(point, key):
    geo = get_geometry()
    ops = get_operating_conditions(geo, **POINTS[point])
    result = solvers.run_solver(key, ops, geo, get_solver_settings(), profile=False)
    expected = GOLDEN[f"{point}/{key}"]
    assert len(result.history_hot) == expected['iterations']
    for name in ('Q_dot', 'T_c_o', 'T_a_o', 'k'):
        assert getattr(result, name) == pytest.approx(expected[name], rel=1e-9), name


def test_ntu_table_values_are_unchanged():
    geo, settings = get_geometry(), get_solver_settings()
    table = ntu_table.get_ntu_table(geo.N_r, settings.ntu_n_elements)
    for key, expected in GOLDEN['table'].items():
        NTU, R = map(float, key.split('/'))
        assert float(table(NTU, R)[0, 0]) == pytest.approx(expected, rel=1e-9), key


def test_benchmark_table_is_the_first_table():
    """The solver benchmarks build their NTU table over NTU_a <= 6 (ntu_table.NTU_TABLE_NTU_BENCHMARK): it reproduces, bit for bit up to the
    interpolation, the values recorded when that was the only table -- the build times of the thesis belong to it."""
    geo, settings = get_geometry(), get_solver_settings()
    axis = ntu_table.NTU_TABLE_NTU_BENCHMARK
    from src.benchmark_solvers import preserved_ntu_tables
    with preserved_ntu_tables():                                    # the narrow table must not stay installed for the other tests
        ntu_table.install_table(geo.N_r, settings.ntu_n_elements, ntu_table.build_ntu_table(geo.N_r, settings.ntu_n_elements, axis), axis)
        table = ntu_table.tables[(geo.N_r, settings.ntu_n_elements)]
        for key, expected in GOLDEN['table_benchmark'].items():
            NTU, R = map(float, key.split('/'))
            assert float(table(NTU, R)[0, 0]) == pytest.approx(expected, rel=1e-9), key


def test_ntu_element_field_matches_the_literature():
    from src.validation_ntu import run_validation
    assert run_validation()


@pytest.mark.parametrize("key", ['lmtd', 'ntu_field', 'cell_2d'])
def test_profile_heat_flow_adds_up_to_q(key):
    """The dQ/dL curve of the cooler-results figure integrates to the solver's Q (Cell's once did not: 89 %)."""
    import numpy as np
    from src.analysis import calc_cell_profile, calc_lmtd_profile, calc_ntu_profile
    geo, settings = get_geometry(), get_solver_settings()
    ops = get_operating_conditions(geo, **POINTS['design'])
    result = solvers.run_solver(key, ops, geo, settings, profile=False)
    profile = {'lmtd': lambda: calc_lmtd_profile(result, ops, geo),
               'ntu_field': lambda: calc_ntu_profile(result, geo, settings.ntu_n_elements),
               'cell_2d': lambda: calc_cell_profile(result, geo, settings.cell_n_segments)}[key]()
    assert np.trapezoid(profile.dQdL, profile.x_frac) * geo.l_tot == pytest.approx(result.Q_dot, rel=2e-3)


def test_whole_number_inputs_give_the_same_result_as_floats():
    """A temperature typed as 80 instead of 80.0 once made Cell's grids integer (every temperature truncated)."""
    geo, settings = get_geometry(), get_solver_settings()
    results = [solvers.run_solver('cell_2d', get_operating_conditions(geo, T_c_i=T_c_i, T_a_i=T_a_i, phi=0.6, p_a=101325.0),
                                  geo, settings, profile=False) for T_c_i, T_a_i in ((60, 20), (60.0, 20.0))]
    assert results[0].Q_dot == results[1].Q_dot and results[0].T_c_o == results[1].T_c_o


def test_registry_and_reference_order():
    assert set(solvers.SOLVER_LABELS) == set(solvers.scenario.SOLVERS) == set(solvers.SOLVER_FAMILY)
    assert solvers.scenario_keys(['cell_2d', 'ntu_table', 'lmtd']) == ['lmtd', 'ntu_table', 'cell_2d']
    with pytest.raises(ValueError, match="Unknown solver"):
        solvers.run_solver('nonsense', None, None, None)
