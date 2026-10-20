"""The figure styles are contexts, not global state: drawing a figure must not change what the next one (or the other workflow in
the same process, the GUI) gets."""

import matplotlib.pyplot as plt
import pytest

from src.plots.plot_style import PRINT_STYLE, print_context


def test_print_context_sets_the_print_sizes_and_restores_them():
    before = plt.rcParams['font.size']
    with print_context():
        assert plt.rcParams['font.size'] == PRINT_STYLE['font.size']
        assert plt.rcParams['axes.grid']                       # seaborn's whitegrid look
    assert plt.rcParams['font.size'] == before


@pytest.mark.parametrize("analysis", ['geometry', 'rows'])
def test_comparison_figures_do_not_leak_their_style(analysis, outputs):
    from comparison.runner import run_analyses
    import src.parameters as parameters
    parameters.COMPARE_YEAR_STRIDE = 48
    before = dict(plt.rcParams)
    run_analyses([analysis])
    assert plt.rcParams['font.size'] == before['font.size'] and plt.rcParams['axes.facecolor'] == before['axes.facecolor']


def test_figures_land_in_the_chapter_folders(outputs):
    """Chapter 4: solver / cooler_performance / geometry; the module drawings are not mixed into the cooler folder."""
    from comparison.runner import run_analyses
    run_analyses(['geometry', 'rows'])
    root = outputs / "figures" / "ch4_results"
    assert (root / "geometry" / "geometry_print.png").exists() and (root / "geometry" / "geometries.png").exists()
    assert (root / "cooler_performance" / "rows_analysis.png").exists()
    assert not (root / "cooler_performance" / "geometry_print.png").exists()


def test_talk_proportions_are_an_option_of_the_print_context():
    with print_context():
        plain = plt.rcParams['xtick.major.size']
    with print_context(talk=True):
        talk = plt.rcParams['xtick.major.size']
        assert plt.rcParams['font.size'] == PRINT_STYLE['font.size']       # the print sizes still win
    assert talk > plain


def test_importing_the_workflows_sets_no_global_style():
    """Importing the solver pipeline (and the comparison) must leave matplotlib's defaults alone; every figure sets its own context."""
    import subprocess
    import sys
    code = ("import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; before = plt.rcParams['font.size']; "
            "import src.workflow.scenario_pipeline, comparison.analyses; print(before == plt.rcParams['font.size'])")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "True"
