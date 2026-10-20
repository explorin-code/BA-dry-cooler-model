"""Shared test setup: figures are rendered off-screen, the parameters module is restored after every test, and the
cooler comparison writes its CSVs and figures into a temporary folder (never over the real results)."""

import copy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

import src.parameters as parameters


@pytest.fixture(autouse=True)
def restore_parameters(tmp_path, monkeypatch):
    from src.plots import figures
    monkeypatch.setattr(figures, "OUTPUT_DIR", tmp_path / "figures")      # no test writes into thesis_figures/
    saved = {name: copy.deepcopy(value) for name, value in vars(parameters).items() if name.isupper()}
    yield
    for name, value in saved.items():
        setattr(parameters, name, value)
    plt.close("all")


@pytest.fixture
def outputs(tmp_path, monkeypatch):
    """Redirects the comparison's CSV output and every figure (src/plots/figures.py)."""
    from comparison import runner
    from src.plots import figures
    monkeypatch.setattr(runner, "DATA_DIR", tmp_path / "csv")
    monkeypatch.setattr(figures, "OUTPUT_DIR", tmp_path / "figures")
    return tmp_path
