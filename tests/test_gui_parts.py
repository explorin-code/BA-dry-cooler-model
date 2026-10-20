"""The GUI's building blocks without a browser: the parameters.py reader, the labels and groups, the presets, the tables."""

import pandas as pd
import pytest

import src.parameters as parameters
from comparison.analyses import ANALYSES
from src.workflow.outputs import OUTPUTS
from src.solvers import SOLVER_LABELS
from ui import labels
from ui.parameter_file import read_parameters
from ui.presets import cell_defaults, flow_equivalents
from ui.results import markdown_table

PARAMETERS = read_parameters()
DEFAULTS = {p.name: p.default for p in PARAMETERS}


def test_every_constant_of_parameters_py_is_read_once_with_its_default():
    names = [p.name for p in PARAMETERS]
    assert len(names) == len(set(names))
    assert all(DEFAULTS[name] == getattr(parameters, name) for name in names)
    assert {p.tag for p in PARAMETERS} <= {"solver", "compare", "both"}


def test_help_texts_come_from_the_comments():
    help_of = {p.name: p.help for p in PARAMETERS}
    assert "face velocity" in help_of["COMPARE_W_FR_DESIGN"] or "fan speed" in help_of["COMPARE_W_FR_DESIGN"]
    assert all(help_of[name] for name in ("SOLVERS", "T_COOLANT_IN", "COMPARE_ANALYSES"))
    assert all(len(p.help) < 260 for p in PARAMETERS), "a help text that long does not fit the sidebar tooltip"


def test_groups_labels_and_flow_triples_are_consistent():
    for mode in ("solver", "compare"):
        shown = {p.group for p in PARAMETERS if p.tag in (mode, "both")}
        assert shown <= set(labels.GROUP_ORDER[mode]), f"a group missing from GROUP_ORDER[{mode}]: {shown - set(labels.GROUP_ORDER[mode])}"
        assert set(labels.OPEN_GROUPS[mode]) <= set(labels.GROUP_ORDER[mode])
    selectors = {"GEOMETRY_SET", "COMPARE_SENSITIVITY_SET"}
    assert [n for n in DEFAULTS if n not in labels.LABELS and n not in selectors] == []
    assert set(labels.LABELS) <= set(DEFAULTS)
    assert set(labels.FLOW_MEMBERS) <= set(DEFAULTS)
    for triple in labels.FLOW_TRIPLES.values():                  # exactly one of each triple is nonzero in the file
        assert sum(bool(DEFAULTS[name]) for name, _, _ in triple) == 1


def test_catalogue_and_solver_selections_use_valid_keys():
    assert set(DEFAULTS["SOLVER_OUTPUTS"]) <= set(OUTPUTS)
    assert set(DEFAULTS["COMPARE_ANALYSES"]) <= set(ANALYSES)
    assert set(DEFAULTS["SOLVERS"]) <= set(SOLVER_LABELS)


def test_presets_show_the_standard_coolers_values():
    p3 = cell_defaults("P3")
    assert p3["N_TUBES"] == 30 and p3["D_TUBE_OUTER_MM"] == 16.3 and p3["FIN_PITCH_MM"] == pytest.approx(4.2, abs=0.01)
    assert cell_defaults("P1")["N_TUBES"] == 34


def test_flow_equivalents_describe_the_same_flow():
    flows = flow_equivalents("P2", DEFAULTS)
    assert flows["M_COOLANT"] == pytest.approx(1.0) and flows["V_COOLANT"] == pytest.approx(1.0 / 993, rel=0.02)
    assert flows["W_O"] == pytest.approx(DEFAULTS["W_O"])


def test_markdown_table_formats_numbers_and_gaps():
    table = markdown_table(pd.DataFrame({"a": [1.2345, None], "b": ["x", "y"], "c": [3, 4]}), decimals=2)
    assert table.splitlines() == ["| a | b | c |", "|---|---|---|", "| 1.23 | x | 3 |", "| – | y | 4 |"]
