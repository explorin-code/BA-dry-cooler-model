"""
gui.py
======
Streamlit GUI for the dry-cooler model -- the interactive alternative to main.py. Start with:

    streamlit run gui.py

Two modes (toggle at the top of the sidebar), each with its own parameters and run button:

  Solver comparison   one cooler (standard cell P1/P2/P3, any value overridable) in an ambient state; the selected
                      solvers (LMTD, NTU field / table, Cell 2D / 3D) compared, ambient + precooled scenario,
                      solver benchmarks.
  Cooler comparison   the standard coolers raced against each other over the Munich climate year with one solver:
                      geometry, fan control, sensitivity, ...

The forms are built from src/parameters.py itself (values as defaults, inline comments as help, section headers
"# --- Name [solver|compare|both]: ... ---" as groups and mode assignment), so new parameters appear
automatically. A run sets the edited values on the parameters module for this run only -- parameters.py is
never written.

The page is this file; its parts are in ui/ (parameters.py reader, labels, widgets, run, results).
"""

import matplotlib
matplotlib.use("Agg")                     # figures are rendered into the page, never shown
import streamlit as st

import src.parameters as parameters
from src.inputs.cooler_sets import CELL_SETS
from ui.labels import MODES
from ui.parameter_file import read_parameters
from ui.presets import cell_defaults
from ui.results import show_output
from ui.runner import run_model
from ui.widgets import PRESET_SELECTORS, flow_choices, parameter_form

st.set_page_config(page_title="Dry-cooler model", layout="wide", initial_sidebar_state="expanded")
st.title("Dry-cooler model")

parameter_list = read_parameters()
defaults = {p.name: p.default for p in parameter_list}

with st.sidebar:
    mode_name = st.radio("Mode", list(MODES), horizontal=True, label_visibility="collapsed")
    mode = MODES[mode_name]
    preset = sensitivity_set = None
    if mode == "solver":
        st.caption("One cooler in one ambient state: the selected solvers compared, plus the solver benchmarks.")
        preset = st.selectbox("Cooler (standard cell)", list(CELL_SETS), index=list(CELL_SETS).index(parameters.GEOMETRY_SET),
                              format_func=lambda k: CELL_SETS[k].label, key="geometry_set",
                              help="The elementary cell of the run (src/inputs/cooler_sets.py). Every value below can be changed for this run.")
    else:
        st.caption("The standard coolers raced against each other over the climate year with one solver.")
        sensitivity_set = st.selectbox("Cooler for the sensitivity study", list(CELL_SETS),
                                       index=list(CELL_SETS).index(parameters.COMPARE_SENSITIVITY_SET),
                                       format_func=lambda k: CELL_SETS[k].label, key="sensitivity_set",
                                       help="`COMPARE_SENSITIVITY_SET` -- the standard cell whose variables the sensitivity "
                                            "analysis varies one at a time.")
    flow_bases = flow_choices(mode, defaults)
    values, run_pressed = parameter_form(parameter_list, mode, preset, flow_bases, defaults)
if mode == "solver":
    values["GEOMETRY_SET"] = preset
else:
    values["COMPARE_SENSITIVITY_SET"] = sensitivity_set

key = f"output_{mode}"
if run_pressed:
    cell_default = cell_defaults(preset) if mode == "solver" else {}
    changed = {n: v for n, v in values.items() if v != cell_default.get(n, defaults[n]) and n not in PRESET_SELECTORS}
    with st.spinner("Running ... (solver benchmarks ~30 s, year analyses 1-3 min, the stride parameter speeds them up)"):
        try:
            st.session_state[key] = (*run_model(mode, values, defaults), changed)
        except Exception as error:                     # show model errors in the page
            st.session_state[key] = None
            st.error(f"Run failed: {error}")
            st.exception(error)

if st.session_state.get(key) is None:
    st.info(f"{mode_name}: set parameters in the sidebar and press **Run**.")
else:
    show_output(mode, st.session_state[key])
