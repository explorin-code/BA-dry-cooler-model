"""
gui.py
======
Streamlit GUI for the dry-cooler model -- the interactive alternative to
main.py. Start with:

    streamlit run gui.py

The form is built from src/parameters.py itself (values as defaults, inline
comments as help, "# --- Group: ... ---" headers as groups), so new
parameters appear automatically. A run sets the edited values on the
parameters module for this session only -- parameters.py is never written.
"""

import ast
import io
import logging
import math
import re
import threading
import tokenize
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                     # figures are rendered into the page, never shown
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import src.parameters as parameters
from src.dry_cooler_physics import SOLID_CONDUCTIVITIES
from src.fluid_properties import COOLANTS
from src.run_modes import get_run_modes
from src.scenario_pipeline import run_scenarios

PARAMETERS_FILE = Path(__file__).parent / "src" / "parameters.py"
SECTION = re.compile(r"^# --- (.*?) -*$")


# =============================================================================
# Reading parameters.py
# =============================================================================

def read_parameters():
    """[(group, section, name, default, help)] in file order. Group = text
    before ':' in the section header ("Geometry: tube/fin" -> "Geometry")."""
    source = PARAMETERS_FILE.read_text()
    lines = source.splitlines()

    comments = {}                                  # line number -> comment text
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok.type == tokenize.COMMENT:
            comments[tok.start[0]] = tok.string.lstrip("# ").strip()

    entries = []
    for node in ast.parse(source).body:
        if not (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)):
            continue
        # Section: nearest "# --- ... ---" header above the assignment.
        section = next((m.group(1).strip() for l in reversed(lines[:node.lineno - 1])
                        if (m := SECTION.match(l))), "Other")
        # Help: inline comment plus indented comment-only continuation lines.
        help_parts = [comments.get(node.lineno, "")]
        n = node.end_lineno + 1
        while n <= len(lines) and lines[n - 1].startswith(" ") and lines[n - 1].strip().startswith("#"):
            help_parts.append(comments[n])
            n += 1
        entries.append((section.split(":")[0].strip(), section, node.targets[0].id,
                        ast.literal_eval(node.value), " ".join(p for p in help_parts if p)))
    return entries


# =============================================================================
# Widgets
# =============================================================================

def widget(name, default, help_text):
    """One input widget, chosen by the default's type. Returns the value."""
    key = f"param_{name}"
    if isinstance(default, bool):
        return st.checkbox(name, value=default, help=help_text, key=key)
    if name.endswith("_MATERIAL"):
        options = list(SOLID_CONDUCTIVITIES)
        return st.selectbox(name, options, index=options.index(default), help=help_text, key=key)
    if name == "COOLANT_TYPE":
        options = list(COOLANTS) + ([default] if default not in COOLANTS else [])   # keep a custom name
        return st.selectbox(name, options, index=options.index(default), help=help_text, key=key,
                            format_func=lambda fluid: COOLANTS.get(fluid, fluid))
    if isinstance(default, int):
        return int(st.number_input(name, value=default, step=1, help=help_text, key=key))
    if isinstance(default, float):
        # step: one decimal place below the default's leading digit (9.0 -> 0.1, 0.12 -> 0.01)
        step = 10.0 ** (math.floor(math.log10(abs(default))) - 1) if default else 0.1
        return float(st.number_input(name, value=default, format="%g", step=step, help=help_text, key=key))
    # None, str, list: free text, parsed as a Python literal (empty = None)
    text = st.text_input(name, value="" if default is None else repr(default),
                         help=f"{help_text}  (Python literal, empty = None)", key=key)
    return parse_literal(name, text)


def parse_literal(name, text):
    if not text.strip():
        return None
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        st.sidebar.error(f"{name}: '{text}' is not a valid value (use e.g. 'Water', [5, 20], 0.3)")
        return None


# =============================================================================
# Running the model
# =============================================================================

class _Capture(logging.Handler):
    """Collects the model's log output (summary / insight lines, warnings)."""
    def __init__(self):
        super().__init__()
        self.lines = []

    def emit(self, record):
        self.lines.append(self.format(record))


@st.cache_resource
def _model_lock():
    """One lock for the whole server: all browser tabs share the model's
    modules (parameter overrides, the benchmark's temporary patches), so
    only one run may execute at a time."""
    return threading.Lock()


def run_model(values):
    """Runs the model with the GUI values (one run at a time across tabs)."""
    lock = _model_lock()
    if not lock.acquire(blocking=False):
        st.warning("Another run (other browser tab) is still in progress -- waiting for it ...")
        lock.acquire()
    try:
        return _run_model(values)
    finally:
        lock.release()


def _run_model(values):
    """Sets the GUI values on the parameters module, runs the pipeline,
    returns (results, log text, [(title, figure)])."""
    for name, value in values.items():
        setattr(parameters, name, value)
    modes = get_run_modes()

    capture = _Capture()
    capture.setFormatter(logging.Formatter("%(message)s"))
    src_logger, warn_logger = logging.getLogger("src"), logging.getLogger("py.warnings")
    src_logger.setLevel(logging.DEBUG if modes.insight_mode else logging.INFO)
    logging.captureWarnings(True)
    for logger in (src_logger, warn_logger):
        logger.addHandler(capture)

    plt.close("all")
    try:
        results = run_scenarios(modes)
    finally:
        for logger in (src_logger, warn_logger):
            logger.removeHandler(capture)

    figures = []
    for number in plt.get_fignums():
        fig = plt.figure(number)
        figures.append((fig.canvas.manager.get_window_title() if fig.canvas.manager else f"Figure {number}", fig))
    return results, "\n".join(capture.lines), figures


def summary_table(results):
    """One row per scenario and solver: the key numbers of the run."""
    rows = []
    for key in ("ambient", "precooled"):
        sc = results.get(key)
        if not sc:
            continue
        for name in ("lmtd", "ntu", "cell"):
            r = getattr(sc["result"], name)
            rows.append({
                "Scenario": sc["label"], "Solver": name.upper() if name != "cell" else "Cell",
                "Q [kW]": r.Q_dot / 1000, "T_c,o [°C]": r.T_c_o, "T_a,o [°C]": r.T_a_o,
                "ΔT_pinch [K]": r.T_c_o - sc["ops"].T_a_i, "k [W/m²K]": r.k,
                "Iterations": len(r.history_hot), "Time [ms]": (r.solve_time or 0) * 1e3,
            })
    return pd.DataFrame(rows)


def economics_table(results):
    rows = []
    for key in ("ambient", "precooled"):
        sc = results.get(key)
        if sc:
            rows.append({"Scenario": sc["label"], "Pump [W]": sc["P_pump"], "Fan [W]": sc["P_fan"],
                         "Water [g/s]": None if sc["m_dot_ev"] is None else sc["m_dot_ev"] * 1000})
    return pd.DataFrame(rows)


# =============================================================================
# Page
# =============================================================================

st.set_page_config(page_title="Dry-cooler model", layout="wide", initial_sidebar_state="expanded")
st.title("Dry-cooler model")

entries = read_parameters()
groups = list(dict.fromkeys(group for group, *_ in entries))
defaults = {name: default for _, _, name, default, _ in entries}

with st.sidebar:
    st.header("Parameters")
    st.caption("Defaults from src/parameters.py -- changes apply to this session only.")
    values = {}
    with st.form("parameters"):
        run_top = st.form_submit_button("▶ Run model", type="primary", width="stretch", key="run_top")
        for group in ["Run modes"] + [g for g in groups if g != "Run modes"]:
            with st.expander(group, expanded=group in ("Run modes", "Operating conditions")):
                sections = {section for g, section, *_ in entries if g == group}
                current_section = None
                for g, section, name, default, help_text in entries:
                    if g != group:
                        continue
                    if section != current_section and len(sections) > 1:   # subheadings only if needed
                        st.markdown(f"**{section.split(':', 1)[-1].strip()}**")
                        current_section = section
                    values[name] = widget(name, default, help_text)
        run_bottom = st.form_submit_button("▶ Run model", type="primary", width="stretch", key="run_bottom")

if run_top or run_bottom:
    changed = {n: v for n, v in values.items() if v != defaults[n]}
    with st.spinner("Running the model ... (benchmark / resolution modes take ~30 s)"):
        try:
            st.session_state["output"] = (*run_model(values), changed)
        except Exception as error:                     # show model errors in the page
            st.session_state["output"] = None
            st.error(f"Run failed: {error}")
            st.exception(error)

output = st.session_state.get("output")
if output is None:
    st.info("Set parameters in the sidebar and press **Run model**.")
else:
    results, log_text, figures, changed = output
    if changed:
        st.caption("Changed vs. parameters.py: " + ", ".join(f"{n} = {v!r}" for n, v in changed.items()))

    left, right = st.columns([3, 2])
    left.subheader("Results")
    left.dataframe(summary_table(results), hide_index=True, width="stretch",
                   column_config={c: st.column_config.NumberColumn(format="%.2f")
                                  for c in ("Q [kW]", "T_c,o [°C]", "T_a,o [°C]", "ΔT_pinch [K]", "k [W/m²K]")} |
                                 {"Time [ms]": st.column_config.NumberColumn(format="%.1f")})
    right.subheader("Economics")
    right.dataframe(economics_table(results), hide_index=True, width="stretch",
                    column_config={c: st.column_config.NumberColumn(format="%.2f")
                                   for c in ("Pump [W]", "Fan [W]", "Water [g/s]")})

    with st.expander("Terminal output", expanded=False):
        st.code(log_text or "(no output)", language=None)

    if figures:
        for tab, (title, fig) in zip(st.tabs([title for title, _ in figures]), figures):
            with tab:
                st.pyplot(fig, width="stretch")
    else:
        st.caption("No figures -- enable PLOT_RESULTS / PLOT_CONVERGENCE / BENCHMARK_MODE / RESOLUTION_MODE.")
