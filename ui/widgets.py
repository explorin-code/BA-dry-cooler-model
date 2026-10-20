"""
widgets.py -- one input widget per parameter (chosen by its name / the type of its default) and the sidebar form of a mode.
"""

import ast
import math

import streamlit as st

from comparison.analyses import ANALYSES
from src.inputs import weather
from src.inputs.cooler_sets import CELL_SETS
from src.inputs.geometry import SOLID_CONDUCTIVITIES
from src.thermo.fluid_properties import COOLANTS
from src.workflow.outputs import OUTPUTS
from src.solvers import SOLVER_LABELS
from ui.labels import FLOW_MEMBERS, FLOW_TRIPLES, GROUP_ORDER, OPEN_GROUPS, label
from ui.presets import PRESET_PARAMETERS, cell_defaults, flow_equivalents

CATALOGUES = {"SOLVER_OUTPUTS": OUTPUTS, "COMPARE_ANALYSES": ANALYSES}      # list parameters chosen from a catalogue
PRESET_SELECTORS = ("GEOMETRY_SET", "COMPARE_SENSITIVITY_SET")              # chosen above the form, not in it


def weather_files():
    """Hourly climate files of climate_data/ with a full year (8760 rows) -- choices of the annual benchmark's weather."""
    return sorted(p for p in weather.CLIMATE_DIR.glob('munich_*.csv') if sum(1 for _ in p.open()) - 1 == weather.N_HOURS)


def widget(parameter, preset=None):
    """One input widget for a parameter, chosen by name / type of its default. Returns the value."""
    name, default = parameter.name, parameter.default
    text_label = label(name)
    help_text = f"`{name}` -- {parameter.help}" if parameter.help else f"`{name}`"
    key = f"param_{name}" + (f"_{preset}" if name in PRESET_PARAMETERS else "")
    if name in PRESET_PARAMETERS and default is None:             # None = the preset cell's value (tube count: its reference count)
        default = cell_defaults(preset)[name]
    if name == "COMPARE_SETS":
        return st.multiselect(text_label, list(CELL_SETS), default=default, help=help_text, key=key)
    if name in ("COMPARE_SENSITIVITY_SET", "GEOMETRY_SET"):
        return st.selectbox(text_label, list(CELL_SETS), index=list(CELL_SETS).index(default), help=help_text, key=key,
                            format_func=lambda k: CELL_SETS[k].label)
    if name == "SOLVERS":
        return st.multiselect(text_label, list(SOLVER_LABELS), default=default, help=help_text, key=key,
                              format_func=SOLVER_LABELS.get)
    if name in CATALOGUES:                                        # which outputs / analyses run: one multiselect
        catalogue = CATALOGUES[name]
        return st.multiselect(text_label, list(catalogue), default=default, help=help_text, key=key,
                              format_func=lambda k: catalogue[k].title)
    if name == "COMPARE_SOLVER":
        return st.selectbox(text_label, list(SOLVER_LABELS), index=list(SOLVER_LABELS).index(default), help=help_text,
                            key=key, format_func=SOLVER_LABELS.get)
    if name == "ANNUAL_WEATHER_CSV":                              # None = the past year; or one of the shipped full years
        files = {None: "past year (default)", **{str(p): p.stem for p in weather_files()}}
        return st.selectbox(text_label, list(files), format_func=files.get, help=help_text, key=key)
    if name == "COMPARE_N_TUBES":                                 # None = each cooler's own reference count
        return int(st.number_input(text_label, value=0, min_value=0, step=1, key=key,
                                   help=f"{help_text}  (0 = each cooler's own reference count)")) or None
    if isinstance(default, bool):
        return st.toggle(text_label, value=default, help=help_text, key=key)
    if name.endswith("_MATERIAL"):
        options = list(SOLID_CONDUCTIVITIES)
        return st.selectbox(text_label, options, index=options.index(default), help=help_text, key=key)
    if name == "COOLANT_TYPE":
        options = list(COOLANTS) + ([default] if default not in COOLANTS else [])
        return st.selectbox(text_label, options, index=options.index(default), help=help_text, key=key,
                            format_func=lambda fluid: COOLANTS.get(fluid, fluid))
    if isinstance(default, int):
        return int(st.number_input(text_label, value=default, step=1, help=help_text, key=key))
    if isinstance(default, float):
        step = 10.0 ** (math.floor(math.log10(abs(default))) - 1) if default else 0.1     # one decimal below the leading digit
        if name.endswith("_MM"):                                  # lengths: tenths of a millimetre (3 digits below 1 mm); the value itself stays exact
            display = "%.1f" if abs(default) >= 1 else "%.3g"
        else:
            display = "%g"
        return float(st.number_input(text_label, value=default, format=display, step=step, help=help_text, key=key))
    text = st.text_input(text_label, value="" if default is None else repr(default),            # None, str, list: Python literal
                         help=f"{help_text}  (Python literal, empty = None)", key=key)
    if not text.strip():
        return None
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        st.sidebar.error(f"{text_label}: '{text}' is not a valid value (use e.g. 'Water', [5, 20], 0.3)")
        return None


def flow_choices(mode, defaults):
    """Above the form (so the number fields follow the choice immediately): in what form each flow is given."""
    if mode != "solver":
        return {}
    choices = {}
    for title, triple in FLOW_TRIPLES.items():
        plain = {name: text for name, text, _ in triple}
        given = next((n for n in plain if defaults[n]), next(iter(plain)))
        choices[title] = st.selectbox(f"{title} given as", list(plain), index=list(plain).index(given),
                                      format_func=plain.get, key=f"flow_basis_{title}",
                                      help="parameters.py: " + ", ".join(plain) + " -- exactly one of them is nonzero")
    return choices


def flow_widget(title, basis, preset, defaults):
    """The number of the chosen flow form (inside the form); the two other forms become 0 (parameters.py: exactly one nonzero)."""
    labels = {name: latex for name, _, latex in FLOW_TRIPLES[title]}
    start = defaults[basis] or flow_equivalents(preset, defaults)[basis]
    value = st.number_input(labels[basis], value=float(start), format="%g", key=f"param_{basis}_{preset}")
    return {name: (float(value) if name == basis else 0.0) for name in labels}


def parameter_form(parameters, mode, preset, flow_bases, defaults):
    """The sidebar form of one mode: groups as expanders (in GROUP_ORDER), one run button. Returns (values, run pressed)."""
    values, rendered = {}, set()
    with st.form(f"form_{mode}"):
        run_top = st.form_submit_button("▶ Run", type="primary", width="stretch", key=f"run_top_{mode}")
        shown = [p for p in parameters if p.tag in (mode, "both") and p.name not in PRESET_SELECTORS]
        groups = list(dict.fromkeys(p.group for p in shown))
        groups.sort(key=lambda g: GROUP_ORDER[mode].index(g) if g in GROUP_ORDER[mode] else len(GROUP_ORDER[mode]))
        for group in groups:
            with st.expander(group, expanded=group in OPEN_GROUPS[mode]):
                current = None
                members = [p for p in shown if p.group == group]
                sections = {p.section for p in members}
                for p in members:
                    if p.section != current and len(sections) > 1:           # subheadings only if needed
                        st.markdown(f"**{p.section.split(' -- ')[0]}**")
                        current = p.section
                    if p.name in FLOW_MEMBERS:                               # one control for the whole triple
                        title = FLOW_MEMBERS[p.name]
                        if title not in rendered:
                            rendered.add(title)
                            values.update(flow_widget(title, flow_bases[title], preset, defaults))
                        continue
                    values[p.name] = widget(p, preset)
        run_bottom = st.form_submit_button("▶ Run", type="primary", width="stretch", key=f"run_bottom_{mode}")
    return values, run_top or run_bottom
