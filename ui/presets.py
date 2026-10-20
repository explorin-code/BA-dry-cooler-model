"""
presets.py -- what "by default" means on the chosen standard cooler: the cell values of the preset (a None cell parameter), its
reference tube count, and the equivalents of the file's flows in every form of FLOW_TRIPLES.
"""

import streamlit as st

import src.parameters as parameters
from src.inputs.cooler_sets import CELL_SETS
from src.inputs.geometry import CELL_PARAMETERS, get_cell_values

PRESET_PARAMETERS = (*CELL_PARAMETERS, 'N_TUBES')           # None = what the standard cooler says; shown with that value


def cell_defaults(preset):
    """{parameters constant: value} of the preset cell (what a None cell parameter means)."""
    saved = parameters.GEOMETRY_SET, {c: getattr(parameters, c) for c in CELL_PARAMETERS}
    parameters.GEOMETRY_SET = preset
    for c in CELL_PARAMETERS:
        setattr(parameters, c, None)
    try:
        cell = get_cell_values()
    finally:
        parameters.GEOMETRY_SET = saved[0]
        for c, v in saved[1].items():
            setattr(parameters, c, v)
    return {**{c: cell[attribute] for c, attribute in CELL_PARAMETERS.items()}, 'N_TUBES': CELL_SETS[preset].N_t_ref}


@st.cache_data
def flow_equivalents(preset, _defaults):
    """The flows of parameters.py's defaults on the preset cooler, in every form of FLOW_TRIPLES: {constant: value} --
    what the number field shows when another form is chosen. _defaults: {name: default} of parameters.py (not hashed)."""
    from src.inputs.geometry import get_geometry
    from src.inputs.operating_conditions import get_operating_conditions
    saved = {name: getattr(parameters, name) for name in _defaults}
    try:
        for name, value in _defaults.items():
            setattr(parameters, name, value)
        parameters.GEOMETRY_SET = preset
        ops = get_operating_conditions(geo=get_geometry())
    finally:
        for name, value in saved.items():
            setattr(parameters, name, value)
    return {"M_COOLANT": ops.m_dot_c, "W_COOLANT": ops.w_c, "V_COOLANT": ops.V_dot_c,
            "W_O": ops.w_fr, "M_O": ops.m_dot_a, "V_O": ops.V_dot_a}
