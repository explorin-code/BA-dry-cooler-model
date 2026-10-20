"""
common.py -- what every solver does first: fill in the defaults, evaluate the states at the mean temperatures.
"""

from src.inputs.geometry import get_geometry
from src.thermo.fluid_properties import get_air_properties, get_fluid_properties
from src.inputs.operating_conditions import get_operating_conditions
from src.inputs.solver_settings import get_solver_settings


def fill_defaults(ops, geo, settings):
    """Fills in whatever the caller didn't pass."""
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    return ops, geo, settings


def mean_states(ops, T_c_o, T_a_o):
    """Coolant and air states at their mean temperatures (inlet + outlet)/2."""
    coolant_state = get_fluid_properties(ops, (ops.T_c_i + T_c_o) / 2.0, ops.p_c)
    air_state = get_air_properties(ops, (ops.T_a_i + T_a_o) / 2.0, ops.p_a)
    return coolant_state, air_state
