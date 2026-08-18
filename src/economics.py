"""
economics.py
=============
Cost-relevant physical quantities: pump/fan power consumption and water
usage from precooling. Consumes already-solved OperatingConditions/
FluidState/Geometry -- no solver/precooling internals.
"""


# AI-REVIEW: pump/fan power report hydraulic power only (m_dot * v * delta_p),
# no pump/fan efficiency factor -- confirm this is the intended basis before
# using these numbers for sizing/cost estimates. See CLAUDE.md.
def calc_pump_power(ops, coolant_state, delta_p: float) -> float:
    """Coolant pump power consumption [W]."""
    v = 1.0 / coolant_state.rho
    return ops.m_coolant * v * delta_p  # W


def calc_fan_power(ops, air_state, delta_p: float) -> float:
    """Fan power consumption [W]."""
    v = 1.0 / air_state.rho
    return ops.m_o * v * delta_p  # W


def calc_total_power(W_pump, W_fan):
    """Total electrical power consumption [W]. Returns None if either
    component isn't available yet."""
    if W_pump is None or W_fan is None:
        return None
    return W_pump + W_fan


def calc_water_usage(ops_ambient, ops_precooled) -> float:
    """Water consumption rate [kg/s] of the precooling pad."""
    m_a = ops_ambient.m_o
    X_E = ops_ambient.X_air
    X_A = ops_precooled.X_air
    return m_a * (X_A - X_E)