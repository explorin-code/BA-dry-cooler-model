"""
economics.py
=============
Cost-relevant physical quantities: pump/fan power consumption and water
usage from precooling. Consumes already-solved OperatingConditions/
FluidState/Geometry -- no solver/precooling internals.
"""


def calc_pump_power(ops, coolant_state, Delta_p: float) -> float:
    """Coolant pump power consumption [W]."""
    eta_pump = 0.9  # reciprocating pump efficiency -- Towler, Chemical Engineering Design, 3rd ed., Sec. 20.7
                 # (written with m_dot/rho instead of the volume flow rate -- same quantity)
    return (ops.m_dot_c * Delta_p) / (eta_pump * coolant_state.rho)  # W
    

def calc_fan_power(ops, air_in, air_out, Delta_p: float, fan_position: str = 'forced') -> float:
    """Fan power consumption [W]. fan_position decides which air state the
    fan actually moves (sets the volume flow m_dot/rho):
      'forced'  -- fan upstream of the bundle, pushes air at the bundle inlet state
      'induced' -- fan downstream of the bundle, pulls air at the bundle outlet state"""
    eta_fan = 0.7  # typical fan efficiency -- Towler, Chemical Engineering Design, 3rd ed., Sec. 19.16
                 # (written with m_dot/rho instead of face velocity * face area -- same quantity)
    if fan_position == 'forced':
        rho = air_in.rho
    elif fan_position == 'induced':
        rho = air_out.rho
    else:
        raise ValueError(f"fan_position must be 'forced' or 'induced', got {fan_position!r}")
    return (ops.m_dot_a * Delta_p) / (eta_fan * rho)  # W


def calc_total_power(P_pump, P_fan):
    """Total electrical power consumption [W]. Returns None if either
    component isn't available yet."""
    if P_pump is None or P_fan is None:
        return None
    return P_pump + P_fan


def calc_water_usage(ops_ambient, ops_precooled) -> float:
    """Water consumption rate [kg/s] of the precooling pad: dry-air mass
    flow times the humidity-ratio increase (X is per kg DRY air; ops.m_dot_a
    is the humid-air mass flow, so the dry-air share is m_dot_a / (1 + X)).
    Source: mass balance, VDI-Wärmeatlas (2019) M8 (degree of
    humidification definition). Fixed 2026-10: previously multiplied the
    humid-air mass flow by delta X (~0.4 % too high)."""
    m_dot_da = ops_ambient.m_dot_a / (1.0 + ops_ambient.Y)
    Y_i = ops_ambient.Y
    Y_pc = ops_precooled.Y
    return m_dot_da * (Y_pc - Y_i)
