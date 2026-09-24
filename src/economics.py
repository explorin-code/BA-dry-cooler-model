"""
economics.py
=============
Cost-relevant physical quantities: pump/fan power consumption and water
usage from precooling. Consumes already-solved OperatingConditions/
FluidState/Geometry -- no solver/precooling internals.
"""


def calc_pump_power(ops, coolant_state, delta_p: float) -> float:
    """Coolant pump power consumption [W]."""
    eta_p = 0.9 # pump efficiency -- Towler2022 - 20.7 (mass * specific volume instead of volume flow rate, but same units)
    return (ops.m_dot_1 * delta_p) / (eta_p * coolant_state.rho)  # W
    

def calc_fan_power(ops, air_in, air_out, delta_p: float, fan_position: str = 'forced') -> float:
    """Fan power consumption [W]. fan_position decides which air state the
    fan actually moves (sets the volume flow m_dot/rho):
      'forced'  -- fan upstream of the bundle, pushes air at the bundle inlet state
      'induced' -- fan downstream of the bundle, pulls air at the bundle outlet state"""
    eta_f = 0.7  # Trafan efficiency -- Towler2022 - 19.16 (mass * specific volume instead of velocity * area, but same units)
    if fan_position == 'forced':
        rho = air_in.rho
    elif fan_position == 'induced':
        rho = air_out.rho
    else:
        raise ValueError(f"fan_position must be 'forced' or 'induced', got {fan_position!r}")
    return (ops.m_dot_2 * delta_p) / (eta_f * rho)  # W


def calc_total_power(P_p, P_f):
    """Total electrical power consumption [W]. Returns None if either
    component isn't available yet."""
    if P_p is None or P_f is None:
        return None
    return P_p + P_f


def calc_water_usage(ops_ambient, ops_precooled) -> float:
    """Water consumption rate [kg/s] of the precooling pad."""
    m_a = ops_ambient.m_dot_2
    X_E = ops_ambient.X
    X_A = ops_precooled.X
    return m_a * (X_A - X_E)
