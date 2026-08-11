"""
fluid_properties.py
====================
CoolProp property lookups for the coolant (PropsSI) and humid air (HAPropsSI).
"""

from dataclasses import dataclass
import CoolProp.CoolProp as CP


@dataclass
class FluidState:
    rho: float
    cp: float
    lambda_: float
    eta: float
    Pr: float
    phi_air: float = None


def get_fluid_properties(ops, T_celsius: float, P: float) -> FluidState:
    T_kelvin = T_celsius + 273.15
    fluid = ops.coolant_type

    rho     = CP.PropsSI('D',       'P', P, 'T', T_kelvin, fluid)
    cp      = CP.PropsSI('C',       'P', P, 'T', T_kelvin, fluid)
    lambda_ = CP.PropsSI('L',       'P', P, 'T', T_kelvin, fluid)
    eta     = CP.PropsSI('V',       'P', P, 'T', T_kelvin, fluid)
    Pr      = CP.PropsSI('Prandtl', 'P', P, 'T', T_kelvin, fluid)

    return FluidState(rho, cp, lambda_, eta, Pr)


def get_air_properties(ops, T_celsius: float, P: float) -> FluidState:
    T_kelvin = T_celsius + 273.15
    X = ops.X_air

    v_ha    = CP.HAPropsSI('Vha', 'T', T_kelvin, 'P', P, 'W', X)
    rho     = 1.0 / v_ha
    cp      = CP.HAPropsSI('Cha', 'T', T_kelvin, 'P', P, 'W', X)
    lambda_ = CP.HAPropsSI('K',   'T', T_kelvin, 'P', P, 'W', X)
    eta     = CP.HAPropsSI('M',   'T', T_kelvin, 'P', P, 'W', X)
    Pr      = cp * eta / lambda_
    phi_air = CP.HAPropsSI('R',   'T', T_kelvin, 'P', P, 'W', X)

    return FluidState(rho, cp, lambda_, eta, Pr, phi_air)


def relative_to_absolute_humidity(T_celsius: float, P: float, phi_air: float) -> float:
    T_kelvin = T_celsius + 273.15
    return CP.HAPropsSI('W', 'T', T_kelvin, 'P', P, 'R', phi_air)

def absolute_to_relative_humidity(T_celsius: float, P: float, X_air: float) -> float:
    T_kelvin = T_celsius + 273.15
    return CP.HAPropsSI('R', 'T', T_kelvin, 'P', P, 'W', X_air)