"""
fluid_properties.py
====================
CoolProp property lookups for the coolant (PropsSI) and humid air
(HAPropsSI), both cached and both keyed off OperatingConditions.
"""

from dataclasses import dataclass
from functools import lru_cache
import CoolProp.CoolProp as CP


@dataclass
class FluidState:
    rho: float                     # density                       [kg/m3]
    cp: float                      # specific heat capacity         [J/kg-K]
    lambda_: float                 # thermal conductivity           [W/m-K]
    eta: float                     # dynamic viscosity              [Pa-s]
    Pr: float                      # Prandtl number                 [-]
    phi_air: float = None          # local relative humidity [0-1] -- set by
                                    # get_air_properties; None for the coolant path


@lru_cache(maxsize=None)
def _get_fluid_properties_cached(fluid: str, T_rounded: float, P: float) -> FluidState:
    T_kelvin = T_rounded + 273.15
    rho     = CP.PropsSI('D',       'P', P, 'T', T_kelvin, fluid)
    cp      = CP.PropsSI('C',       'P', P, 'T', T_kelvin, fluid)
    lambda_ = CP.PropsSI('L',       'P', P, 'T', T_kelvin, fluid)
    eta     = CP.PropsSI('V',       'P', P, 'T', T_kelvin, fluid)
    Pr      = CP.PropsSI('Prandtl', 'P', P, 'T', T_kelvin, fluid)
    return FluidState(rho, cp, lambda_, eta, Pr)


def get_fluid_properties(ops, T_celsius: float, P: float) -> FluidState:
    """Coolant-side properties at (P, T), using ops.coolant_type."""
    return _get_fluid_properties_cached(ops.coolant_type, round(T_celsius, 2), P)


# AI-REVIEW: HAPropsSI key conventions (Vha/Cha = per kg humid air, not per
# kg dry air) never independently verified against a reference psychrometric
# chart/table. See CLAUDE.md.
@lru_cache(maxsize=None)
def _get_air_properties_cached(X_rounded: float, T_rounded: float, P: float) -> FluidState:
    T_kelvin = T_rounded + 273.15
    v_ha    = CP.HAPropsSI('Vha', 'T', T_kelvin, 'P', P, 'W', X_rounded)
    rho     = 1.0 / v_ha
    cp      = CP.HAPropsSI('Cha', 'T', T_kelvin, 'P', P, 'W', X_rounded)
    lambda_ = CP.HAPropsSI('K',   'T', T_kelvin, 'P', P, 'W', X_rounded)
    eta     = CP.HAPropsSI('M',   'T', T_kelvin, 'P', P, 'W', X_rounded)
    Pr      = cp * eta / lambda_
    phi_air = CP.HAPropsSI('R',   'T', T_kelvin, 'P', P, 'W', X_rounded)
    return FluidState(rho, cp, lambda_, eta, Pr, phi_air)


def get_air_properties(ops, T_celsius: float, P: float) -> FluidState:
    """Humid-air properties at (P, T), using the fixed inlet humidity ratio
    ops.X_air (resolved once in OperatingConditions.__post_init__)."""
    return _get_air_properties_cached(round(ops.X_air, 6), round(T_celsius, 2), P)


def relative_to_absolute_humidity(T_celsius: float, P: float, phi_air: float) -> float:
    """Converts relative humidity [0-1] to humidity ratio [kg water/kg dry
    air] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve phi_air into the canonical X_air."""
    T_kelvin = T_celsius + 273.15
    return CP.HAPropsSI('W', 'T', T_kelvin, 'P', P, 'R', phi_air)


def absolute_to_relative_humidity(T_celsius: float, P: float, X_air: float) -> float:
    """Converts humidity ratio [kg water/kg dry air] to relative humidity
    [0-1] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve X_air into phi_air for display/diagnostics."""
    T_kelvin = T_celsius + 273.15
    return CP.HAPropsSI('R', 'T', T_kelvin, 'P', P, 'W', X_air)