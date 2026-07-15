"""
fluid_properties.py
====================
Thin wrappers around CoolProp for pulling fluid-state properties
(density, cp, thermal conductivity, viscosity, Prandtl number) at a
given (pressure, temperature) state. Used by both solvers for the
air-side and coolant-side property lookups at the mean bulk temperature.
"""

from dataclasses import dataclass
import CoolProp.CoolProp as CP


@dataclass
class FluidState:
    rho: float                     # density                       [kg/m3]
    cp: float                      # specific heat capacity         [J/kg-K]
    lambda_: float                 # thermal conductivity           [W/m-K]
    eta: float                     # dynamic viscosity              [Pa-s]
    Pr: float                      # Prandtl number                 [-]
    R: float = 0.0                 # relative humidity (set by get_humid_air_properties;
                                    # left at 0.0 for the plain single-phase path above)


def get_fluid_properties(fluid: str, T_celsius: float, P: float) -> FluidState:
    """Single-phase fluid properties at (P, T) via CoolProp's PropsSI."""
    T_kelvin = T_celsius + 273.15

    rho     = CP.PropsSI('D',       'P', P, 'T', T_kelvin, fluid)
    cp      = CP.PropsSI('C',       'P', P, 'T', T_kelvin, fluid)
    lambda_ = CP.PropsSI('L',       'P', P, 'T', T_kelvin, fluid)
    eta     = CP.PropsSI('V',       'P', P, 'T', T_kelvin, fluid)
    Pr      = CP.PropsSI('Prandtl', 'P', P, 'T', T_kelvin, fluid)

    return FluidState(rho, cp, lambda_, eta, Pr)


def get_humid_air_properties(T_celsius: float, P: float, R: float) -> FluidState:
    """Humid-air properties at (P, T, relative humidity) via CoolProp's
    HAPropsSI, for the adiabatic/hybrid cooler modes.

    NOTE: I couldn't run CoolProp in the environment where I wrote this
    (no network access to install it), so this is untested -- please
    sanity-check the numbers against a known point (e.g. HAPropsSI('Vha',
    'T', 293.15, 'P', 101325, 'R', 0.5) should come out around 0.845
    m3/kg dry air) before trusting it in the solver.

    Property mapping, per the CoolProp HAPropsSI key reference
    (http://www.coolprop.org/fluid_properties/HumidAir.html):
      - 'Vha'  : humid-air specific volume, per kg of humid (moist) air [m3/kg]
                 -> density is just the reciprocal, 1/Vha
      - 'Cha'  : humid-air specific heat, per kg of humid air           [J/kg-K]
      - 'K'    : mixture thermal conductivity                          [W/m-K]
      - 'M'    : mixture dynamic viscosity                             [Pa-s]
    HAPropsSI has no direct Prandtl output, so Pr is computed by hand
    from the definition Pr = cp * eta / lambda.

    The previous (broken) attempt at this function is kept below as a
    commented-out reference -- see the note there for what was wrong
    with it.
    """
    T_kelvin = T_celsius + 273.15

    v_ha    = CP.HAPropsSI('Vha', 'T', T_kelvin, 'P', P, 'R', R)  # m3 per kg humid air
    rho     = 1.0 / v_ha
    cp      = CP.HAPropsSI('Cha', 'T', T_kelvin, 'P', P, 'R', R)  # J/kg-K, per kg humid air
    lambda_ = CP.HAPropsSI('K',   'T', T_kelvin, 'P', P, 'R', R)  # W/m-K
    eta     = CP.HAPropsSI('M',   'T', T_kelvin, 'P', P, 'R', R)  # Pa-s
    Pr      = cp * eta / lambda_                                 # no direct HAPropsSI output for this

    return FluidState(rho, cp, lambda_, eta, Pr, R)

# -----------------------------------------------------------------------
# REFERENCE (inactive) -- original attempt at get_humid_air_properties().
# Kept for comparison only -- do not use as-is. It called plain PropsSI
# (not HAPropsSI) for cp/lambda_/eta/Pr, and passed it the 'D' (density)
# key for all four, so it silently returned density under four different
# names instead of the actual properties.
# -----------------------------------------------------------------------
# def get_humid_air_properties(T_celsius: float, P: float, R: float) -> FluidState:
#
#     T_kelvin = T_celsius + 273.15
#
#     rho = CP.HAPropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
#     cp = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
#     lambda_ = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
#     eta = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
#     Pr = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
#
#     return FluidState(rho, cp, lambda_, eta, Pr)