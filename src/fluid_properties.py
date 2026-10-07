"""
fluid_properties.py
====================
CoolProp property lookups for the coolant (PropsSI) and humid air
(HAPropsSI), both cached and both keyed off OperatingConditions.
Source: Bell, Wronski, Quoilin & Lemort (2014), "Pure and Pseudo-pure
Fluid Thermophysical Property Evaluation and the Open-Source
Thermophysical Property Library CoolProp", Ind. Eng. Chem. Res. 53(6),
2498-2508, DOI 10.1021/ie4033999.
"""

from dataclasses import dataclass
from functools import lru_cache
import CoolProp.CoolProp as CP


# Coolants offered in the GUI (CoolProp fluid names -> description). Any other
# CoolProp name also works in parameters.py. Glycol mixtures are CoolProp
# incompressible fluids, mass fraction in %; all checked at 5-60 °C.
COOLANTS = {
    'Water':           'Water',
    'INCOMP::MEG-20%': 'Ethylene glycol 20 %',
    'INCOMP::MEG-30%': 'Ethylene glycol 30 %',
    'INCOMP::MEG-40%': 'Ethylene glycol 40 %',
    'INCOMP::MEG-50%': 'Ethylene glycol 50 %',
    'INCOMP::MPG-20%': 'Propylene glycol 20 %',
    'INCOMP::MPG-30%': 'Propylene glycol 30 %',
    'INCOMP::MPG-40%': 'Propylene glycol 40 %',
    'INCOMP::MPG-50%': 'Propylene glycol 50 %',
}


@dataclass
class FluidState:
    rho: float                     # density                       [kg/m3]
    cp: float                      # specific heat capacity         [J/kg-K]
    lambda_: float                 # thermal conductivity           [W/m-K]
    mu: float                      # dynamic viscosity              [Pa-s]
    Pr: float                      # Prandtl number                 [-]
    phi: float = None              # local relative humidity [0-1] -- set by
                                    # get_air_properties; None for the coolant path


@lru_cache(maxsize=None)
def _get_fluid_properties_cached(fluid: str, T_rounded: float, P: float) -> FluidState:
    T_kelvin = T_rounded + 273.15
    rho     = CP.PropsSI('D',       'P', P, 'T', T_kelvin, fluid)
    cp      = CP.PropsSI('C',       'P', P, 'T', T_kelvin, fluid)
    lambda_ = CP.PropsSI('L',       'P', P, 'T', T_kelvin, fluid)
    mu      = CP.PropsSI('V',       'P', P, 'T', T_kelvin, fluid)
    Pr      = CP.PropsSI('Prandtl', 'P', P, 'T', T_kelvin, fluid)
    return FluidState(rho, cp, lambda_, mu, Pr)


# Property-cache key resolution: temperatures rounded to 0.01 K. The rounding
# step is a noise floor -- a state on a rounding boundary flips between two
# property sets, a period-2 limit cycle (Cell ~1.5e-5 K; LMTD/NTU up to
# ~1.5e-3 K at cold, humid hours, i.e. above the 1e-3 K threshold). LMTD/NTU's
# outer loop damps such stalled oscillations (solvers._relax_lmtd_ntu); 0.001 K
# would also remove them, but costs Cell ~+31 % in annual runs (CLAUDE.md).
CACHE_T_DECIMALS = 2


def get_fluid_properties(ops, T: float, P: float) -> FluidState:
    """Coolant-side properties at (P, T), using ops.coolant_type."""
    return _get_fluid_properties_cached(ops.coolant_type, round(T, CACHE_T_DECIMALS), P)


# AI-REVIEW: HAPropsSI key conventions (Vha/Cha = per kg humid air, not per
# kg dry air) never independently verified against a reference psychrometric
# chart/table. See CLAUDE.md.
@lru_cache(maxsize=None)
def _get_air_properties_cached(Y_rounded: float, T_rounded: float, P: float) -> FluidState:
    T_kelvin = T_rounded + 273.15
    v_ha    = CP.HAPropsSI('Vha', 'T', T_kelvin, 'P', P, 'W', Y_rounded)
    rho     = 1.0 / v_ha
    cp      = CP.HAPropsSI('Cha', 'T', T_kelvin, 'P', P, 'W', Y_rounded)
    lambda_ = CP.HAPropsSI('K',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    mu      = CP.HAPropsSI('mu',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    Pr      = cp * mu / lambda_
    phi     = CP.HAPropsSI('R',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    return FluidState(rho, cp, lambda_, mu , Pr, phi)


def get_air_properties(ops, T: float, P: float) -> FluidState:
    """Humid-air properties at (P, T), using the fixed inlet humidity ratio
    ops.Y (resolved once in OperatingConditions.__post_init__)."""
    return _get_air_properties_cached(round(ops.Y, 6), round(T, CACHE_T_DECIMALS), P)


def relative_to_absolute_humidity(T: float, P: float, phi: float) -> float:
    """Converts relative humidity [0-1] to humidity ratio [kg water/kg dry
    air] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve phi into the canonical X."""
    T_kelvin = T + 273.15
    return CP.HAPropsSI('W', 'T', T_kelvin, 'P', P, 'R', phi)


def absolute_to_relative_humidity(T: float, P: float, Y: float) -> float:
    """Converts humidity ratio [kg water/kg dry air] to relative humidity
    [0-1] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve X into phi for display/diagnostics."""
    T_kelvin = T + 273.15
    return CP.HAPropsSI('R', 'T', T_kelvin, 'P', P, 'W', Y)