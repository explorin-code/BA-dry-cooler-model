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

from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
import CoolProp.CoolProp as CP

from src.profiling import timed

_props = timed('coolprop')(CP.PropsSI)            # every CoolProp call of the model goes through these two
_ha_props = timed('coolprop')(CP.HAPropsSI)       # (counted by profiling.record, nothing else)


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
    rho     = _props('D',       'P', P, 'T', T_kelvin, fluid)
    cp      = _props('C',       'P', P, 'T', T_kelvin, fluid)
    lambda_ = _props('L',       'P', P, 'T', T_kelvin, fluid)
    mu      = _props('V',       'P', P, 'T', T_kelvin, fluid)
    Pr      = _props('Prandtl', 'P', P, 'T', T_kelvin, fluid)
    return FluidState(rho, cp, lambda_, mu, Pr)


# Property-cache key resolution: temperatures rounded to 0.01 K. The rounding
# step is a noise floor -- a state on a rounding boundary flips between two
# property sets, a period-2 limit cycle (Cell ~1.5e-5 K; LMTD/NTU up to
# ~1.5e-3 K at cold, humid hours, i.e. above the 1e-3 K threshold). LMTD/NTU's
# outer loop damps such stalled oscillations (solvers/relaxation.py); 0.001 K
# would also remove them, but costs Cell ~+31 % in annual runs (docs/model-notes.md).
CACHE_T_DECIMALS = 2


WATER_T_MIN = 0.01      # [°C] CoolProp's water has no supercooled liquid below the melting point


@timed('lookup')
def get_fluid_properties(ops, T: float, P: float) -> FluidState:
    """Coolant-side properties at (P, T), using ops.coolant_type. Pure water is
    clamped to >= WATER_T_MIN (properties of the 0 °C liquid are used below) so a
    solver that iterates through/below 0 °C still runs -- the caller must treat
    such a result as 'water freezes' (studies/evaluate.py flags T_c,o < 0 °C),
    the numbers there are not physical."""
    if ops.coolant_type == 'Water':
        T = max(T, WATER_T_MIN)
    return _fluid_lookup(ops.coolant_type, round(T, CACHE_T_DECIMALS), P)


# HAPropsSI key conventions: Vha/Cha are per kg humid air, not per kg dry air (docs/model-notes.md).
@lru_cache(maxsize=None)
def _get_air_properties_cached(Y_rounded: float, T_rounded: float, P: float) -> FluidState:
    T_kelvin = T_rounded + 273.15
    v_ha    = _ha_props('Vha', 'T', T_kelvin, 'P', P, 'W', Y_rounded)
    rho     = 1.0 / v_ha
    cp      = _ha_props('Cha', 'T', T_kelvin, 'P', P, 'W', Y_rounded)
    lambda_ = _ha_props('K',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    mu      = _ha_props('mu',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    Pr      = cp * mu / lambda_
    phi     = _ha_props('R',   'T', T_kelvin, 'P', P, 'W', Y_rounded)
    return FluidState(rho, cp, lambda_, mu , Pr, phi)


# What get_fluid_properties / get_air_properties call: the cached lookups, or (caching(False)) the plain functions behind them.
_fluid_lookup = _get_fluid_properties_cached
_air_lookup = _get_air_properties_cached


def clear_cache():
    """Empties both property caches (a 'cold' start)."""
    _get_fluid_properties_cached.cache_clear()
    _get_air_properties_cached.cache_clear()


@contextmanager
def caching(enabled: bool):
    """Switches the property cache on or off for the block (off: every lookup calls CoolProp; benchmarks)."""
    global _fluid_lookup, _air_lookup
    previous = _fluid_lookup, _air_lookup
    if not enabled:
        _fluid_lookup, _air_lookup = _get_fluid_properties_cached.__wrapped__, _get_air_properties_cached.__wrapped__
    try:
        yield
    finally:
        _fluid_lookup, _air_lookup = previous


@timed('lookup')
def get_air_properties(ops, T: float, P: float) -> FluidState:
    """Humid-air properties at (P, T), using the fixed inlet humidity ratio
    ops.Y (resolved once in OperatingConditions.__post_init__)."""
    return _air_lookup(round(ops.Y, 6), round(T, CACHE_T_DECIMALS), P)


def relative_to_absolute_humidity(T: float, P: float, phi: float) -> float:
    """Converts relative humidity [0-1] to humidity ratio [kg water/kg dry
    air] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve phi into the canonical X."""
    T_kelvin = T + 273.15
    return _ha_props('W', 'T', T_kelvin, 'P', P, 'R', phi)


def absolute_to_relative_humidity(T: float, P: float, Y: float) -> float:
    """Converts humidity ratio [kg water/kg dry air] to relative humidity
    [0-1] via HAPropsSI. Used once by OperatingConditions.__post_init__ to
    resolve X into phi for display/diagnostics."""
    T_kelvin = T + 273.15
    return _ha_props('R', 'T', T_kelvin, 'P', P, 'W', Y)