from dataclasses import dataclass
import CoolProp.CoolProp as CP

# Functions that return fluid properties based on their state (P, T) and optional humidity

@dataclass
class FluidState:
    rho: float          # density
    cp: float           # specific heat capacity
    lambda_: float      # thermal conductivity
    eta: float          # dynamic viscosity
    Pr: float           # Prandtl number
    R: float = 0.0      # relative humidity

def get_fluid_properties(fluid: str, T_celsius: float, P: float) -> FluidState:
    
    T_kelvin = T_celsius + 273.15

    rho = CP.PropsSI('D', 'P', P, 'T', T_kelvin, fluid)
    cp = CP.PropsSI('C', 'P', P, 'T', T_kelvin, fluid)
    lambda_ = CP.PropsSI('L', 'P', P, 'T', T_kelvin, fluid)
    eta = CP.PropsSI('V', 'P', P, 'T', T_kelvin, fluid)
    Pr = CP.PropsSI('Prandtl', 'P', P, 'T', T_kelvin, fluid)
    
    return FluidState(rho, cp, lambda_, eta, Pr)

def get_humid_air_properties(T_celsius: float, P: float, R: float) -> FluidState:

    T_kelvin = T_celsius + 273.15

    rho = CP.HAPropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
    cp = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
    lambda_ = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
    eta = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)
    Pr = CP.PropsSI('D', 'P', P, 'T', T_kelvin, 'R', R)

    return FluidState(rho, cp, lambda_, eta, Pr)
