from CoolProp.CoolProp import PropsSI

# Functions that return fluid properties based on their state (P, T) and optional humidity

def get_coolant_properties(P: float, T: float, coolant: str):
   
    rho = PropsSI('D', 'P', P, 'T', T, coolant)             # Density of coolant at given pressure and temperature [kg/m³]
    c_p = PropsSI('C', 'P', P, 'T', T, coolant)             # Heat capacity of coolant at given pressure and temperature [J/kg-K]
    lambda_coolant = PropsSI('L', 'P', P, 'T', T, coolant)  # Thermal conductivity of coolant at given pressure and temperature [W/m-K]
    eta = PropsSI('V', 'P', P, 'T', T, coolant)             # Dynamic viscosity of coolant at given pressure and temperature [Pa-s]
    Pr_fluid = c_p * eta / lambda_coolant                   # Prandtl number of coolant [-]

    return rho, c_p, lambda_coolant, eta, Pr_fluid


def get_air_properties(P: float, T: float, RH: float = None):

    if RH is None:
        fluid = 'Air'
        args = ('P', P, 'T', T, fluid)
    else:
        fluid = 'HumidAir'
        args = ('P', P, 'T', T, 'RH', RH, fluid)

    rho = PropsSI('D', *args)                               # Density of air or humid air at given pressure, temperature, and optional relative humidity [kg/m³]
    c_p = PropsSI('C', *args)                               # Heat capacity of air or humid air at given pressure, temperature, and optional relative humidity [J/kg-K]                 
    lambda_air = PropsSI('L', *args)                        # Thermal conductivity of air or humid air at given pressure, temperature, and optional relative humidity [W/m-K]
    eta = PropsSI('V', *args)                               # Dynamic viscosity of air or humid air at given pressure, temperature, and optional relative humidity [Pa-s]
    Pr_air = c_p * eta / lambda_air                         # Prandtl number of air or humid air [-]

    return rho, c_p, lambda_air, eta, Pr_air
