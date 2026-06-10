from CoolProp.CoolProp import PropsSI

# Functions that return fluid properties based on their state (P, T) and optional humidity

def get_water_properties(P: float, T: float):
    """Return water liquid properties at pressure P [Pa] and temperature T [K].

    Returns density [kg/m³], specific heat [J/kg-K], thermal conductivity [W/m-K], dynamic viscosity [Pa-s].
    """
    rho = PropsSI('D', 'P', P, 'T', T, 'Water')
    cp = PropsSI('C', 'P', P, 'T', T, 'Water')
    k = PropsSI('L', 'P', P, 'T', T, 'Water')
    mu = PropsSI('V', 'P', P, 'T', T, 'Water')
    return rho, cp, k, mu


def get_air_properties(P: float, T: float, RH: float = None):
    """Return air properties at pressure P [Pa] and temperature T [K].

    If RH is None, returns dry air properties.
    If RH is given, returns humid air properties using relative humidity in [0, 1].
    """
    if RH is None:
        fluid = 'Air'
        args = ('P', P, 'T', T, fluid)
    else:
        fluid = 'HumidAir'
        args = ('P', P, 'T', T, 'RH', RH, fluid)

    rho = PropsSI('D', *args)
    cp = PropsSI('C', *args)
    k = PropsSI('L', *args)
    mu = PropsSI('V', *args)
    return rho, cp, k, mu
