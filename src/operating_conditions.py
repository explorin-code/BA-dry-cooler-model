"""
operating_conditions.py
========================
Operating conditions for a solver run (temperatures, flows, pressures,
humidity) -- derives whichever inputs weren't explicitly given. Raw
values come from parameters.py.
"""

from dataclasses import dataclass, InitVar

from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties, relative_to_absolute_humidity, absolute_to_relative_humidity
from src.param_loader import from_parameters, PERCENT

# Used only when BOTH phi and X come back None -- can't be a field
# default (see get_operating_conditions() call site for why).
DEFAULT_PHI_AIR_FALLBACK = 0.30  # Assumption: inlet relative humidity if neither PHI_AIR_PERCENT nor Y_AIR is given (fraction)


@dataclass
class OperatingConditions:
    # --- Raw input -------------------------------------------------------
    T_c_i: float            # hot coolant entering [°C]
    T_a_i: float                # cold air entering [°C]

    # coolant side: exactly one nonzero
    w_c: float                      # velocity inside tubes [m/s]
    m_dot_c: float                  # mass flow rate [kg/s]
    V_dot_c: float                # volumetric flow rate [m3/s]

    # air side: exactly one nonzero
    w_fr: float                      # air approach velocity [m/s]
    m_dot_a: float                  # air mass flow rate [kg/s]
    V_dot_a: float                      # air volumetric flow rate [m3/s]

    # --- Fallbacks/defaults ------------------------------------------------
    coolant_type: str = 'Water'
    p_c: float = 101325       # [Pa]
    p_a: float = 101325           # [Pa]
    phi: float = None               # inlet relative humidity [0-1] (input in %) -- give at most one of phi/Y
    Y: float = None                 # inlet humidity ratio [kg water/kg dry air]

    # --- Construction-only: reuse an already-built Geometry -----------------
    # Not stored as a field. Pass the SAME Geometry object the caller already
    # built (e.g. main.py's `geo`) instead of triggering a second, independent
    # get_geometry() call inside __post_init__ below.
    geo: InitVar[object] = None

    def __post_init__(self, geo=None):
        # --- Resolve inlet humidity into the canonical X ---------------
        if self.phi is not None and self.Y is not None:
            raise ValueError("Give at most one of PHI_AIR_PERCENT / Y_AIR (phi, Y) -- not both.")
        if self.Y is None:
            phi = self.phi if self.phi is not None else DEFAULT_PHI_AIR_FALLBACK
            self.Y = relative_to_absolute_humidity(self.T_a_i, self.p_a, phi)
            self.phi = phi
        else:
            self.phi = absolute_to_relative_humidity(self.T_a_i, self.p_a, self.Y)

        # --- Validate: exactly one of w_c/m_dot_c/V and w_fr/m_dot_a/V given per fluid -------------
        coolant_inputs = (self.w_c, self.m_dot_c, self.V_dot_c)
        air_inputs = (self.w_fr, self.m_dot_a, self.V_dot_a)

        if sum(1 for x in coolant_inputs if x != 0) != 1:
            raise ValueError("Exactly one of W_COOLANT / M_COOLANT / V_COOLANT (w_c, m_dot_c, V_dot_c) must be nonzero.")
        if sum(1 for x in air_inputs if x != 0) != 1:
            raise ValueError("Exactly one of W_O / M_O / V_O (w_fr, m_dot_a, V_dot_a) must be nonzero.")

        # --- Pull the relevant inflow areas from the geometry ------------
        if geo is None:
            geo = get_geometry()
        A_coolant = geo.A_cs_c
        A_air = geo.A_fr

        # --- Pull the relevant inlet densities from fluid_properties -----
        rho_c_i = get_fluid_properties(self, self.T_c_i, self.p_c).rho
        rho_a_i = get_air_properties(self, self.T_a_i, self.p_a).rho

        # --- Coolant side: fill in whichever two were not given ----------
        if self.w_c != 0:
            self.V_dot_c = self.w_c * A_coolant
            self.m_dot_c = rho_c_i * self.V_dot_c
        elif self.m_dot_c != 0:
            self.V_dot_c = self.m_dot_c / rho_c_i
            self.w_c = self.V_dot_c / A_coolant
        else:
            self.m_dot_c = rho_c_i * self.V_dot_c
            self.w_c = self.V_dot_c / A_coolant

        # --- Air side: fill in whichever two were not given ---------------
        if self.w_fr != 0:
            self.V_dot_a = self.w_fr * A_air
            self.m_dot_a = rho_a_i * self.V_dot_a
        elif self.m_dot_a != 0:
            self.V_dot_a = self.m_dot_a / rho_a_i
            self.w_fr = self.V_dot_a / A_air
        else:
            self.m_dot_a = rho_a_i * self.V_dot_a
            self.w_fr = self.V_dot_a / A_air


def get_operating_conditions(geo=None, **overrides) -> OperatingConditions:
    """Current operating point -- edit parameters.py to change the numbers.
    Pass geo to resolve flow areas against an already-built Geometry
    instead of constructing a second, independent one internally.
    overrides (e.g. T_a_i=..., phi=...) replace single parameters.py values
    for this call only (the annual benchmark's hourly weather)."""
    return from_parameters(OperatingConditions, {
        'T_c_i': 'T_COOLANT_IN',
        'T_a_i': 'T_AIR_IN',
        'w_c': 'W_COOLANT',
        'm_dot_c': 'M_COOLANT',
        'V_dot_c': 'V_COOLANT',
        'w_fr': 'W_O',
        'm_dot_a': 'M_O',
        'V_dot_a': 'V_O',
        'coolant_type': 'COOLANT_TYPE',
        'p_c': 'P_COOLANT',
        'p_a': 'P_AIR',
        'phi': ('PHI_AIR_PERCENT', PERCENT),
        'Y': 'Y_AIR',
    }, geo=geo, **overrides)
