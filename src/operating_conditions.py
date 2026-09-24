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
from src.param_loader import from_parameters

# Used only when BOTH phi and X come back None -- can't be a field
# default (see get_operating_conditions() call site for why).
DEFAULT_PHI_AIR_FALLBACK = 0.30  # TEMP: placeholder inlet RH until a real site/design value is wired in


@dataclass
class OperatingConditions:
    # --- Raw input -------------------------------------------------------
    T_coolant_in: float            # hot coolant entering [°C]
    T_air_in: float                # cold air entering [°C]

    # coolant side: exactly one nonzero
    u_i: float                      # velocity inside tubes [m/s]
    m_dot_1: float                  # mass flow rate [kg/s]
    V_coolant: float                # volumetric flow rate [m3/s]

    # air side: exactly one nonzero
    w_f: float                      # air approach velocity [m/s]
    m_dot_2: float                  # air mass flow rate [kg/s]
    V_o: float                      # air volumetric flow rate [m3/s]

    # --- Fallbacks/defaults ------------------------------------------------
    coolant_type: str = 'Water'
    P_coolant: float = 101325       # [Pa]
    P_air: float = 101325           # [Pa]
    phi: float = None               # inlet relative humidity [0-1] -- give at most one of phi/X
    X: float = None                 # inlet humidity ratio [kg water/kg dry air]

    # --- Construction-only: reuse an already-built Geometry -----------------
    # Not stored as a field. Pass the SAME Geometry object the caller already
    # built (e.g. main.py's `geo`) instead of triggering a second, independent
    # get_geometry() call inside __post_init__ below.
    geo: InitVar[object] = None

    def __post_init__(self, geo=None):
        # --- Resolve inlet humidity into the canonical X ---------------
        if self.phi is not None and self.X is not None:
            raise ValueError("Give at most one of phi, X -- not both.")
        if self.X is None:
            phi = self.phi if self.phi is not None else DEFAULT_PHI_AIR_FALLBACK
            self.X = relative_to_absolute_humidity(self.T_air_in, self.P_air, phi)
            self.phi = phi
        else:
            self.phi = absolute_to_relative_humidity(self.T_air_in, self.P_air, self.X)

        # --- Validate: exactly one of u_i/m_dot_1/V and w_f/m_dot_2/V given per fluid -------------
        coolant_inputs = (self.u_i, self.m_dot_1, self.V_coolant)
        air_inputs = (self.w_f, self.m_dot_2, self.V_o)

        if sum(1 for x in coolant_inputs if x != 0) != 1:
            raise ValueError("Exactly one of u_i, m_dot_1, V_coolant must be nonzero.")
        if sum(1 for x in air_inputs if x != 0) != 1:
            raise ValueError("Exactly one of w_f, m_dot_2, V_o must be nonzero.")

        # --- Pull the relevant inflow areas from the geometry ------------
        if geo is None:
            geo = get_geometry()
        A_coolant = geo.A_flow_coolant
        A_air = geo.inflow_cross_section

        # --- Pull the relevant inlet densities from fluid_properties -----
        rho_coolant_in = get_fluid_properties(self, self.T_coolant_in, self.P_coolant).rho
        rho_air_in = get_air_properties(self, self.T_air_in, self.P_air).rho

        # --- Coolant side: fill in whichever two were not given ----------
        if self.u_i != 0:
            self.V_coolant = self.u_i * A_coolant
            self.m_dot_1 = rho_coolant_in * self.V_coolant
        elif self.m_dot_1 != 0:
            self.V_coolant = self.m_dot_1 / rho_coolant_in
            self.u_i = self.V_coolant / A_coolant
        else:
            self.m_dot_1 = rho_coolant_in * self.V_coolant
            self.u_i = self.V_coolant / A_coolant

        # --- Air side: fill in whichever two were not given ---------------
        if self.w_f != 0:
            self.V_o = self.w_f * A_air
            self.m_dot_2 = rho_air_in * self.V_o
        elif self.m_dot_2 != 0:
            self.V_o = self.m_dot_2 / rho_air_in
            self.w_f = self.V_o / A_air
        else:
            self.m_dot_2 = rho_air_in * self.V_o
            self.w_f = self.V_o / A_air


def get_operating_conditions(geo=None) -> OperatingConditions:
    """Current operating point -- edit parameters.py to change the numbers.
    Pass geo to resolve flow areas against an already-built Geometry
    instead of constructing a second, independent one internally."""
    return from_parameters(OperatingConditions, {
        'T_coolant_in': 'T_COOLANT_IN',
        'T_air_in': 'T_AIR_IN',
        'u_i': 'W_COOLANT',
        'm_dot_1': 'M_COOLANT',
        'V_coolant': 'V_COOLANT',
        'w_f': 'W_O',
        'm_dot_2': 'M_O',
        'V_o': 'V_O',
        'coolant_type': 'COOLANT_TYPE',
        'P_coolant': 'P_COOLANT',
        'P_air': 'P_AIR',
        'phi': 'PHI_AIR',
        'X': 'X_AIR',
    }, geo=geo)
