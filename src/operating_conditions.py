"""
operating_conditions.py
========================
Operating conditions for a solver run (temperatures, flows, pressures,
humidity) -- derives whichever inputs weren't explicitly given. Raw
values come from parameters.py.
"""

from dataclasses import dataclass

from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties, relative_to_absolute_humidity, absolute_to_relative_humidity
from src.param_loader import from_parameters

# Used only when BOTH phi_air and X_air come back None -- can't be a field
# default (see get_operating_conditions() call site for why).
DEFAULT_PHI_AIR_FALLBACK = 0.30  # TEMP: placeholder inlet RH until a real site/design value is wired in


@dataclass
class OperatingConditions:
    # --- Raw input -------------------------------------------------------
    T_coolant_in: float            # hot coolant entering [°C]
    T_air_in: float                # cold air entering [°C]

    # coolant side: exactly one nonzero
    w_coolant: float                # velocity inside tubes [m/s]
    m_coolant: float                # mass flow rate [kg/s]
    V_coolant: float                # volumetric flow rate [m3/s]

    # air side: exactly one nonzero
    w_o: float                      # air approach velocity [m/s]
    m_o: float                      # air mass flow rate [kg/s]
    V_o: float                      # air volumetric flow rate [m3/s]

    # --- Fallbacks/defaults ------------------------------------------------
    coolant_type: str = 'Water'
    P_coolant: float = 101325       # [Pa]
    P_air: float = 101325           # [Pa]
    phi_air: float = None           # inlet relative humidity [0-1] -- give at most one of phi_air/X_air
    X_air: float = None             # inlet humidity ratio [kg water/kg dry air]

    def __post_init__(self):
        # --- Resolve inlet humidity into the canonical X_air ---------------
        if self.phi_air is not None and self.X_air is not None:
            raise ValueError("Give at most one of phi_air, X_air -- not both.")
        if self.X_air is None:
            phi = self.phi_air if self.phi_air is not None else DEFAULT_PHI_AIR_FALLBACK
            self.X_air = relative_to_absolute_humidity(self.T_air_in, self.P_air, phi)
            self.phi_air = phi
        else:
            self.phi_air = absolute_to_relative_humidity(self.T_air_in, self.P_air, self.X_air)

        # --- Validate: exactly one of w/m/V given per fluid -------------
        coolant_inputs = (self.w_coolant, self.m_coolant, self.V_coolant)
        air_inputs = (self.w_o, self.m_o, self.V_o)

        if sum(1 for x in coolant_inputs if x != 0) != 1:
            raise ValueError("Exactly one of w_coolant, m_coolant, V_coolant must be nonzero.")
        if sum(1 for x in air_inputs if x != 0) != 1:
            raise ValueError("Exactly one of w_o, m_o, V_o must be nonzero.")

        # --- Pull the relevant inflow areas from the geometry ------------
        geo = get_geometry()
        A_coolant = geo.A_flow_coolant
        A_air = geo.inflow_cross_section

        # --- Pull the relevant inlet densities from fluid_properties -----
        rho_coolant_in = get_fluid_properties(self, self.T_coolant_in, self.P_coolant).rho
        rho_air_in = get_air_properties(self, self.T_air_in, self.P_air).rho

        # --- Coolant side: fill in whichever two were not given ----------
        if self.w_coolant != 0:
            self.V_coolant = self.w_coolant * A_coolant
            self.m_coolant = rho_coolant_in * self.V_coolant
        elif self.m_coolant != 0:
            self.V_coolant = self.m_coolant / rho_coolant_in
            self.w_coolant = self.V_coolant / A_coolant
        else:
            self.m_coolant = rho_coolant_in * self.V_coolant
            self.w_coolant = self.V_coolant / A_coolant

        # --- Air side: fill in whichever two were not given ---------------
        if self.w_o != 0:
            self.V_o = self.w_o * A_air
            self.m_o = rho_air_in * self.V_o
        elif self.m_o != 0:
            self.V_o = self.m_o / rho_air_in
            self.w_o = self.V_o / A_air
        else:
            self.m_o = rho_air_in * self.V_o
            self.w_o = self.V_o / A_air


def get_operating_conditions() -> OperatingConditions:
    """Current operating point -- edit parameters.py to change the numbers."""
    return from_parameters(OperatingConditions, {
        'T_coolant_in': 'T_COOLANT_IN',
        'T_air_in': 'T_AIR_IN',
        'w_coolant': 'W_COOLANT',
        'm_coolant': 'M_COOLANT',
        'V_coolant': 'V_COOLANT',
        'w_o': 'W_O',
        'm_o': 'M_O',
        'V_o': 'V_O',
        'coolant_type': 'COOLANT_TYPE',
        'P_coolant': 'P_COOLANT',
        'P_air': 'P_AIR',
        'phi_air': 'PHI_AIR',
        'X_air': 'X_AIR',
    })