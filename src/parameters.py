"""
parameters.py
==============
Operating conditions (inlet temperatures, flow rates, pressures) for a
solver run. For each fluid, specify exactly ONE of velocity / mass flow
/ volumetric flow -- the other two are derived automatically in
__post_init__ using the geometry's flow area and the inlet density.
"""

from dataclasses import dataclass

from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties


@dataclass
class OperatingConditions:
    # --- Fluids --------------------------------------------------------
    coolant_type: str = 'Water'
    air_type: str = 'Air'

    # --- Static inlet temperatures [°C] ---------------------------------
    T_coolant_in: float = 37.0     # hot coolant entering  -- Konrad: 37 °C, 25 °C target output
    T_air_in: float = 15.0         # cold air entering

    # --- Coolant side: specify EXACTLY ONE of the following three ------
    # (the other two are derived automatically from geometry + inlet density)
    w_coolant: float = 0.0         # velocity inside tubes         [m/s]
    m_coolant: float = 0.28        # mass flow rate                [kg/s]  -- Konrad: 0.28 kg/s
    V_coolant: float = 0.0         # volumetric flow rate          [m3/s]

    # --- Air side: specify EXACTLY ONE of the following three ----------
    w_o: float = 2.0                # air approach velocity        [m/s]  -- sweet spot ~2-2.5 m/s
    m_o: float = 0.0                # air mass flow rate           [kg/s]
    V_o: float = 0.0                # air volumetric flow rate     [m3/s]

    # --- Pressures [Pa] --------------------------------------------------
    P_coolant: float = 101325
    P_air: float = 101325

    def __post_init__(self):
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
        rho_coolant_in = get_fluid_properties(self.coolant_type, self.T_coolant_in, self.P_coolant).rho
        rho_air_in = get_fluid_properties(self.air_type, self.T_air_in, self.P_air).rho

        # --- Coolant side: fill in whichever two were not given ----------
        if self.w_coolant != 0:
            self.V_coolant = self.w_coolant * A_coolant
            self.m_coolant = rho_coolant_in * self.V_coolant
        elif self.m_coolant != 0:
            self.V_coolant = self.m_coolant / rho_coolant_in
            self.w_coolant = self.V_coolant / A_coolant
        else:  # V_coolant given
            self.m_coolant = rho_coolant_in * self.V_coolant
            self.w_coolant = self.V_coolant / A_coolant

        # --- Air side: fill in whichever two were not given ---------------
        if self.w_o != 0:
            self.V_o = self.w_o * A_air
            self.m_o = rho_air_in * self.V_o
        elif self.m_o != 0:
            self.V_o = self.m_o / rho_air_in
            self.w_o = self.V_o / A_air
        else:  # V_o given
            self.m_o = rho_air_in * self.V_o
            self.w_o = self.V_o / A_air