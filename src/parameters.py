from dataclasses import dataclass

@dataclass
class OperatingConditions:
    # Fluids
    coolant_type: str = 'Water'
    air_type: str = 'Air'
    
    # Static Inlet Temperatures [°C]
    T_coolant_in: float = 37.0  # Hot coolant entering  20 - 40 °C
    T_air_in: float = 20.0      # Cold air entering

    
    # Velocities [m/s]
    w_coolant: float = 0.5     # Velocity inside tubes    0.28 kg/s for 15 kw cooling -> 25 °C target output, 37 °C input
    w_o: float = 2.0            # Air approach velocity     output temp vs air velocity
    
    # Pressures [Pa]
    P_coolant: float = 101325
    P_air: float = 101325