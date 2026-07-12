from dataclasses import dataclass

@dataclass
class OperatingConditions:
    # Fluids
    coolant_type: str = 'Water'
    air_type: str = 'Air'
    
    # Static Inlet Temperatures [°C]
    T_coolant_in: float = 90.0  # Hot coolant entering
    T_air_in: float = 20.0      # Cold air entering
    
    # Velocities [m/s]
    w_coolant: float = 0.5      # Velocity inside tubes
    w_o: float = 2.0            # Air approach velocity
    
    # Pressures [Pa]
    P_coolant: float = 101325
    P_air: float = 101325