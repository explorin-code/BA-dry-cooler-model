from CoolProp.CoolProp import PropsSI

# Cooler topologies, flow arrangements, nusselt correlations, fin and tube geometries

def calc_alpha_i(v: float, d_i: float, coolant: str):
    rho = PropsSI('D', 'P', 101325, 'T', 293.15, coolant)             # Density of the coolant fluid at ambient conditions [kg/m³]
    eta = PropsSI('V', 'P', 101325, 'T', 293.15, coolant)             # Dynamic viscosity of the coolant fluid at ambient conditions [Pa-s]
    lambda_fluid = PropsSI('L', 'P', 101325, 'T', 293.15, coolant)    # Thermal conductivity of the coolant fluid at ambient conditions [W/m-K]
    mu = eta / rho                                                    # Kinematic viscosity [m²/s]

    Re = (v * d_i) / mu                                               # Reynolds number
    Nu = -1                                                           # Nusselt number

    if Re < 2300:                   # Laminar flow (G1, 3.1)
        Nu = 1 ########################
       
    elif 2300 <= Re < 4*10**3:      # Transitional flow (G1, 4.2)
        Nu = 1 ########################
        
    else:                           # turbulent flow (G1 4.1)
        Nu = 3 ########################

    alpha_i = (Nu * lambda_fluid) / d_i  # Convective heat transfer coefficient [W/m²-K]



    return alpha_i


def get_geometry():
    D = 0.056       # Outer diameter of the fin [m]
    d = 0.0254      # Inner diameter of the tube [m]
    s = 0.0004      # Fin thickness [m]
    a = 0.00242     # Fin spacing [m]
    d_i = 0.0254    # Inner diameter of the tube [m]
    
    material = 'Aluminum'  # Fin material
    lambda_R = PropsSI('L', 'T', 293.15, 'P', 101325, material) # Thermal conductivity of the fin material at room temperature [W/m-K]

    fin_density = 9*2.54*100   # Number of fins per meter
    t_R = 0.00282   # fin spacing [m] (from fin left edge to left edge of next fin)
    t_q = 0.06      # tube spacing [m] (from tube center to tube center)
    inflow_cross_section = 1 # Cross-sectional area for the inflow [m²]

    inflow_velocity = 2.0  # Inflow velocity [m/s]
    
    coolant = 'Water'  # Coolant fluid
    coolant_velocity = 0.5  # Cooling velocity [m/s]

    alpha_i = calc_alpha_i(coolant_velocity, d_i, coolant)  # Convective heat transfer coefficient at the inner tube wall [W/m²-K]





    return D, d, s, a, d_i, material, lambda_R, fin_density, t_R, t_q, inflow_cross_section, inflow_velocity, coolant_velocity, alpha_i





