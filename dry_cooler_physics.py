from cmath import log10, pi
from CoolProp.CoolProp import PropsSI

# Cooler topologies, flow arrangements, nusselt correlations, fin and tube geometries

def calc_Nu_turbulent(Re: float, Pr: float, d_i: float, L: float):

    psi = (1.8 * log10(Re) - 1.5)**(-2)                                                             # Correction factor for turbulent flow (G1 27)
    Nu_mT = (((psi/8)*(Re - 1000)*Pr)/(1 + 12.7*(psi/8)**(1/2)*(Pr**(2/3) - 1)))*(1+(d_i/L)**(2/3)) # Gnielinski equation for turbulent flow with correction factor (G1 26)

    return Nu_mT

def calc_Nu_laminar(Re: float, Pr: float, d_i: float, L: float):

    Nu_mq1 = 4.364                                                      # Nusselt number for fully developed laminar flow with constant heat flux (G1 17)
    Nu_mq2 = 1.953 * (Re * Pr * (d_i/L))**(1/3)                         # Nusselt number for developing laminar flow with constant heat flux (G1 18)
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i/L))**(1/2)                  # Nusselt number for developing laminar flow with constant heat flux and Prandtl number correction (G1 24)

    Nu_mL = (Nu_mq1**3 + 0.6**3 + (Nu_mq2-0.6)**3 + Nu_mq3**3)**(1/3)   # Combined Nusselt number for laminar flow with constant heat flux (G1 25)

    return Nu_mL

def calc_alpha_i(w: float, d_i: float, L: float, coolant: str):
    rho = PropsSI('D', 'P', 101325, 'T', 293.15, coolant)             # Density of the coolant fluid at ambient conditions [kg/m³]
    eta = PropsSI('V', 'P', 101325, 'T', 293.15, coolant)             # Dynamic viscosity of the coolant fluid at ambient conditions [Pa-s]
    lambda_fluid = PropsSI('L', 'P', 101325, 'T', 293.15, coolant)    # Thermal conductivity of the coolant fluid at ambient conditions [W/m-K]
    c_p = PropsSI('C', 'P', 101325, 'T', 293.15, coolant)             # Specific heat capacity of the coolant fluid at ambient conditions [J/kg-K]
    nu = eta / rho                                                    # Kinematic viscosity [m²/s]
    Pr = c_p * eta / lambda_fluid                                     # Prandtl number

    Re = (w * d_i) / eta                                              # Reynolds number
    Nu                                                                # Nusselt number

    if Re < 2300:                                                                                           # Laminar flow (G1, 3.1)
        Nu = calc_Nu_laminar(Re, Pr, d_i, L)
       
    elif 2300 <= Re < 4000:                                                                                 # Transitional flow (G1, 4.2)
        gamma = (Re - 2300) / (4000 - 2300)                                                                 # Transition parameter (G1 30)
        Nu = (1-gamma) * calc_Nu_laminar(2300, Pr, d_i, L) + gamma * calc_Nu_turbulent(4000, Pr, d_i, L)    # Linear interpolation between laminar and turbulent Nusselt numbers (G1 29)
        
    else:                                                                                                   # turbulent flow (G1 4.1)
        Nu = calc_Nu_turbulent(Re, Pr, d_i, L)

    alpha_i = (Nu * lambda_fluid) / d_i                                                                     # Convective heat transfer coefficient [W/m²-K]

    return alpha_i


def get_geometry():                     # exemplary values for geometry from M1 4
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
    w0_coolant = 0.5  # Cooling velocity [m/s]

    alpha_i = calc_alpha_i(w0_coolant, d_i, coolant)  # Convective heat transfer coefficient at the inner tube wall [W/m²-K]

    n_tubes = 17

    width = n_tubes * t_q                  # Width of the cooler [m]
    height = inflow_cross_section / width  # Height of the cooler [m]

    fins_per_pipe = 348

    A_fin_seg = 2 * (pi/4) * (D**2 - d**2)  # Outer surface area of a segment of the tube between two rips [m²]
    A_pipe = fins_per_pipe * A_fin_seg                 # Outer surface area of the tube per meter length [m²/m]

    A_min = (t_q - d) * a + (t_q - D) * s  # Minimum cross-sectional area for the airflow between the fins and tubes [m²]




    return D, d, s, a, d_i, material, lambda_R, fin_density, t_R, t_q, inflow_cross_section, inflow_velocity, coolant_velocity, alpha_i ##################





