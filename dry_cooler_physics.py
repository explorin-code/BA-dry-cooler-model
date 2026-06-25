from dataclasses import dataclass
from math import pi
import CoolProp.CoolProp as CP
import numpy as np

# Cooler topologies, flow arrangements, nusselt correlations, fin and tube geometries

@dataclass
class Geometry:
    
    # 1 Raw Input
    
    D: float                   # Outer diameter of the fin [m]
    d: float                   # Outer diameter of the tube [m]
    s: float                   # Fin thickness [m]
    a: float                   # Fin spacing [m]
    d_i: float                 # Inner diameter of the tube [m]
    material: str              # Fin material
    fins_per_pipe: int         # Number of fins per pipe
    n_tubes: int               # Number of tubes
    n_rows: int                # Number of tube rows
    t_q: float                 # Tube pitch/spacing (center-to-center) [m]
    inflow_cross_section: float # Cross-sectional area for the inflow [m²]
    
    # 2 Intermediate Properties

    @property
    def fin_density(self) -> float:
        """Number of fins per meter."""
        return 9 * 2.54 * 100

    @property
    def t_R(self) -> float:
        """Fin pitch [m] (left edge to next left edge)."""
        return self.a + self.s

    @property
    def width(self) -> float:
        """Width of the cooler [m]. An intermediate variable, safe here!"""
        return self.n_tubes * self.t_q

    @property
    def height(self) -> float:
        """Height of the cooler [m]. Derived from inflow area and width."""
        return self.inflow_cross_section / self.width

    @property
    def lambda_R(self) -> float:
        """Thermal conductivity of the fin material at room temperature [W/m-K]."""
        # CoolProp only does fluids! We use a static lookup for solid metals.
        solid_conductivities = {
            'Aluminum': 237.0,   # W/m-K (Standard pure aluminum)
            'Copper': 401.0,     # W/m-K
            'Carbon Steel': 50.0 # W/m-K
        }
        
        # Look up the material, default to 237.0 if there's a typo
        return solid_conductivities.get(self.material, 237.0)
    
    # 3 Outputs

    @property
    def A_R(self) -> float:
        """Outer surface area of fin segments on one tube [m²]."""
        return 2 * (np.pi / 4) * (self.D**2 - self.d**2) * self.fins_per_pipe

    @property
    def A_G(self) -> float:
        """Exposed base tube area between fins [m²]."""
        return (self.fins_per_pipe + 1) * np.pi * self.d * self.a

    @property
    def A(self) -> float:
        """Total outer surface area product [m²]."""
        return self.A_R + self.A_G

    @property
    def A_i(self) -> float:
        """Inner surface area of one tube [m²]. Uses intermediate height."""
        return self.height * self.d_i * np.pi

    @property
    def Ao_Ae_ratio(self) -> float:
        """Ratio of total to minimum cross-sectional area for airflow [-]."""
        numerator = self.t_q * (self.a + self.s)
        denominator = ((self.t_q - self.d) * self.a) + ((self.t_q - self.D) * self.s)
        return numerator / denominator
    
    @property
    def A_Go(self) -> float:
        """Bare tube surface area per element [m²] (from VDI: pi * d * height)"""
        return np.pi * self.d * self.height
    
    @property
    def l(self) -> float:
        """length of single tube while taking n_row passes through the cooler, ignores the bending radiuses for now"""
        return self.n_rows * self.height

def get_geometry() -> Geometry:
    return Geometry(
        D=0.056,
        d=0.0254,
        s=0.0004,
        a=0.00242,
        d_i=0.0254,
        material='Aluminum',
        fins_per_pipe=348,
        n_tubes=17,
        n_rows=6,
        t_q=0.06,
        inflow_cross_section=1.0,
    )
