import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

# 1. Define our constant pressure (1 bar = 100,000 Pa)
pressure = 100000 

# 2. Create an array of temperatures from 5°C to 95°C
# We calculate in Kelvin (add 273.15) because CoolProp requires SI units
temperatures_C = np.linspace(5, 95, 100)
temperatures_K = temperatures_C + 273.15

# 3. Calculate water density for every single temperature in our array
densities = [CP.PropsSI('D', 'P', pressure, 'T', T, 'Water') for T in temperatures_K]

# 4. Create the plot
plt.figure(figsize=(8, 5))
plt.plot(temperatures_C, densities, color='blue', linewidth=2, label='Water at 1 bar')

# 5. Customize the chart (labels, title, grid)
plt.title('Water Density vs. Temperature (CoolProp Data)', fontsize=14, fontweight='bold')
plt.xlabel('Temperature (°C)', fontsize=12)
plt.ylabel('Density (kg/m³)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# 6. Show the plot on your screen
plt.show()