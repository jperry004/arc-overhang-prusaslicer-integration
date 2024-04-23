import numpy as np

# Constants and material properties for PLA
E = 3.5e9  # Young's modulus in Pascals (Pa) for PLA at room temperature
rho = 1250  # Density in kg/m^3 for PLA
g = 9.81  # Acceleration due to gravity in m/s^2
speed_xy = 150  # mm/minute, speed in XY plane
layer_thickness = 0.2  # mm
width = 0.6  # mm for extrusion width

# Viscosity for molten PLA
mu = 1000  # Viscosity in Pascal-seconds (Pas)
U = 0.025  # Characteristic velocity in m/s (extrusion speed)

# Cross-sectional area and second moment of area calculations
area_m2_06mm = (width * 1e-3) * (width * 1e-3)
I_06mm = (area_m2_06mm * (width * 1e-3)**2) / 12

# Load per unit length (w) calculation
w_06mm = rho * g * area_m2_06mm

# Beam length for droop calculations (10 mm)
length_m_10mm = 0.01  # Length in meters

# Maximum deflection (δ) calculation for a 10 mm section with 0.6 mm extrusion
delta_10mm_06mm = (5 * w_06mm * length_m_10mm**4) / (384 * E * I_06mm)
delta_10mm_06mm_mm = delta_10mm_06mm * 1e3  # Convert meters to millimeters

# Calculate the maximum length L before significant drooping occurs for molten PLA
L = mu * U / (rho * g)  # Maximum length in meters

# Calculate print time for a 10 mm arc at given radius (approximation for full circle)
arc_length = 2 * np.pi * 10  # Arc length in mm for a 10 mm radius circle
print_time = arc_length / (speed_xy / 60)  # Convert speed to mm/second

import numpy as np

# Constants and properties for PLA
rho = 1250  # Density in kg/m^3 for PLA
g = 9.81  # Acceleration due to gravity in m/s^2
mu = 1000  # Viscosity in Pascal-seconds (Pas)
U = 0.025  # Characteristic velocity in m/s (extrusion speed)
speed_xy = 150  # mm/minute, speed in XY plane

# Kinematic viscosity (ν)
nu = mu / rho

# Calculate droop distance over time
time_seconds = np.linspace(0, 10, num=100)  # Evaluate from 0 to 10 seconds
droop_distance = 0.5 * g * (time_seconds**2) / nu

# Print time for a 10 mm arc at given radius (approximation for full circle)
arc_length = 2 * np.pi * 10  # Arc length in mm for a 10 mm radius circle
print_time_seconds = arc_length / (speed_xy / 60)  # Convert speed to mm/second

# Time to droop 1mm
time_to_droop_1mm = np.sqrt((2 * 1 * nu) / g)

time_to_droop_1mm  # Output the time in seconds

print(f'{time_to_droop_1mm = }')
# Output
# print("Droop distance over time:", droop_distance)


# Output the results
print(f"Maximum length before drooping for molten PLA: {L * 1000:.2f} mm")
print(f"Print time for a 10 mm arc: {print_time:.2f} seconds")


import numpy as np

# Constants and properties for PLA at room temperature (25°C)
rho = 1250  # Density in kg/m^3 for PLA
g = 9.81  # Acceleration due to gravity in m/s^2
mu_180 = 1000  # Viscosity in Pascal-seconds (Pas) for PLA at 180°C
mu_25 = 3000  # Viscosity in Pascal-seconds (Pas) for PLA at room temperature
U = 0.025  # Characteristic velocity in m/s (extrusion speed)
speed_xy = 30  # mm/minute, speed in XY plane

# Linearly interpolate viscosity based on temperature
def linear_interpolate_viscosity(temperature):
    # Linear interpolation formula: mu = mu_25 + (mu_180 - mu_25) * (T - 25) / (180 - 25)
    mu = mu_25 + (mu_180 - mu_25) * (temperature - 25) / (180 - 25)
    return mu

# Input temperature
# temperature = float(input("Enter temperature in Celsius (e.g., 180): "))
temperature = 180
# Calculate viscosity based on linear interpolation
mu = linear_interpolate_viscosity(temperature)

# Calculate kinematic viscosity based on temperature
nu = mu / rho

# Calculate time to droop 1mm
time_to_droop_1mm = np.sqrt((2 * 1 * nu) / g)

# Adjust the ratio of X and Y movement to extrusion rate (alpha)
alpha = 1  # Example: move 2 mm in X and Y for every 1 mm of droop in E

# Calculate the distance traveled in X and Y axes for the same time as droop
distance_xy_mm = time_to_droop_1mm * speed_xy / 60 * alpha  # Convert speed to mm/second

# Output the distance in millimeters
print("Distance traveled in X and Y axes:", distance_xy_mm, "mm")

# Calculate the minimum speed required in XY plane to prevent droop
min_speed_xy = 1 / (time_to_droop_1mm * alpha)  # Convert speed to mm/second

# Output the minimum speed in millimeters per second
print("Minimum speed required in XY plane to prevent droop:", min_speed_xy, "mm/second")


