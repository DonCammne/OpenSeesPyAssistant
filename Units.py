# Module for units conversion to the default ones (m, N, kg, s)
#   Carmine Schipani, 2021

"""Module with the units conversion and the definition of the units used as default (m, N, kg, s)
"""
m = 1
N = 1
kg = 1
s = 1

# Distance
mm = m*1e-3
cm = m*1e-2
dm = m*1e-1
km = m*1e3
inch = m*0.0254
ft = m*0.3048
mile = m*1609.34

# Area
mm2 = mm*mm
cm2 = cm*cm
dm2 = dm*dm
m2 = m*m
inch2 = inch*inch
ft2 = ft*ft

# Volume
mm3 = mm*mm*mm
cm3 = cm*cm*cm
dm3 = dm*dm*dm
m3 = m*m*m
inch3 = inch*inch*inch
ft3 = ft*ft*ft

# Moment of inertia
mm4 = mm3*mm
cm4 = cm3*cm
dm4 = dm3*dm
m4 = m3*m
inch4 = inch3*inch
ft4 = ft3*ft

# Force
kN = N*1e3
MN = N*1e6
GN = N*1e9
kip = N*4448.2216

# Mass
t = kg*1e3
pound = kg*0.45359237

# Pressure/Stress
Pa = N/m2
kPa = Pa*1e3
MPa = Pa*1e6
GPa = Pa*1e9
psi = Pa*6894.76
ksi = psi*1000

# Time
min = s*60
hours = min*60


