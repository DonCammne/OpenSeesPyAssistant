"""Module with the units conversion and the definition of the units used as default (m, N, s). \n
Note that the decision of which unit for each measure (distance, force, mass, time) is equal to 1 is not arbitrary: 
for example the natural frequency is computed behind the scene by the OpenSeesPy framework, thus the stiffness of the structure divided by the mass should result in a unit of 1 (thus seconds). \n
Furthermore, there are constants like the gravitational one g that is dependent on this decision. If the units are used in a consistent way (using this library), these issues can be avoided. \n
Carmine Schipani, 2021
"""


# Fundamental
m_unit = 1.0
length_unit = "m" # It's the length unit associated with 1 (fundamental)
N_unit = 1.0
force_unit = "N" # It's the force unit associated with 1 (fundamental)
s_unit = 1.0
time_unit = "s" # It's the time unit associated with 1 (fundamental)

# Distance
mm_unit = m_unit*1e-3
cm_unit = m_unit*1e-2
dm_unit = m_unit*1e-1
km_unit = m_unit*1e3
inch_unit = m_unit*0.0254
ft_unit = m_unit*0.3048
mile_unit = m_unit*1609.34

# Area
mm2_unit = mm_unit*mm_unit
cm2_unit = cm_unit*cm_unit
dm2_unit = dm_unit*dm_unit
m2_unit = m_unit*m_unit
inch2_unit = inch_unit*inch_unit
ft2_unit = ft_unit*ft_unit

# Volume
mm3_unit = mm_unit*mm_unit*mm_unit
cm3_unit = cm_unit*cm_unit*cm_unit
dm3_unit = dm_unit*dm_unit*dm_unit
m3_unit = m_unit*m_unit*m_unit
inch3_unit = inch_unit*inch_unit*inch_unit
ft3_unit = ft_unit*ft_unit*ft_unit

# Moment of inertia
mm4_unit = mm3_unit*mm_unit
cm4_unit = cm3_unit*cm_unit
dm4_unit = dm3_unit*dm_unit
m4_unit = m3_unit*m_unit
inch4_unit = inch3_unit*inch_unit
ft4_unit = ft3_unit*ft_unit

# Force
kN_unit = N_unit*1e3
MN_unit = N_unit*1e6
GN_unit = N_unit*1e9
kip_unit = N_unit*4448.2216

# Moment (and rotational stiffnes (moment-rotation))
Nm_unit = N_unit*m_unit
kNm_unit = kN_unit*m_unit
MNm_unit = MN_unit*m_unit
Nmm_unit = N_unit*mm_unit
kNmm_unit = kN_unit*mm_unit
MNmm_unit = MN_unit*mm_unit

# Mass
kg_unit = N_unit*s_unit**2/m_unit
t_unit = kg_unit*1e3
pound_unit = kg_unit*0.45359237

# Pressure/Stress
Pa_unit = N_unit/m2_unit
kPa_unit = Pa_unit*1e3
MPa_unit = Pa_unit*1e6
GPa_unit = Pa_unit*1e9
psi_unit = Pa_unit*6894.76
ksi_unit = psi_unit*1000

# Time
min_unit = s_unit*60
hours_unit = min_unit*60

