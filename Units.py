# Module for units conversion to the default ones (m, N, kg, s)
#   Carmine Schipani, 2021

"""Module with the units conversion and the definition of the units used as default (m, N, kg, s)
"""
m_unit = 1
N_unit = 1
kg_unit = 1
s_unit = 1

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

# Moment
Nm_unit = N_unit*m_unit
kNm_unit = kN_unit*m_unit
MNm_unit = MN_unit*m_unit
Nmm_unit = N_unit*mm_unit
kNmm_unit = kN_unit*mm_unit
MNmm_unit = MN_unit*mm_unit

# Mass
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

#TODO: check with Andy

