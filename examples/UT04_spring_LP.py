##################################################################################################################################
#                                 						  Example - UT04
#                                                Loading Protocol - Spring-Based Members
#                                                                                                                                
# Date Created: January 7, 2022
# 	Carmine Schipani
#                       
##################################################################################################################################

# Import libraries
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import math
from copy import copy, deepcopy
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.Fibers import *
from OpenSeesPyAssistant.Connections import *
from OpenSeesPyAssistant.FunctionalFeatures import *
from OpenSeesPyAssistant.MemberModel import *
from OpenSeesPyAssistant.AnalysisAndPostProcessing import *
from OpenSeesPyAssistant.GeometryTemplate import *

####################################################################################################
#             								INITIALIZATION						                   #
####################################################################################################
# Set up
data_dir = "UT04_LP"
Initialize2DModel(data_dir)

# Define geometric transformation
geo_transf_ID = 1
geomTransf('PDelta', geo_transf_ID)

# Show info
show_info = True
show_plots = True

####################################################################################################
#                                        INPUT VARIABLE                                            #
####################################################################################################
### Member parameters and sections
HalfBayWidth = 178.188*inch_unit
HalfHStory1 = 109*inch_unit
HalfHStory2 = 92.625*inch_unit
AxisY3 = HalfHStory1+HalfHStory2

# NB: the depth must remain continuous between same type of elements (columns and beams)
d_Col = 18.25*inch_unit
d_Beam = 35.875*inch_unit
r_Col = (4.125-2.85)*inch_unit
r_Beam = (1.875-0.94)*inch_unit
E = 200.0*GPa_unit

# Panel Zone
Fy_pz = 337.8*MPa_unit;     # Fy wed of the column
pzwidth = d_Col/2           # Half of the PZ width
pzheight = d_Beam/2         # Half of the PZ height

# Bottom Column Section Properties (W14x398)
col_bottom_section = SteelIShape("Col", d_Col, 16.5*inch_unit, 2.86*inch_unit, 1.784*inch_unit, HalfHStory1-pzheight, r_Beam, E, 351.6*MPa_unit, Fy_pz, "Bottom Column")
if show_info: col_bottom_section.ShowInfo()

# Bottom Column Section Properties (W14x398)
col_top_section = SteelIShape("Col", d_Col, 16.5*inch_unit, 2.86*inch_unit, 1.784*inch_unit, HalfHStory2-pzheight, r_Beam, E, 351.6*MPa_unit, Fy_pz, "Top Column")
if show_info: col_top_section.ShowInfo()

# East Beam Section Properties (W36x150) (right)
beam_east_section = SteelIShape("Beam", d_Beam, 11.875*inch_unit, 0.892*inch_unit, 0.650*inch_unit, HalfBayWidth-pzwidth, r_Beam, E, 358.5*MPa_unit, name_tag="East Beam")
if show_info: beam_east_section.ShowInfo()

# West Beam Section Properties (W36x150) (left)
beam_west_section = SteelIShape("Beam", d_Beam, 11.875*inch_unit, 0.890*inch_unit, 0.643*inch_unit, HalfBayWidth-pzwidth, r_Beam, E, 365.4*MPa_unit, name_tag="West Beam")
if show_info: beam_west_section.ShowInfo()


###################################################################################################
#                                    ACTIONS (FORCES, DISPLACEMENT, LP)                           #
###################################################################################################
# Forces

# Loading protocol
SDR_LP = np.array([0.00375, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06])
nr_cycles_LP = np.array([6, 6, 6, 4, 2, 2, 2, 2, 2, 2])
discr = DiscretizeLoadProtocol(SDR_LP, nr_cycles_LP, 21, show_plots, True)
discr_LP = discr * AxisY3


####################################################################################################
#                                             MODEL                                                #
####################################################################################################
# Geometry (grid nodes)
DefineSubassemblageNodes(beam_west_section.L, beam_east_section.L, col_top_section.L, col_bottom_section.L, d_Col, d_Beam, show_plot=show_plots)
ctrl_node = 23

# Panel zone
panel_zone = PanelZoneSteelIShapeSkiadopoulos2021(22, col_top_section, beam_west_section, geo_transf_ID)
panel_zone.CreateMember()
if show_info: panel_zone.ShowInfo(show_plots, False)

# Beam left
beam_west_spring = ModifiedIMKSteelIShape(220832208, beam_west_section)
beam_west_spring.Bilin()
if show_info: beam_west_spring.ShowInfo(show_plots, False)
beam_west = SpringBasedElementSteelIShape(12, 2208, beam_west_section, geo_transf_ID, mat_ID_j=beam_west_spring.ID)
# beam_west = SpringBasedElementModifiedIMKSteelIShape(12, 2208, beam_west_section, geo_transf_ID, new_mat_ID_j=220832208)
beam_west.CreateMember()
if show_info: beam_west.ShowInfo(show_plots, False)

# Beam right
beam_east_spring = ModifiedIMKSteelIShape(220222022, beam_east_section)
beam_east_spring.Bilin()
if show_info: beam_east_spring.ShowInfo(show_plots, False)
beam_east = SpringBasedElementSteelIShape(2202, 32, beam_east_section, geo_transf_ID, mat_ID_i=beam_east_spring.ID)
# beam_east = SpringBasedElementModifiedIMKSteelIShape(2202, 32, beam_east_section, geo_transf_ID, new_mat_ID_i=220222022)
beam_east.CreateMember()
if show_info: beam_east.ShowInfo(show_plots, False)

# Column top
col_top_spring = ModifiedIMKSteelIShape(22622, col_top_section)
col_top_spring.Bilin()
if show_info: col_top_spring.ShowInfo(show_plots, False)
col_top = SpringBasedElementSteelIShape(22, 23, col_top_section, geo_transf_ID, mat_ID_i=col_top_spring.ID)
# col_top = SpringBasedElementModifiedIMKSteelIShape(22, 23, col_top_section, geo_transf_ID, new_mat_ID_i=22622)
col_top.CreateMember()
if show_info: col_top.ShowInfo(show_plots, False)

# Column bottom
col_bottom_spring = ModifiedIMKSteelIShape(220522057, col_bottom_section)
col_bottom_spring.Bilin()
if show_info: col_bottom_spring.ShowInfo(show_plots, False)
col_bottom = SpringBasedElementSteelIShape(21, 2205, col_bottom_section, geo_transf_ID, mat_ID_j=col_bottom_spring.ID)
# col_bottom = SpringBasedElementModifiedIMKSteelIShape(21, 2205, col_bottom_section, geo_transf_ID, new_mat_ID_j=220522057)
col_bottom.CreateMember()
if show_info: col_bottom.ShowInfo(show_plots, False)

# Model complete
print("")
print("Model Built")

# Show plots (if any)
plt.show(block = False)
plt.pause(0.1)

# Display structure with node and element IDs
opsplt.plot_model("nodes", "elements")


###################################################################################################
#                                             RECORDERS                                           #
###################################################################################################
# Roof lateral displacement
recorder("Node", "-file", '{}/LatDispl.txt'.format(data_dir), "-time", "-node", ctrl_node, "-dof", 1, "disp")

# Force
col_top.Record("element", "TipColumnForce", data_dir, True, False, True)

# Deformation
panel_zone.Record("PZDef", data_dir, False, True, False)
col_top.Record("spring_i", "ColSpringT", data_dir, False, True, False)
col_bottom.Record("spring_j", "ColSpringB", data_dir, False, True, False)
beam_west.Record("spring_j", "BeamSpringW", data_dir, False, True, False)
beam_east.Record("spring_i", "BeamSpringE", data_dir, False, True, False)


###################################################################################################
#                                    ANALYSIS AND POSTPROCESSING                                  #
###################################################################################################
# Initialize analysis
analysis_opt = Analysis(data_dir, data_dir)

# Gravity analysis
# analysis_opt.Gravity([ctrl_node], [-N_G], 1, 1, show_plot=show_plots)

# Pushover analysis
# analysis_opt.LateralForce([ctrl_node], [1*kN_unit], 2, 2, show_plot=show_plots, block = True)
# analysis_opt.Pushover(ctrl_node, (HalfHStory1+HalfHStory2)*0.02, 0.04*mm_unit, 2, 2,
    # show_plot=show_plots)
analysis_opt.LoadingProtocol(ctrl_node, discr_LP, 2, 2,
    show_plot=show_plots)

# Postprocessing
analysis_opt.DeformedShape(animate = True)

# Save info of sections, elements, material models,...
print("Saving infos for postprocessing analysis")
with open("{}/SavedInfos.txt".format(data_dir), 'w') as f:
    col_bottom_section.SaveData(f)
    col_top_section.SaveData(f)
    beam_east_section.SaveData(f)
    beam_west_section.SaveData(f)
    col_bottom_spring.SaveData(f)
    col_top_spring.SaveData(f)
    beam_east_spring.SaveData(f)
    beam_west_spring.SaveData(f)
    panel_zone.SaveData(f)

