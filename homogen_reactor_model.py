###########################################################################################################################################################
###########################################################################################################################################################
############################################# File for Generating Materials, Geometry, and Settings xml Files #############################################
#############################################              Also runs the eigenvalue calculations              #############################################
###########################################################################################################################################################
###########################################################################################################################################################

###########################################################################################################################################################
#################################################################### Relevant Libraries ###################################################################
###########################################################################################################################################################
import numpy as np
import openmc
import openmc.deplete
import argparse

###########################################################################################################################################################
############################################################         Reactor Parameters         ###########################################################
############################################################ Density and Temperature Parameters ###########################################################
############################################################    Temp. in K, legnths in cm      ############################################################
###########################################################################################################################################################
N_blocks = 4 # number of reactor blocks in one vertical column
N_flooded = 0 # number of reactor blocks where the coolant channels are filled with water
h_reflector = 20 # height of the top/bottom reflector, 20 cm for BeO, 50 cm for graphite
###########################################################################################################################################################
'''
Temperatures may have to be changed several times for the flooding scenario. Just comment out the default ones and write in new values
Temperature values (except T_FPB, T_coolant_average, T_SS, and T_water) were adapted from:
    "Variable Reactivity Control in Small Modular High Temperature Reactors Using Moderation Manipulation Techniques" paper by Atkinson et al (2018)
Other temperature values:
    T_coolant_average from:
        "Small modular high temperature reactor optimisation part 2: Reactivity control for prismatic core high temperature small modular reactor, including
         fixed burnable poisons, spectrum hardening and control rods" by Atkinson et al (2018)
    T_FBP approximated as being the same as the graphite block material
'''


""" Original temp set, with explenation found in the summer project report. 
T_fuel = 1023.15
T_FBP = 973.5
T_graphite = 973.5
T_coolant_average = 674.15
T_helium = 600
T_reflector = 873.15
T_insulation = 973.5
T_barrel = 673.5
T_RPV = 300
T_water = 273 # temperature for the water flooding the reactor
"""



# Temp set for Normal Operation case (chatGPT suggestion)
T_fuel            = 1373.15   # 1100 C
T_graphite        = 1323.15   # 1050 C
T_FBP             = 1323.15   # 1050 C

T_coolant_average = 1123.15   # 850 C
T_helium          = 1123.15   # 850 C

T_reflector       = 1023.15   # 750 C
T_insulation      = 873.15    # 600 C
T_barrel          = 773.15    # 500 C
T_RPV             = 573.15    # 300 C

# accident/flooding case:
T_water           = 273.15    # 0 C


"""
# Temp set for DBA case (chatGPT suggestion)
T_fuel            = 1873.15   # 1600 C
T_graphite        = 1623.15   # 1323 + 300
T_FBP             = 1623.15

T_coolant_average = 1323.15   # 1123 + 200
T_helium          = 1323.15

T_reflector       = 1173.15   # 1023 + 150
T_insulation      = 973.15    # 873 + 100
T_barrel          = 848.15    # 773 + 75
T_RPV             = 623.15    # 573 + 50
T_water           = 273.15    # 0 C
"""
"""
# Temp set for BDBA case (chatGPT suggestion)
T_fuel            = 2073.15   # 1800 C
T_graphite        = 1743.15   # 1323 + 420
T_FBP             = 1743.15

T_coolant_average = 1403.15   # 1123 + 280
T_helium          = 1403.15

T_reflector       = 1233.15   # 1023 + 210
T_insulation      = 1013.15   # 873 + 140
T_barrel          = 878.15    # 773 + 105
T_RPV             = 643.15    # 573 + 70
T_water           = 273.15    # 0 C
"""

###########################################################################################################################################################
P_coolant_averarge = 4 # in MPa 

###########################################################################################################################################################
################################################################ TRISO particle parameters ################################################################
################################################################        lengths in cm      ################################################################
###########################################################################################################################################################
enrichment = 0.20
r_fuel_kernel = 0.025
r_buffer = r_fuel_kernel + 100e-4
r_IPyC = r_buffer + 35e-4
r_SiC = r_IPyC + 35e-4
r_OPyC = r_SiC + 40e-4
r_array = [r_fuel_kernel, r_buffer, r_IPyC, r_SiC, r_OPyC]

###########################################################################################################################################################
################################################################# Fuel Compact parameters #################################################################
#################################################################      lengths in cm      #################################################################
###########################################################################################################################################################
packing_fraction = 0.37
r_fuel_compact = 0.6225
r_fuel_channel = 0.6350
t_fuel_gap = r_fuel_channel - r_fuel_compact
h_fuel_compact = 4.93

###########################################################################################################################################################
############################################################### Homogeneous fuel parameters ###############################################################
#################################################################      lengths in cm      #################################################################
###########################################################################################################################################################
'''
For details see the fuel_testing.py file. And fuel homogenization text in report
'''
r_fuel_RPT = 0.4819
r_fuel_GR = r_fuel_channel - t_fuel_gap # outer radius of the graphite ring around the homoegeneous fuel
###########################################################################################################################################################
N_particles = 4416
V_particles = N_particles * 4/3 * np.pi * r_OPyC**3
V_fuel = N_particles * 4/3 * np.pi * r_fuel_kernel**3
V_buffer = N_particles * 4/3 * np.pi * (r_buffer**3 - r_fuel_kernel**3)
V_IPyC = N_particles * 4/3 * np.pi * (r_IPyC**3 - r_buffer**3)
V_SiC = N_particles * 4/3 * np.pi * (r_SiC**3 - r_IPyC**3)
V_OPyC = N_particles * 4/3 * np.pi * (r_OPyC**3 - r_SiC**3)
V_compact = h_fuel_compact * np.pi * r_fuel_RPT**2 - V_particles
V_total = h_fuel_compact * np.pi * r_fuel_RPT**2
V_ratios = np.array([V_fuel, V_buffer, V_IPyC, V_SiC, V_OPyC, V_compact])/V_total

###########################################################################################################################################################
################################################################# Fuel Assembly parameters ################################################################
#################################################################      lengths in cm       ################################################################
###########################################################################################################################################################
pitch_assembly = 36.10
d_flat_to_flat = 36
h_block = 80
h_plug = 3.025 # height of the graphite plugs at the top and bottom of the fuel channel
pitch_holes = 1.88 # pitch of the fuel/FBP/coolant channels in the block

###########################################################################################################################################################
################################################################ Burnable poison parameters ###############################################################
#################################################################      lengths in cm       ################################################################
###########################################################################################################################################################
FBP_packing_fraction = 0.14
r_FBP_kernel = 0.01
r_FBP_buffer = r_FBP_kernel + 0.0018
r_FBP_PyC = r_FBP_buffer + 0.0023
###########################################################################################################################################################
h_FBP_channel = 72.14
N_FBP = 10 # axial divisions for FBP channel, arbitrary
h_FBP_compact = h_FBP_channel/N_FBP # just to have some axial divisions for better calculations
r_FBP_channel = 0.6350
r_FBP_compact = 1.143/2
t_FBP_gap = r_FBP_channel - r_FBP_compact
###########################################################################################################################################################
r_FBP_RPT = 0.5020
r_FBP_GR = r_FBP_channel - t_FBP_gap
###########################################################################################################################################################
N_particles_FBP = 48935
V_particles_FBP = N_particles_FBP * 4/3 * np.pi * r_FBP_PyC**3
V_B4C = N_particles_FBP * 4/3 * np.pi * r_FBP_kernel**3
V_buffer_FBP = N_particles_FBP * 4/3 * np.pi * (r_FBP_buffer**3 - r_FBP_kernel**3)
V_PyC = N_particles_FBP * 4/3 * np.pi * (r_FBP_PyC**3 - r_FBP_buffer**3)
V_compact_FBP = h_FBP_compact * np.pi * r_FBP_RPT**2 - V_particles_FBP
V_total_FBP = h_FBP_compact * np.pi * r_FBP_RPT**2
V_ratios_FBP = np.array([V_B4C, V_buffer_FBP, V_PyC, V_compact_FBP])/V_total_FBP

###########################################################################################################################################################
################################################################ Other component parameters ###############################################################
#################################################################      lengths in cm       ################################################################
###########################################################################################################################################################
r_coolant_channel = 0.794 # for large coolant channels
r_inner_reflector = 68 # 68 cm for BeO reflector, 83.5 cm for graphite reflector
r_insulation = r_inner_reflector + 5
r_barrel = r_insulation + 2
r_airgap = r_barrel + 5
r_RPV = r_airgap + 10
h_reactor = N_blocks * h_block + 2 * h_reflector
###########################################################################################################################################################
V_graphite  = ((np.sqrt(3)/2) * h_block * d_flat_to_flat**2) * N_blocks * 7 # Graphite volume of one block without any channels times number of blocks (7 * 4)
V_graphite -= N_blocks * 6 * 102 * h_block * np.pi * r_coolant_channel**2 # subtracting large coolant channel volume
V_graphite -= N_blocks * 6 * 6 * h_block * np.pi * r_fuel_channel**2 # subtracting small coolant channel volume
V_graphite -= N_blocks * 6 * 210 * 15 * h_fuel_compact * np.pi * r_fuel_channel**2 # subtracting fuel channel volume
V_graphite += N_blocks * 6 * 210 * 15 * h_fuel_compact * np.pi * (r_fuel_GR**2 - r_fuel_RPT**2) # adding fuel channel graphite ring
V_graphite -= N_blocks * 6 * 6 * h_FBP_channel * np.pi * r_FBP_channel**2 # subtracting FBP channel volume
V_graphite += N_blocks * 6 * 6 * h_FBP_channel * np.pi * (r_FBP_GR**2 - r_FBP_RPT**2) # adding FBP channel graphite ring

###########################################################################################################################################################
######################################################################## Materials ########################################################################
###########################################################################################################################################################
materials_list = []
'''
The .add_s_alpha_beta cross-sections have to be removed from all materials that used them when creating the mixed fuel/FBP because the mix materials function
cannot mix materials that have them added. Can be added seperately to the mix material. 
Discrete representations (in other py scripts) were modelled with these cross-sections
'''
###########################################################################################################################################################
# Representative AGR-inspired kernel chemistry
C_over_U = 0.3253   # AGR-1-like
O_over_U = 1.3613   # AGR-1-like
rho_uco  = 10.923   # g/cm3, AGR-1-like kernel density
UCO = openmc.Material(name="AGR-like UCO kernel")
UCO.set_density("g/cm3", rho_uco)
# Uranium vector
UCO.add_nuclide("U235", enrichment, percent_type="ao")
UCO.add_nuclide("U238", 1.0 - enrichment, percent_type="ao")
# Add oxygen and carbon relative to 1 uranium atom
UCO.add_element("O", O_over_U, percent_type="ao")
UCO.add_element("C", C_over_U, percent_type="ao")
#materials_list.append(UCO)
UCO.depletable = False 
###########################################################################################################################################################
coolant = openmc.Material(name = "Coolant")
coolant.add_element("He", 1)
coolant.set_density("kg/m3", 2.86) # Used the ideal gas law, P/(R_S*T), with specific gas constant R_S = 2077 J/kg-K
coolant.temperature = T_coolant_average
materials_list.append(coolant)
###########################################################################################################################################################
'''helium material is used for gaps and other volumes that don't correspond to coolant'''
helium = openmc.Material(name = "Helium")
helium.add_element("He", 1)
helium.set_density("g/cm3", 0.00200) # from graphite vs BeO paper
helium.temperature = T_helium
materials_list.append(helium)
###########################################################################################################################################################
graphite = openmc.Material(name = "Graphite")
graphite.add_element("C", 1)
graphite.set_density("g/cm3", 1.74) # from general atomics report
graphite.add_s_alpha_beta("c_Graphite")
graphite.temperature = T_graphite
graphite.volume = V_graphite
graphite.depletable = False
materials_list.append(graphite)
###########################################################################################################################################################
matrix = openmc.Material(name = "Fuel Compact Matrix")
matrix.add_element("C", 1)
matrix.set_density("g/cm3", 1.2)
# matrix.add_s_alpha_beta('c_Graphite')
materials_list.append(matrix)
###########################################################################################################################################################
buffer = openmc.Material(name = 'Buffer')
buffer.set_density('g/cm3', 1.0)
buffer.add_element('C', 1.0)
# buffer.add_s_alpha_beta('c_Graphite')
materials_list.append(buffer)
###########################################################################################################################################################
PyC1 = openmc.Material(name = 'IPyC')
PyC1.set_density('g/cm3', 1.9)
PyC1.add_element('C', 1.0)
# PyC1.add_s_alpha_beta('c_Graphite')
materials_list.append(PyC1)
###########################################################################################################################################################
PyC2 = openmc.Material(name = 'OPyC')
PyC2.set_density('g/cm3', 1.87)
PyC2.add_element('C', 1.0)
# PyC2.add_s_alpha_beta('c_Graphite')
materials_list.append(PyC2)
###########################################################################################################################################################
SiC = openmc.Material(name = 'SiC')
SiC.set_density('g/cm3', 3.2)
SiC.add_element('C', 0.5)
SiC.add_element('Si', 0.5)
materials_list.append(SiC)
###########################################################################################################################################################
fuel = openmc.Material.mix_materials([UCO, buffer, PyC1, SiC, PyC2, matrix], V_ratios, 'vo', name = "Homogeneous Fuel") # homogeneous fuel
fuel.id = 10
#fuel.add_s_alpha_beta('c_Graphite')
#fuel.add_s_alpha_beta("c_U_in_UO2")
fuel.temperature = T_fuel
fuel.depletable = True
fuel.volume = 6 * N_blocks * 210 * 15 * h_fuel_compact * np.pi * r_fuel_RPT**2 # 6 fuel block columns, 4 blocks per column, 210 fuel channels per block, 15 fuel compacts per channel
materials_list.append(fuel)
###########################################################################################################################################################
B4C = openmc.Material(name = "Boron Carbide")
B4C.add_element("B", 4)
B4C.add_element("C", 1)
B4C.set_density('g/cm3', 2.47)
materials_list.append(B4C)
###########################################################################################################################################################
"""
burnable_poison = openmc.Material.mix_materials([B4C, buffer, PyC2, matrix], V_ratios_FBP, 'vo', name = "Homogeneous FBP") # homogeneous burnable poison
burnable_poison.add_s_alpha_beta('c_Graphite')
burnable_poison.temperature = T_FBP
burnable_poison.volume = 6 * N_blocks * 6 * h_FBP_channel * np.pi * r_FBP_RPT**2 # 6 fuel block columns, 4 blocks per column, 6 FBP channels per block with total length h_FBP_channel
burnable_poison.depletable = False
materials_list.append(burnable_poison)
"""
burnable_poison = openmc.Material(name="FBP replaced by graphite")
burnable_poison.add_element("C", 1)
burnable_poison.set_density("g/cm3", 1.74)
burnable_poison.add_s_alpha_beta("c_Graphite")
burnable_poison.temperature = T_FBP
burnable_poison.volume = (6 * N_blocks * 6 *h_FBP_channel *np.pi * r_FBP_RPT**2)
burnable_poison.depletable = False
materials_list.append(burnable_poison)
###########################################################################################################################################################
BeO = openmc.Material(name = "Beryllium Oxide")
BeO.add_element("Be", 1)
BeO.add_element("O", 1)
BeO.set_density('g/cm3', 3.02)
BeO.temperature = T_reflector
materials_list.append(BeO)
###########################################################################################################################################################
insulation = openmc.Material(name = 'Insulation SiC')
insulation.set_density('g/cm3', 3.2)
insulation.add_element('C', 0.5)
insulation.add_element('Si', 0.5)
insulation.temperature = T_insulation
materials_list.append(insulation)
###########################################################################################################################################################
steel = openmc.Material(name = "Steel") # composition from SA508 Gr.3 from Characterization of high strength and high toughness Ni–Mo–Cr low alloy steels for nuclear application
steel.add_element('C', 0.25, 'wo')
steel.add_element('Mn', 1.35, 'wo')
steel.add_element('P', 0.025, 'wo')
steel.add_element('S', 0.025, 'wo')
steel.add_element('Si', 0.30, 'wo')
steel.add_element('Ni', 0.80,'wo')
steel.add_element('Cr', 0.25, 'wo')
steel.add_element('Mo', 0.525, 'wo')
steel.add_element('Cu', 0.20, 'wo')
steel.add_element('V', 0.05, 'wo')
steel.add_element('Fe', 96.225, 'wo')
steel.set_density('g/cm3', 8)
steel.depletable = False
steel.volume = np.pi * h_reactor * (r_barrel**2 - r_insulation**2 + r_RPV**2 - r_airgap**2)
materials_list.append(steel)
###########################################################################################################################################################
water = openmc.Material(name = "Water")
water.add_element("H", 2.0)
water.add_element("O", 1.0)
water.set_density("g/cm3", 1.0)
# water.add_s_alpha_beta('c_H_in_H2O')
water.temperature = T_water
materials_list.append(water)
###########################################################################################################################################################
materials = openmc.Materials(materials_list)

###########################################################################################################################################################
#################################################################### Defining Geometry ####################################################################
###########################################################################################################################################################
def FBP_compact():
    '''
    Output:
        OpenMC universe of a cylindrical homogeneous FBP compact surrounded by a graphite ring and helium gap in an infinite graphite volume (infinite in x and y)
    ''' 
    # Defining surfaces
    compact_cylinder = openmc.ZCylinder(r = r_FBP_RPT, boundary_type = "transmission")
    gap_cylinder = openmc.ZCylinder(r = r_FBP_channel-t_FBP_gap, boundary_type = "transmission")
    channel_cylinder = openmc.ZCylinder(r = r_FBP_channel, boundary_type = "transmission")
    z_bottom = openmc.ZPlane(-h_FBP_compact/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_FBP_compact/2, boundary_type = "transmission")
    # Defining regions
    fbp_region = -compact_cylinder &+ z_bottom &- z_top
    graphite_region = +compact_cylinder &- gap_cylinder &+ z_bottom &- z_top
    gap_region = +gap_cylinder &- channel_cylinder &+ z_bottom &- z_top
    outer_region = +channel_cylinder &+ z_bottom &- z_top
    # creating cells
    compact_cell = openmc.Cell(fill = burnable_poison, region = fbp_region)
    graphite_cell = openmc.Cell(fill = graphite, region = graphite_region)
    gap_cell = openmc.Cell(fill = helium, region = gap_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    # returning universe
    return openmc.Universe(cells = [compact_cell, graphite_cell, gap_cell, outer_cell])
###########################################################################################################################################################
def fuel_compact():
    '''
    Output:
        OpenMC universe of a cylindrical homogeneous fuel compact surrounded by a graphite ring and helium gap in an infinite graphite volume (infinite in x and y)
    ''' 
    # Defining surfaces
    compact_cylinder = openmc.ZCylinder(r = r_fuel_RPT, boundary_type = "transmission")
    gap_cylinder = openmc.ZCylinder(r = r_fuel_channel-t_fuel_gap, boundary_type = "transmission")
    channel_cylinder = openmc.ZCylinder(r = r_fuel_channel, boundary_type = "transmission")
    z_bottom = openmc.ZPlane(-h_fuel_compact/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_fuel_compact/2, boundary_type = "transmission")
    # Defining regions
    fuel_region = -compact_cylinder &+ z_bottom &- z_top
    graphite_region = +compact_cylinder &- gap_cylinder &+ z_bottom &- z_top
    gap_region = +gap_cylinder &- channel_cylinder &+ z_bottom &- z_top
    outer_region = +channel_cylinder &+ z_bottom &- z_top
    # creating cells
    compact_cell = openmc.Cell(fill = fuel, region = fuel_region)
    graphite_cell = openmc.Cell(fill = graphite, region = graphite_region)
    gap_cell = openmc.Cell(fill = helium, region = gap_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    # returning universe
    return openmc.Universe(cells = [compact_cell, graphite_cell, gap_cell, outer_cell])
###########################################################################################################################################################
def coolant_channel(size, flooded):
    '''
    Input:
        size: if "large" the larger radius of the coolant channel is used, while for "small" the fuel channel radius is used
        flooded: if True the coolant channel is filled with water instead of coolant
    Output:
        OpenMC universe corresponding to a coolant channel running through the fuel block
    '''
    # Defining surfaces
    if size == "large":
        channel_cylinder = openmc.ZCylinder(r = r_coolant_channel, boundary_type = "transmission")
    else:
        channel_cylinder = openmc.ZCylinder(r = r_fuel_channel, boundary_type = "transmission")
    z_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    # Defining regions
    coolant_region = - channel_cylinder &+ z_bottom &- z_top
    outer_region = + channel_cylinder  &+ z_bottom &- z_top
    # Creating cells
    if flooded == False:
        coolant_cell = openmc.Cell(fill = coolant, region = coolant_region)
    else:
        coolant_cell = openmc.Cell(fill = water, region = coolant_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    # Returning universe
    return openmc.Universe(cells = [coolant_cell, outer_cell])
###########################################################################################################################################################
def graphite_pin():
    '''
    Output:
        OpenMC universe of a graphite pin for the center of the fuel blocks
    '''
    # Defining surfaces
    channel_cylinder = openmc.ZCylinder(r = r_fuel_channel, boundary_type = "transmission")
    z_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    # Defining regions
    pin_region = - channel_cylinder &+ z_bottom &- z_top
    outer_region = + channel_cylinder  &+ z_bottom &- z_top
    # Creating cells
    pin_cell = openmc.Cell(fill = graphite, region = pin_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    # Returning universe
    return openmc.Universe(cells = [pin_cell, outer_cell])
###########################################################################################################################################################
def FBP_channel():
    '''
    Output:
        OpenMC universe containing several FBP compacts stacked ontop of oneanother and capped off with graphite plugs. Overall length corresponds to fuel
        block height
    '''
    # loading in the fuel compact
    compact = FBP_compact()
    # Cell for filling the assembly/outer universe
    outer_cell = openmc.Cell(fill = graphite)
    outer = openmc.Universe(cells = [outer_cell])
    # Defining rectangular lattice to stack fuel compacts
    lat = openmc.HexLattice()
    lat.center = (0, 0, 0)
    lat.orientation = 'y'
    lat.pitch = (pitch_holes, h_FBP_compact)
    lat.outer = outer
    lat.universes = [[[compact]]]*N_FBP
    # Defining surfaces
    z_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_bottom_plug = openmc.ZPlane(-h_FBP_compact*N_FBP/2, boundary_type = "transmission")
    z_top_plug = openmc.ZPlane(h_FBP_compact*N_FBP/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    channel_cylinder = openmc.ZCylinder(r = r_FBP_channel, boundary_type = "transmission")
    # Defining regions
    bottom_plug_region = - channel_cylinder &+ z_bottom &- z_bottom_plug
    FBP_channel_region = - channel_cylinder &+ z_bottom_plug &- z_top_plug
    top_plug_region = - channel_cylinder &+ z_top_plug &- z_top
    outer_region = + channel_cylinder &+ z_bottom &- z_top
    # Defining cells and returning universe
    bottom_plug = openmc.Cell(fill = graphite, region = bottom_plug_region)
    FBP_channel = openmc.Cell(fill = lat, region = FBP_channel_region)
    top_plug = openmc.Cell(fill = graphite, region = top_plug_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    return openmc.Universe(cells = [bottom_plug, FBP_channel, top_plug, outer_cell])
###########################################################################################################################################################
def fuel_channel():
    '''
    Output:
        OpenMC universe containing 15 fuel compacts stacked on top of another, capped off with two graphite plugs. Overall length corresponds to fuel
        block height
    '''
    # loading in the fuel compact
    compact = fuel_compact()
    # Cell for filling the assembly/outer universe
    outer_cell = openmc.Cell(fill = graphite)
    outer = openmc.Universe(cells = [outer_cell])
    # Defining rectangular lattice to stack fuel compacts
    lat = openmc.HexLattice()
    lat.center = (0, 0, 0)
    lat.orientation = 'y'
    lat.pitch = (pitch_holes, h_fuel_compact)
    lat.outer = outer
    lat.universes = [[[compact]]]*15
    # Defining surfaces
    z_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_bottom_plug = openmc.ZPlane(-h_fuel_compact*15/2, boundary_type = "transmission")
    z_top_plug = openmc.ZPlane(h_fuel_compact*15/2, boundary_type = "transmission")
    z_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    channel_cylinder = openmc.ZCylinder(r = r_fuel_channel, boundary_type = "transmission")
    # Defining regions
    bottom_plug_region = - channel_cylinder &+ z_bottom &- z_bottom_plug
    fuel_channel_region = - channel_cylinder &+ z_bottom_plug &- z_top_plug
    top_plug_region = - channel_cylinder &+ z_top_plug &- z_top
    outer_region = + channel_cylinder &+ z_bottom &- z_top
    # Defining cells and returning universe
    bottom_plug = openmc.Cell(fill = graphite, region = bottom_plug_region)
    fuel_channel = openmc.Cell(fill = lat, region = fuel_channel_region)
    top_plug = openmc.Cell(fill = graphite, region = top_plug_region)
    outer_cell = openmc.Cell(fill = graphite, region = outer_region)
    return openmc.Universe(cells = [bottom_plug, fuel_channel, top_plug, outer_cell])
###########################################################################################################################################################
def fuel_block(flooded, FBP):
    '''
    Function for creating a single graphite 'assembly' block containing 210 fuel channels with fuel compacts, 108 coolant channels (6 of which are smaller), 
    and (potentially) 6 FBP channels
    Input:
        flooded: if True, the coolant channels are filled with water rather than coolant
        FBP: if True, adds 6 FBP channels in the corners of the hexagonal block. If false, adds only 0
    Output:
        OpenMC universe of a hexagonal graphite 'assembly' block
    '''
    # Making the pins for the block
    empty_pin = graphite_pin()
    coolant_pin_L = coolant_channel('large', flooded)
    coolant_pin_S = coolant_channel('small', flooded)
    fuel_pin = fuel_channel()
    FBP_pin = FBP_channel()
    if flooded == True:
        gap_material = water
    else:
        gap_material = coolant
    # Cell for filling the aasembly/outer universe
    outer_cell = openmc.Cell(fill = graphite)
    outer = openmc.Universe(cells = [outer_cell])
    # Making the hex lattice
    lat = openmc.HexLattice()
    lat.center = (0, 0)
    lat.pitch = (pitch_holes, )
    lat.orientation = 'x'
    lat.outer = outer
    # filling the hex lattice
    ring_list = [[empty_pin], [empty_pin]*6]
    ring3 = []
    for n in range(1,13):
        if n % 2 == 0:
            ring3.append(coolant_pin_S)
        else:
            ring3.append(fuel_pin)
    ring_list.append(ring3)
    ring4 = []
    for n in range(1,19):
        if n in [x for x in range(1, 19, 3)]:
            ring4.append(coolant_pin_L)
        else:
            ring4.append(fuel_pin)
    ring_list.append(ring4)
    ring5 = []
    for n in range(1,25):
        if n in [x for x in range(3, 25, 4)]:
            ring5.append(coolant_pin_L)
        else:
            ring5.append(fuel_pin)
    ring_list.append(ring5)
    ring6 = []
    for n in range(1,31):
        if n in [2,5,7,10,12,15,17,20,22,25,27,30]:
            ring6.append(coolant_pin_L)
        else:
            ring6.append(fuel_pin)
    ring_list.append(ring6)
    ring7 = []
    for n in range(1,37):
        if n in [x for x in range(1, 37, 3)]:
            ring7.append(coolant_pin_L)
        else:
            ring7.append(fuel_pin)
    ring_list.append(ring7)
    ring8 = []
    for n in range(1,43):
        if n in [3,6,10,13,17,20,24,27,31,34,38,41]:
            ring8.append(coolant_pin_L)
        else:
            ring8.append(fuel_pin)
    ring_list.append(ring8)
    ring9 = []
    for n in range(1,49):
        if n in [2,5,8,10,13,16,18,21,24,26,29,32,34,37,40,42,45,48]:
            ring9.append(coolant_pin_L)
        else:
            ring9.append(fuel_pin)
    ring_list.append(ring9)
    ring10 = []
    for n in range(1,55):
        if n in [x for x in range(1, 55, 3)]:
            ring10.append(coolant_pin_L)
        else:
            ring10.append(fuel_pin)
    ring_list.append(ring10)
    ring11 = []
    for n in range(1,61):
        if n in [3,6,9,13,16,19,23,26,29,33,36,39,43,46,49,53,56,59]:
            ring11.append(coolant_pin_L)
        elif n in [1,11,21,31,41,51] and FBP == True:
            ring11.append(FBP_pin)
        else:
            ring11.append(fuel_pin)
    ring_list.append(ring11)
    # Creating Lattice universe
    ring_list.reverse()
    lat.universes = ring_list
    # Boundary surfaces
    z_block_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_block_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    hex_inner = openmc.model.HexagonalPrism(d_flat_to_flat/np.sqrt(3), orientation = 'x', boundary_type = "transmission")
    hex_outer = openmc.model.HexagonalPrism(pitch_assembly/np.sqrt(3), orientation = 'x', boundary_type = "transmission")
    # Regions
    inner_region = - hex_inner &+ z_block_bottom &- z_block_top
    outer_region = - hex_outer &+ hex_inner &+ z_block_bottom &- z_block_top
    # Making cells and reutning the block universe
    fuel_block = openmc.Cell(fill = lat, region = inner_region)
    block_gap = openmc.Cell(fill = gap_material, region = outer_region)
    return openmc.Universe(cells = [fuel_block, block_gap])
###########################################################################################################################################################
def reflector_block(flooded):
    '''
    Input:
        flooded: if True, the coolant regions (gap between blocks) are filled with water rather than coolant
    Output:
        OpenMC universe of a hexagonal graphite block acting as a central reflector
    '''
    if flooded == True:
        gap_material = water
    else:
        gap_material = coolant
    # Defining surfaces:
    z_block_bottom = openmc.ZPlane(-h_block/2, boundary_type = "transmission")
    z_block_top = openmc.ZPlane(h_block/2, boundary_type = "transmission")
    hex_inner = openmc.model.HexagonalPrism(d_flat_to_flat/np.sqrt(3), orientation = 'x', boundary_type = "transmission")
    hex_outer = openmc.model.HexagonalPrism(pitch_assembly/np.sqrt(3), orientation = 'x', boundary_type = "transmission")
    # Defining regions
    reflector_region = -hex_inner &+ z_block_bottom &- z_block_top
    gap_region = -hex_outer &+ hex_inner &+ z_block_bottom &- z_block_top
    # Definign cells and returning universe
    reflector_cell = openmc.Cell(fill = graphite, region = reflector_region)
    gap_cell = openmc.Cell(fill = gap_material, region = gap_region)
    return openmc.Universe(cells = [reflector_cell, gap_cell])
###########################################################################################################################################################
def make_reactor(N_blocks, N_flooded, reflector_material = BeO):
    '''
    Function which creates the OpenMC universe describing the entire reactor geometry
    Input:
        N_blocks: number of blocks per fuel column (z-direction). Integer 
        N_flooded: number of blocks (vertically) that are to be flooded. Should be an integer beween 0 and N_blocks
        h_reflector: height of the top/bottom reflectors. In cm
        reflector_material: material used for the radial and axial reflectors. Default is BeO. openmc.material
    Output:
        reactor_universe: OpenMC universe corresponding to the reactor geometry
        reactor_cells: list of cells that make up the reactor universe. Useful as a domain input for certain tallys/filters
    '''
    # Relevant assemblies
    reflector_flooded = reflector_block(True)
    reflector = reflector_block(False)
    fuel_6_flooded = fuel_block(True,False) # changed second argument to false to remove FBP
    fuel_6 = fuel_block(False, False) # changed second argument to false to remove FBP
    fuel_0_flooded = fuel_block(True, False)
    fuel_0 = fuel_block(False, False)
    # Empty cell/universe for filling the lattice and final reactor universe
    outer_cell = openmc.Cell(fill = BeO) # graphite or BeO
    outer = openmc.Universe(cells = [outer_cell])
    outer_cell2 = openmc.Cell(fill = None) # graphite or BeO
    outer2 = openmc.Universe(cells = [outer_cell2])
    # making the reactor lattice
    lat = openmc.HexLattice()
    lat.center = (0, 0, 0)
    lat.orientation = 'y'
    lat.pitch = (pitch_assembly, h_block)
    lat.outer = outer
    # Filling the reactor lattice
    ring_list = []
    '''
    These lists were the ones that were changed for the FBP k_eff tests
    '''
    if N_flooded == 0:
        ring_list.append([6*[fuel_0], [reflector]])
        for n in range(N_blocks-2):
            ring_list.append([6*[fuel_6], [reflector]])
        ring_list.append([6*[fuel_0], [reflector]])
    else:
        for n in range(N_flooded):
            if n == 0 or n == 3:
                ring_list.append([6*[fuel_0_flooded], [reflector_flooded]])
            else:
                ring_list.append([6*[fuel_6_flooded], [reflector_flooded]])
        for n in range(N_blocks - N_flooded):
            if n == N_blocks - N_flooded - 1:
                ring_list.append([6*[fuel_0], [reflector]])
            else:
                ring_list.append([6*[fuel_6], [reflector]])
    # Setting lattice universe
    lat.universes = ring_list
    # Defining surfaces
    z_reactor_bottom = openmc.ZPlane(-h_reactor/2, boundary_type = "vacuum")
    z_fuel_bottom = openmc.ZPlane(-N_blocks*h_block/2, boundary_type = "transmission")
    z_fuel_top = openmc.ZPlane(N_blocks*h_block/2, boundary_type = "transmission")
    z_reactor_top = openmc.ZPlane(h_reactor/2, boundary_type = "vacuum")
    reflector_cylinder = openmc.ZCylinder(r = r_inner_reflector, boundary_type = "transmission")
    insulation_cylinder = openmc.ZCylinder(r = r_insulation, boundary_type = "transmission")
    barrel_cylinder = openmc.ZCylinder(r = r_barrel, boundary_type = "transmission")
    gap_cylinder = openmc.ZCylinder(r = r_airgap, boundary_type = "transmission")
    RPV_cylinder = openmc.ZCylinder(r = r_RPV, boundary_type = "vacuum")
    # Reator region
    reactor_region = - reflector_cylinder &+ z_fuel_bottom &- z_fuel_top
    UR_region = - reflector_cylinder &+ z_fuel_top &- z_reactor_top
    LR_region = - reflector_cylinder &+ z_reactor_bottom &- z_fuel_bottom
    SR_region = - reflector_cylinder &+ z_reactor_bottom &- z_reactor_top
    insulation_region = +reflector_cylinder &- insulation_cylinder &+ z_reactor_bottom &- z_reactor_top
    barrel_region = + insulation_cylinder &- barrel_cylinder &+ z_reactor_bottom &- z_reactor_top
    gap_region = + barrel_cylinder &- gap_cylinder &+ z_reactor_bottom &- z_reactor_top
    RPV_region = + gap_cylinder &- RPV_cylinder &+ z_reactor_bottom &- z_reactor_top
    # Reactor cells
    reactor_cell = openmc.Cell(fill = lat, region = reactor_region)
    UR_cell = openmc.Cell(fill = reflector_material, region = UR_region)
    LR_cell = openmc.Cell(fill = reflector_material, region = LR_region)
    SR_cell = openmc.Cell(fill = reflector_material, region = SR_region)
    insulation_cell = openmc.Cell(fill = insulation, region = insulation_region)
    barrel_cell = openmc.Cell(fill = steel, region = barrel_region)
    barrel_cell.temperature = T_barrel
    gap_cell = openmc.Cell(fill = helium, region = gap_region)
    RPV_cell = openmc.Cell(fill = steel, region = RPV_region)
    RPV_cell.temperature = T_RPV
    # Reactor universe
    reactor_cells = [reactor_cell, UR_cell, LR_cell, SR_cell, insulation_cell, barrel_cell, gap_cell, RPV_cell]
    reactor_universe = openmc.Universe(cells = reactor_cells)
    reactor_universe.outer = outer2
    return reactor_universe, reactor_cells

###########################################################################################################################################################
############################################################## Creating the Reactor Geomgetry #############################################################
###########################################################################################################################################################
reactor, reactor_domains = make_reactor(N_blocks, N_flooded)
block = fuel_block(False, True)
geometry = openmc.Geometry()
geometry.root_universe = reactor
# geometry.root_universe = block

###########################################################################################################################################################
################################################################# Setting up Tallies File #################################################################
###########################################################################################################################################################
"""tallies_file = openmc.Tallies() 
###########################################################################################################################################################
'''
Since our simulation is being run in neutron-transport only mode, used the heating-local score
'''
tally_heating = openmc.Tally(name = "Heating")
tally_heating.scores = ["heating-local"]
tallies_file.append(tally_heating)
###########################################################################################################################################################
axial_mesh = openmc.RegularMesh()
axial_mesh.dimension = [1, 1, 96]
'''
32 because I'm ignoring the axial refletors (no power production there) which gives 320/20 = 16, and because I wanted a slightly finer mesh I simply multiplied 16 by a whole number
'''
axial_mesh.lower_left = [-r_RPV, -r_RPV, -h_reactor/2 + h_reflector]
axial_mesh.upper_right = [r_RPV, r_RPV, h_reactor/2 - h_reflector]
axial_mesh_filter = openmc.MeshFilter(axial_mesh)
tally_AP = openmc.Tally(name = "Axial Power")
tally_AP.filters = [axial_mesh_filter]
tally_AP.scores = ["fission-q-recoverable"]
tallies_file.append(tally_AP)
###########################################################################################################################################################
r_grid = np.linspace(0, r_inner_reflector, 30)
z_grid = np.array([-N_blocks*h_block/2, N_blocks*h_block/2])
phi_grid = np.array([0, 2*np.pi])
radial_mesh = openmc.CylindricalMesh(r_grid, z_grid, phi_grid, origin = (0,0,0))
radial_mesh_filter = openmc.MeshFilter(radial_mesh)
tally_RP = openmc.Tally(name = "Radial Power")
tally_RP.filters = [radial_mesh_filter]
tally_RP.scores = ["fission-q-recoverable"]
tallies_file.append(tally_RP)
"""
###########################################################################################################################################################
############################################################## Creating the Settings xml File #############################################################
###########################################################################################################################################################
'''
General particle/batch settings unless otherwise specified in report:
300 batches, 50 inactive, 10,000 particles (could use more for single k-eigenvalue calculations)
Turn off photon transport
'''
settings = openmc.Settings()
settings.batches = 40
settings.inactive = 10
settings.particles = 5000
settings.photon_transport = False
settings.temperature = {'method': 'interpolation', 'range': (250, 2500)}
bounds = [-r_RPV,-r_RPV,-h_reactor/2,r_RPV,r_RPV, h_reactor/2] # box which encapsules the reactor pressure vessel and all internal components
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(space = uniform_dist, constraints = {'domains': reactor_domains}) # uniform distribution
###########################################################################################################################################################

def make_model(T_fuel_in: float = 1200.0) -> openmc.Model:
    # Oppdater fuel-temp
    fuel.temperature = T_fuel_in

    # Bygg og returner model (uten tallies)
    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
    model.export_to_xml()
    return model

if __name__ == "__main__":
    import os
    import argparse
    import numpy as np
    import openmc
    import openmc.deplete

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["keff", "deplete"],
        default="keff",
        help="What to run: keff (eigenvalue) or deplete (depletion).",
    )

    parser.add_argument("--T_fuel", type=float, default=T_fuel, help="Fuel temperature in K")
    parser.add_argument("--power", type=float, default=1.0e7, help="Reactor power in W (thermal)")

    # Always interpret these as DAYS (we will convert to seconds internally for robustness)
    parser.add_argument(
        "--timesteps",
        type=float,
        nargs="+",
        default=[0.05],
        help="Depletion timesteps in DAYS, e.g. --timesteps 30 30 30",
    )

    parser.add_argument(
        "--chain",
        type=str,
        default="chain_endfb71_pwr.xml",
        help="Path to depletion chain XML (relative to run dir or absolute).",
    )

    args = parser.parse_args()

    # Build model with chosen fuel temperature
    model = make_model(T_fuel_in=args.T_fuel)

    if args.mode == "keff":
        model.run()

    elif args.mode == "deplete":
        # ----------------------------
        # Idiot-safe timestep handling
        # ----------------------------
        # Interpret CLI timesteps as DAYS
        timesteps_days = np.array(args.timesteps, dtype=float)

        if timesteps_days.ndim != 1 or len(timesteps_days) == 0:
            raise ValueError("No timesteps provided. Example: --timesteps 30 30 30")

        if np.any(timesteps_days <= 0.0):
            raise ValueError(f"All timesteps must be > 0 days. Got: {timesteps_days.tolist()}")

        # Convert DAYS -> SECONDS (robust across OpenMC builds; Results store seconds)
        timesteps_s = (timesteps_days * 86400.0).astype(float).tolist()

        # One power value per timestep
        power_level = [float(args.power)] * len(timesteps_s)

        # Sanity warning: if steps look tiny in seconds, user may have passed seconds by mistake
        if max(timesteps_s) < 3600.0:
            print(
                "WARNING: Largest timestep is < 1 hour (in seconds). "
                "This script interprets --timesteps as DAYS."
            )

        total_days = float(np.sum(timesteps_days))
        total_years = total_days / 365.25
        total_seconds = float(np.sum(timesteps_s))

        # Record what we ran (clear and unambiguous)
        with open("case_info.txt", "w") as f:
            f.write(f"T_fuel={args.T_fuel}\n")
            f.write(f"power_W={float(args.power)}\n")
            f.write("timesteps_input_unit=days\n")
            f.write(f"timesteps_days={timesteps_days.tolist()}\n")
            f.write(f"timesteps_seconds={timesteps_s}\n")
            f.write(f"total_days={total_days}\n")
            f.write(f"total_years={total_years}\n")
            f.write(f"total_seconds={total_seconds}\n")
            f.write(f"chain={args.chain}\n")
            f.write(f"OPENMC_CROSS_SECTIONS={os.environ.get('OPENMC_CROSS_SECTIONS')}\n")

        print("DEBUG: Depletion will run with timestep_units='s' (forced for robustness).")
        print("DEBUG: timesteps_days    =", timesteps_days.tolist())
        print("DEBUG: timesteps_seconds =", timesteps_s)
        print(f"DEBUG: total = {total_days:.6f} days = {total_years:.6f} years")

        operator = openmc.deplete.CoupledOperator(model, chain_file=args.chain)
        integrator = openmc.deplete.PredictorIntegrator(
            operator,
            timesteps_s,
            power_level,
            timestep_units="s",
        )
        integrator.integrate()
