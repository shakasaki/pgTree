import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import time
import pygimli.frameworks as pgf

from functions import strip_first_col
from functions import Tree_geometry
from functions import ERT_position
from functions import error_calc
from functions import create_starting_model_3D
from functions import create_starting_model_2D
from functions import med_u_i
import pyvista as pv

start = time.time()
tree = 'Tree835_1_220511_morning'
hardwood_radius = 0.0891

################################################################################
# 0) Read geometry and data file and Write data for BERT Inversion

# TO DO BEFORE ##################################################################
# geometry.txt: change , with tab/space separation between columns              #
# measured_values.csv: save as .txt + delete first row (the one with letters)   #
################################################################################

GEOM = np.loadtxt('geometry_835.txt')
nel = GEOM.shape[0]  # n° of electrodes
height = 0.5  # height of the tree
zel = 0.  # electrodes z coordinates
DATA = np.loadtxt(strip_first_col('220511_Tree_nr_835_1.txt'))
N = DATA.shape[0]  # n° of data points
med_u, med_i = med_u_i('220511_Tree_nr_835_1.txt', '220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')
err, med_rhoa = error_calc('220511_Tree_nr_835_1.txt','220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')

#create Dataset from geometry and measured data:
with open('Tree_Dataset_835_morning_220511_1.dat', 'w') as f:  #
    f.write('%d' % nel + '\n')  #
    f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
    for i in range(nel):  #
        f.write('%.5f ' % GEOM[i, 1] + '%.5f ' % GEOM[i, 2] + '%.5f ' % zel + '\n')  #
    # write data                                                                                                                                                                            #
    f.write('%d' % N + '\n')  #
    f.write('# ' + 'a ' + 'b ' + 'm ' + 'n ' + 'rhoa ' + 'u ' + 'i ' + 'k ' + '\n')  #
    for i in range(N):  #
        f.write('%d ' % DATA[i, 0] + '%d ' % DATA[i, 1] + '%d ' % DATA[i, 2] + '%d ' % DATA[i, 3] + '%.6f ' % np.abs(DATA[
            i, 6]) + '%.6f ' % np.abs(DATA[i, 5]) + '%.6f ' % np.abs(DATA[i, 4]) + '%.6f ' % 1 + '\n')  #
    f.write('%d' % 0)  #
f.close()  #
#########################################################################################################################################################################################

# 1) Crete the tree geometry (cylinder) from electrodes locations
area = 0.1  # maximum cell size for resulting triangles after mesh generation
quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
#refi_node = 5  # n° of additional nodes for refinement between electrodes
center_pos, Tree_Geom = Tree_geometry(GEOM, area, quality, height)
Tree_Geom.translate([0, 0, -height /2])

# 2) Set the ERT position (along the circumference and equally spaced) and create mesh
Tree, Tree_Geom = ERT_position(nel, Tree_Geom, GEOM, zel, center_pos, height)

# 3) Simulate response for homogeneous resistivities distributions
hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos,  marker=2)
plcMod = Tree_Geom + hardwood_cylinder
HetTree = mt.createMesh(plcMod)

# 4) Set the ERT scheme (from field geometry) and simulate response from a uniform resistivity (to compute geometric
# factor) OBS: The apparent resistivity for a homogeneous model of 1 Ohmm should be 1 Ohmm. Therefore we can take the
# inverse of the modeled resistances for the homogeneous model and use it as geometric factors to find the apparent
# resistivities for the inhomogeneous model.
DataSet = ert.load('Tree_Dataset_835_morning_220511_1.dat')
DataSet.set('u', med_u)
DataSet.set('i', med_i)
DataSet.set('rhoa', med_rhoa)
Homogeneous = ert.simulate(Tree, res=180, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
Homogeneous2 = ert.simulate(Tree, res=180, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
k = 1.0 * Homogeneous('i') / Homogeneous('u')
DataSet.set('k', -k)

# Set different resistivities within the two cylinders (inner and outer one) and noise to the data
#res = [[1, 180.0],
      # [2, 1340.0]]  # This maps markers 1 (outer cylinder) to 10 Ohm m and markers 2 (inner cylinder) to 100 Ohm m
#resistivity values from Bieker and Rust, 2010: Non-Destructive Estimation of Sapwood and Heartwood Width in Scots Pine (Pinus sylvestris L.)
# the error (err) is imported from the error calculation file
DataSet.set('err', err)
DataSet.set('rhoa', med_rhoa)
DataSet.save('Tree_Dataset_GF_morning_220511_1.dat', 'a b m n err rhoa k u i')

# Circular pseudosections:
pg.show(DataSet, circular=True)

##################################################################################################################################################

# 2D Inversion

#load the starting model, and the homogeneous and heterogeneous mesh from the functions file
starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(DataSet, hardwood_radius)

#perform the inverion for a heterogeneous (with hardwood/sapwood boundary) and homogeneous tree
ert_2D_het = ert.ERTManager()
inv_2D_het = ert_2D_het.invert(DataSet, mesh=mesh_2D, startModel=starting_model_2D, verbose=True)
cov_2D_het = ert_2D_het.coverage()
ert_2D_het.showResult(model=inv_2D_het, coverage=cov_2D_het, cMin=60, cMax=110)#, logScale=True)
ert_2D_het.showResult(cMin=60, cMax=110, logScale=True)

#simulation1 = ert_2D_het.simulate(mesh_2D, scheme=DataSet, res=180, sr=True, calcOnly=False, verbose=True)
simulation2 = ert_2D_het.simulate(mesh_2D, scheme=DataSet, res=starting_model_2D, sr=False, calcOnly=True, verbose=True)


