import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import ert
from functions import *

tree = 'Tree835_1_220511_morning'
hardwood_radius = 0.0891

GEOM = np.loadtxt('data/geometry_835.txt')
nel = GEOM.shape[0]  # n° of electrodes
height = 0.5  # height of the tree
zel = 0.  # electrodes z coordinates
DATA = np.loadtxt(strip_first_col('data/220511_Tree_nr_835_1.txt'))
N = DATA.shape[0]  # n° of data points
med_u, med_i = med_u_i('data/220511_Tree_nr_835_1.txt',
                       'data/220511_Tree_nr_835_2.txt',
                       'data/220511_Tree_nr_835_3.txt')
err, med_rhoa = error_calc('data/220511_Tree_nr_835_1.txt',
                           'data/220511_Tree_nr_835_2.txt',
                           'data/220511_Tree_nr_835_3.txt')

# create Dataset from geometry and measured data:
with open('data/Tree_Dataset.dat', 'w') as f:  #
    f.write('%d' % nel + '\n')  #
    f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
    for i in range(nel):  #
        f.write('%.5f ' % GEOM[i, 1] + '%.5f ' % GEOM[i, 2] + '%.5f ' % zel + '\n')  #
    # write data                                                                                                                                                                            #
    f.write('%d' % N + '\n')  #
    f.write('# ' + 'a ' + 'b ' + 'm ' + 'n ' + 'rhoa ' + 'u ' + 'i ' + 'k ' + '\n')  #
    for i in range(N):  #
        f.write(
            '%d ' % DATA[i, 0] + '%d ' % DATA[i, 1] + '%d ' % DATA[i, 2] + '%d ' % DATA[i, 3] + '%.6f ' % np.abs(DATA[
                                                                                                                     i, 6]) + '%.6f ' % np.abs(
                DATA[i, 5]) + '%.6f ' % np.abs(DATA[i, 4]) + '%.6f ' % 1 + '\n')  #
    f.write('%d' % 0)  #
f.close()  #
#########################################################################################################################################################################################

# 1) Crete the tree geometry (cylinder) from electrodes locations
area = 0.1  # maximum cell size for resulting triangles after mesh generation
quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
# refi_node = 5  # n° of additional nodes for refinement between electrodes
center_pos = [GEOM[:, 1].min() + (GEOM[:, 1].max() - GEOM[:, 1].min()) / 2,
              GEOM[:, 2].min() + (GEOM[:, 2].max() - GEOM[:, 2].min()) / 2, 0.0]
El_geom = mt.createPolygon(GEOM[:, 1:3],
                           isClosed=True,
                           area=area,
                           quality=quality,
                           boundaryMarker=1)  # , addNodes=refi_node)
Tree_Geom = mt.extrude(El_geom, z=height)  # marker=2)
Tree_Geom.translate([0, 0, -height / 2])

# 2) Set the ERT position (along the circumference and equally spaced) and create mesh
for i in range(nel):
    # set the electrodes as nodes with marker -99 to the geometry
    # Index = Tree_Geom.findNearestNode((GEOM[i,1], GEOM[i,2], zel))
    # Tree_Geom.node(Index).setMarker(-99)
    Tree_Geom.createNode(pg.Pos(GEOM[i, 1], GEOM[i, 2], zel), marker=-99)
    # For sufficient numerical accuracy it is generally a good idea to refine the mesh in the vicinity of the electrodes positions.
    Tree_Geom.createNode([GEOM[i, 1], GEOM[i, 2], zel - 1e-3 / 2])
# Always need dipole current injection since there can be no current flow out of the closed boundaries.
# Define a reference electrode position inside the PLC, with a marker -999, somewhere away from the electrodes (and refine it).
Tree_Geom.createNode(center_pos, marker=-999)  # center_pos
Tree_Geom.createNode([center_pos[0], center_pos[1], center_pos[2] - 1e-3 / 2])
# The second problem for pure Neumann domains is the non-uniqueness of the partial differential equation (there are only partial derivatives of the electric potential so an arbitrary value might be added, i.e. calibrated).
# Add calibration node with marker -1000 where the potential is fixed , somewhere on the boundary and far from the electrodes.
Tree_Geom.createNode([center_pos[0], center_pos[1], height/2], marker=-1000)
Tree = mt.createMesh(Tree_Geom)

# 3) Simulate response for homogeneous resistivities distributions
hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)
plcMod = Tree_Geom + hardwood_cylinder
HetTree = mt.createMesh(plcMod)

# 4) Set the ERT scheme (from field geometry) and simulate response from a uniform resistivity (to compute geometric
# factor) OBS: The apparent resistivity for a homogeneous model of 1 Ohmm should be 1 Ohmm. Therefore we can take the
# inverse of the modeled resistances for the homogeneous model and use it as geometric factors to find the apparent
# resistivities for the inhomogeneous model.
DataSet = ert.load('data/Tree_Dataset.dat')
DataSet.set('u', med_u)
DataSet.set('i', med_i)
DataSet.set('rhoa', med_rhoa)

# When setting the calcOnly flag on, the potential (Homogeneous['u']) turns out to be negative
# and so does the calculation of k. The apparent resistivity (Homogeneous['rhoa']) is
# identical to the data (so no simulation saved)
Homogeneous = ert.simulate(Tree, res=1, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
k = 1.0 * Homogeneous('i') / Homogeneous('u')

# When the calcOnly flag is false (response saved in Homogeneous['rhoa']) then
# the output seems ok, and we can use the inverse of Rho App to compute k.
# There is no minus sign resulting here
Homogeneous = ert.simulate(Tree, res=1, scheme=DataSet, sr=False, calcOnly=False, verbose=True)
k = 1.0 / Homogeneous('rhoa')
DataSet.set('k', k)

# Simulation output seems to suggest
# Response: min = -239.715 max = -3.35043 mean = -36.3533
# but actual values do not agree
Homogeneous['rhoa'].array().min()
Homogeneous['rhoa'].array().max()

plt.figure()
plt.plot(DataSet['rhoa'].array(), 'x')
plt.plot(Homogeneous['rhoa'].array(),'o')
plt.xlabel('Datapoint (-)')
plt.ylabel('Apparent resistivity (Ohm meter) (-)')
plt.legend(['Measured data', 'Simulated data'])
plt.show()


DataSet.set('err', err)
DataSet.save('data/Tree_Dataset.dat', 'a b m n err rhoa k u i')

# Circular pseudosections:
pg.show(DataSet, circular=True)
plt.show()
