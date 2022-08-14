import pandas as pd
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from core import OUTPUT_DIR, DATA_DIR
import numpy as np
from core.helpers import save_pickle_obj
import os

trees = [616, 671, 835, 836, 870, 875]
tree_dict = {}
dataset_dict = {}
plt.figure()
for tree in trees:
    electrodes = pd.read_table(DATA_DIR + 'geometries' + os.sep + 'geometry_' + str(tree) + '.txt', header=None).loc[:, 1:2].to_numpy()
    tree_dict[tree] = {}
    tree_dict[tree]['electrode positions'] = electrodes
    tree_dict[tree]['center'] = [electrodes[:, 0].min() + (electrodes[:, 0].max() - electrodes[:, 0].min()) / 2,
                                 electrodes[:, 1].min() + (electrodes[:, 1].max() - electrodes[:, 1].min()) / 2, 0.0]
    plt.scatter(tree_dict[tree]['electrode positions'][:, 0],
                tree_dict[tree]['electrode positions'][:, 1])
    plt.scatter(tree_dict[tree]['center'][0],
                tree_dict[tree]['center'][1])
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

geom_16 = pd.read_table(DATA_DIR + 'geometries' + os.sep + '16_electrode_geometry.txt', header=None).loc[:, 1:4].to_numpy()
geom_18 = pd.read_table(DATA_DIR + 'geometries' + os.sep + '18_electrode_geometry.txt', header=None).loc[:, 1:4].to_numpy()

for tree in trees:
    print('Tree ' + str(tree) + ' has ' + str(tree_dict[tree]['electrode positions'].shape[0]) + ' electrodes')
print('Trees with 16 electrodes have ' + str(geom_16.shape[0]) + ' electrode combinations')
print('Trees with 18 electrodes have ' + str(geom_18.shape[0]) + ' electrode combinations')

# Create the geometric factors

for tree in trees:
    electrode_geometry = tree_dict[tree]['electrode positions']
    nel = electrode_geometry.shape[0]  # n° of electrodes
    if nel == 16:
        electrode_arrays = geom_16
    if nel == 18:
        electrode_arrays = geom_18
    height = 1.0  # height of the tree
    zel = 0.  # electrodes z coordinates
    area = 0.1  # maximum cell size for resulting triangles after mesh generation
    quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)

    El_geom = mt.createPolygon(electrode_geometry,
                               isClosed=True,
                               area=area,
                               quality=quality,
                               boundaryMarker=1)  # , addNodes=refi_node)
    Tree_Geom = mt.extrude(El_geom, z=height)  # marker=2)
    Tree_Geom.translate([0, 0, -height / 2])

    for i in range(nel):
        # set the electrodes as nodes with marker -99 to the geometry
        # Index = Tree_Geom.findNearestNode((GEOM[i,1], GEOM[i,2], zel))
        # Tree_Geom.node(Index).setMarker(-99)
        Tree_Geom.createNode(pg.Pos(electrode_geometry[i, 0],
                                    electrode_geometry[i, 1], zel), marker=-99)
        # For sufficient numerical accuracy it is generally a good idea to refine the mesh in the vicinity of the electrodes positions.
        Tree_Geom.createNode([electrode_geometry[i, 0],
                              electrode_geometry[i, 1],
                              zel - 1e-3 / 2])

    # Always need dipole current injection since there can be no current flow out of the closed boundaries.
    # Define a reference electrode position inside the PLC, with a marker -999, somewhere away from the electrodes (and refine it).
    reference_electrode = [electrode_geometry[0, 0],
                           electrode_geometry[0, 1],
                           height / 2]
    Tree_Geom.createNode(reference_electrode, marker=-999)
    # Tree_Geom.createNode([reference_electrode[0],
    #                       reference_electrode[1],
    #                       reference_electrode[2] - 1e-3 / 2])

    # The second problem for pure Neumann domains is the non-uniqueness of the partial differential equation (there are only partial derivatives of the electric potential so an arbitrary value might be added, i.e. calibrated).
    # Add calibration node with marker -1000 where the potential is fixed , somewhere on the boundary and far from the electrodes.

    calibration_electrode = [electrode_geometry[0, 0],
                             electrode_geometry[0, 1],
                             -height / 2]
    Tree_Geom.createNode(calibration_electrode,
                         marker=-1000)  # height/2
    TreeMesh = mt.createMesh(Tree_Geom)
    N_data = electrode_arrays.shape[0]
    # create Dataset from geometry and measured data:
    with open(OUTPUT_DIR + 'temp_dataset.dat', 'w') as f: #
        f.write('%d' % nel + '\n') #
        f.write('# ' + 'x ' + 'y ' + 'z ' + '\n') #
        for i in range(nel): #
            f.write('%.5f ' % electrode_geometry[i, 0] + '%.5f ' % electrode_geometry[i, 1] + '%.5f ' % zel + '\n') #
        # write data
        f.write('%d' % N_data + '\n') #
        f.write('# ' + 'a ' + 'b ' + 'm ' + 'n '+ '\n') #
        for i in range(N_data): # ABMNIposUposRposInegUnegRnegStufe
            A = electrode_arrays[i, 0]
            B = electrode_arrays[i, 1]
            M = electrode_arrays[i, 2]
            N = electrode_arrays[i, 3]
            f.write('%d ' % A + '%d ' % B + '%d ' % M + '%d ' % N + '\n') #

        f.write('%d' % 0) #
    f.close()

    DataSet = ert.load(OUTPUT_DIR + 'temp_dataset.dat')

    # TreeMesh.exportVTK('TreeMesh')

    Homogeneous = ert.simulate(TreeMesh,
                               res=1,
                               scheme=DataSet,
                               sr=False,
                               calcOnly=True,
                               verbose=True)

    DataSet['k'] = 1 * Homogeneous('i') / Homogeneous('u')
    DataSet['i'] = np.zeros(DataSet['a'].shape)
    DataSet['u'] = np.zeros(DataSet['a'].shape)
    DataSet['err'] = np.zeros(DataSet['a'].shape)
    DataSet['rhoa'] = np.zeros(DataSet['a'].shape)
    dataset_dict[tree] = DataSet
    DataSet.save(DATA_DIR + 'Tree_' + str(tree) + '_geometry_layout.dat')
    TreeMesh.exportVTK(DATA_DIR + 'Tree_' + str(tree) + '_mesh.vtk')

save_pickle_obj(tree_dict, directory=DATA_DIR, name='electrode_geometry')

fig, axs = plt.subplots(2,3)
fig.set_figheight(15)
fig.set_figwidth(15)
axs = axs.ravel()
for index, tree in enumerate(trees):
    DataSet = dataset_dict[tree]
    pg.show(DataSet,
            vals=DataSet['k'],
            circular=True,
            cMin=np.min(DataSet['k']),
            cMax=np.max(DataSet['k']),
            ax=axs[index],
            label='Geometric factor')
    axs[index].set_title('Tree ' + str(tree))
fig.savefig(OUTPUT_DIR + 'all_geometric_factors.png')
