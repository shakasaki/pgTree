import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from core.helpers import (load_pickle_obj,
                          create_starting_model_3D,
                          create_3D_plc)
from core import OUTPUT_DIR, DATA_DIR

lambda_start = 20
delta_phi = 1
error_percentage = 0.05
apply_starting_model = False
hardwood_radius = 0.094884
resistivities = [180, 1800]
transform_resistivities = False

# Load the saved geometries for all trees, and saved data schemes
data_dict = load_pickle_obj(directory=DATA_DIR, name='data_dict')
path_dict = load_pickle_obj(directory=DATA_DIR, name='path_dict')
all_electrodes = load_pickle_obj(directory=DATA_DIR, name='electrode_geometry')

# choose any tree to run a synthetic inversion
tree = 875
load_data = pg.load(DATA_DIR + 'Tree_' + str(tree) + '_geometry_layout.dat')
electrodes = all_electrodes[tree]['electrode positions']
center_pos = all_electrodes[tree]['center']

height = 1
nel = electrodes.shape[0]
homogeneous_tree_plc = create_3D_plc(electrodes=electrodes)
homogeneous_tree_mesh = mt.createMesh(homogeneous_tree_plc)

heterogeneous_tree_plc = create_3D_plc(electrodes=electrodes,
                                       hardwood_radius=hardwood_radius)
heterogeneous_tree_mesh = mt.createMesh(heterogeneous_tree_plc)

starting_model_3D, mesh_3D, mesh_3D_hom = create_starting_model_3D(electrodes=electrodes,
                                                                   hardwood_radius=hardwood_radius,
                                                                   height=height,
                                                                   nel=nel,
                                                                   center_pos=center_pos)

pg.show(mesh_3D, starting_model_3D)

synthetic_data = ert.simulate(mesh_3D,
                              res=starting_model_3D,
                              scheme=load_data,
                              sr=False,
                              calcOnly=True,
                              verbose=True)

synthetic_data['err'] = error_percentage * np.ones(synthetic_data['a'].shape)

ert_3D = ert.ERTManager(sr=False)
if transform_resistivities:
    ert_3D.inv.dataTrans = pg.trans.Trans()

if apply_starting_model:
    inv_3D = ert_3D.invert(synthetic_data,
                           mesh=mesh_3D,
                           startModel=starting_model_3D,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)
else:
    inv_3D = ert_3D.invert(synthetic_data,
                           mesh=homogeneous_tree_mesh,
                           startModel=180,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)

plt.plot(ert_3D.inv.response.array(), 'o')
plt.plot(synthetic_data['rhoa'].array(), 'x')
plt.show()

homogeneous_tree_mesh.exportVTK(OUTPUT_DIR + 'HomogeneousMesh_HeterogeneousSyntheticData.vtk',
                                inv_3D)

plt.plot(ert_3D.inv.dataVals.array(), 'o')
plt.plot(ert_3D.inv.response.array(), 'x')
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2)
min_data, max_data = pg.utils.interperc(ert_3D.inv.dataVals)
ax1.set_title("Measured data\n")
pg.show(load_data, vals=ert_3D.inv.dataVals, circular=True, cMin=min_data, cMax=max_data, ax=ax1)
pg.show(load_data, vals=ert_3D.inv.response, circular=True, cMin=min_data, cMax=max_data, ax=ax2)
ax2.set_title("Model response\n")
plt.show()

chi2_2D_het = ert_3D.inv.chi2History
plt.figure()
plt.semilogy(chi2_2D_het)
plt.xlabel('Number of Iterations')
plt.ylabel(r'$\chi^2$ error')
plt.title(r'$\chi^2$ misfit')
plt.show()
