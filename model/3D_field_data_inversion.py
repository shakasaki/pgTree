import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from core.helpers import (load_pickle_obj,
                          create_starting_model_3D,
                          create_3D_plc,
                          load_data_for_tree,
                          plot_field_data,
                          get_data_error,
                          get_pg_dataset,
                          plot_chi,
                          plot_data_fit)
from core import OUTPUT_DIR
import os

# Choose a tree to make the geometry and load the data
tree = 875
tree_data, data_list = load_data_for_tree(tree)
print('All available datasets are : ')
print(*data_list)
chosen_dataset = data_list[0] # change number here or type the dataset name as a string
print('Chosen dataset is ' + chosen_dataset)
dataset = list(tree_data[chosen_dataset].values())[0]

# tree parameters
hardwood_radius = 0.094884
resistivities = [180, 1800]
height = 1

# inversion parameters
apply_starting_model = True
lambda_start = 1
delta_phi = 1

# Load the saved geometries for all trees, and saved data schemes
all_electrodes = load_pickle_obj(directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                                 name='electrode_geometries')
all_data = load_pickle_obj(directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                           name='data_dict')

field_data = all_data[chosen_dataset]
# Plot the field data
plot_field_data(field_data, name='field_data_survey_ ' + chosen_dataset)

center_vals, errors = get_data_error(field_data)


# load electrode scheme
load_scheme = pg.load(OUTPUT_DIR + 'geometry' + os.sep + 'Tree_' + str(tree) + '_geometry_layout.dat')
electrodes = all_electrodes[tree]['electrode positions']
center_pos = all_electrodes[tree]['center']
nel = electrodes.shape[0]

DataSet = get_pg_dataset(field_data, tree=tree)
###########################################################################################################
## Force data to negative or only remove positive u

DataSet.remove(DataSet['u'] > 0)
# voltage = DataSet['u'].array()
# voltage[voltage > 0] = -voltage[voltage > 0]
# DataSet['u'] = voltage


ert_3D = ert.ERTManager(sr=False)

if apply_starting_model:
    starting_model_3D, mesh_3D, mesh_3D_hom = create_starting_model_3D(electrodes=electrodes,
                                                                       hardwood_radius=hardwood_radius,
                                                                       height=height,
                                                                       nel=nel,
                                                                       center_pos=center_pos)
    inv_3D = ert_3D.invert(DataSet,
                           mesh=mesh_3D,
                           startModel=starting_model_3D,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)
    plot_data_fit(ert_3D,
                  DataSet,
                  'misfit_HeterogeneousMesh_FieldData_' + chosen_dataset)
    plot_chi(ert_3D, name='chi_HeterogeneousMesh_FieldData' + chosen_dataset)
else:
    homogeneous_tree_plc = create_3D_plc(electrodes=electrodes)
    homogeneous_tree_mesh = mt.createMesh(homogeneous_tree_plc)
    inv_3D = ert_3D.invert(DataSet,
                           mesh=homogeneous_tree_mesh,
                           startModel=180,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)
    homogeneous_tree_mesh.exportVTK(OUTPUT_DIR + 'meshes' + os.sep + 'HomogeneousMesh_FieldData', inv_3D)
    plot_data_fit(inv_3D,
                  load_scheme,
                  'misfit_HomogeneousMesh_FieldData_' + chosen_dataset)
    plot_chi(inv_3D, name='chi_HomogeneousMesh_FieldData' + chosen_dataset)

