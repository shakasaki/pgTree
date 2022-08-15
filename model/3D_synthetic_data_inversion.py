import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from core.helpers import (load_pickle_obj,
                          create_starting_model_3D,
                          create_3D_plc,
                          plot_data_fit,
                          plot_chi)
from core import OUTPUT_DIR
import os

# Choose a tree to make the geometry
tree = 875

# Run simulation with added noise
add_noise = True
noise_level = 0.1

# tree parameters
hardwood_radius = 0.094884
resistivities = [180, 1800]
height = 1

# inversion parameters
lambda_start = 1
delta_phi = 1

# Load the saved geometries for all trees, and saved data schemes
all_electrodes = load_pickle_obj(directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                                 name='electrode_geometries')

# choose any tree to run a synthetic inversion
load_scheme = pg.load(OUTPUT_DIR + 'geometry' + os.sep + 'Tree_' + str(tree) + '_geometry_layout.dat')
electrodes = all_electrodes[tree]['electrode positions']
center_pos = all_electrodes[tree]['center']
nel = electrodes.shape[0]

# Create mesh and PLC
homogeneous_tree_plc = create_3D_plc(electrodes=electrodes)
homogeneous_tree_mesh = mt.createMesh(homogeneous_tree_plc)
heterogeneous_tree_plc = create_3D_plc(electrodes=electrodes,
                                       hardwood_radius=hardwood_radius)
heterogeneous_tree_mesh = mt.createMesh(heterogeneous_tree_plc)

# Starting model
starting_model_3D, mesh_3D, mesh_3D_hom = create_starting_model_3D(electrodes=electrodes,
                                                                   res_map=[[1, 200.0], [2, 2000.0]],
                                                                   hardwood_radius=hardwood_radius,
                                                                   height=height,
                                                                   nel=nel,
                                                                   center_pos=center_pos)
# Create synthetic data
synthetic_data = ert.simulate(mesh_3D,
                              res=starting_model_3D,
                              scheme=load_scheme,
                              sr=False,
                              calcOnly=True,
                              verbose=True)

synthetic_data['err'] = noise_level * np.ones(synthetic_data['a'].shape)

if add_noise:
    # Add some noise to the data - the noise will be a random number from 1% to 10% of the actual values
    fig, axs = plt.subplots(1)
    axs.plot(synthetic_data['u'].array(), 'x')
    noise_vector = np.random.standard_normal(synthetic_data['u'].shape[0]) * noise_level
    synthetic_data['u'] = synthetic_data['u'] + synthetic_data['u'] * noise_vector
    axs.plot(synthetic_data['u'].array(), 'o')
    axs.set_xlabel('Measurement (-)')
    axs.set_ylabel('Voltage (V)')
    plt.legend(['synthetic data', 'synthetic data with normally distributed noise'])
    plt.show()


# run inversion with heterogeneous starting model
ert_3D_het = ert.ERTManager(sr=False)
inv_3D_heterogeneous = ert_3D_het.invert(synthetic_data,
                                         mesh=mesh_3D,
                                         startModel=starting_model_3D,
                                         verbose=True,
                                         lam=lambda_start,
                                         dPhi=delta_phi)
mesh_3D.exportVTK(OUTPUT_DIR + 'meshes' + os.sep + 'HeterogeneousMesh_HeterogeneousSyntheticData.vtk',
                  inv_3D_heterogeneous)
# Plot and save figures
plot_data_fit(ert_3D_het,
              load_scheme,
              'misfit_HeterogeneousMesh_HeterogeneousSyntheticData')
plot_chi(ert_3D_het, name='chi_HeterogeneousMesh_HeterogeneousSyntheticData')


# run inversion with homogeneous starting model
ert_3D_hom = ert.ERTManager(sr=False)
inv_3D_hom = ert_3D_hom.invert(synthetic_data,
                               mesh=homogeneous_tree_mesh,
                               startModel=180,
                               verbose=True,
                               lam=lambda_start,
                               dPhi=delta_phi)
homogeneous_tree_mesh.exportVTK(OUTPUT_DIR + 'HomogeneousMesh_HeterogeneousSyntheticData.vtk',
                                inv_3D_hom)
# Plot and save figures
plot_data_fit(ert_3D_hom,
              load_scheme,
              'misfit_HomogeneousMesh_HeterogeneousSyntheticData')
plot_chi(ert_3D_hom, name='chi_HomogeneousMesh_HeterogeneousSyntheticData')
