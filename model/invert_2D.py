import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import pygimli.frameworks as pgf
import pyvista as pv
from core.helpers import (create_tree_mesh,
                          error_calc,
                          load_datasets,
                          create_starting_model_3D,
                          create_starting_model_2D,
                          get_filepaths)


tree = 'Tree835_1_220511_morning'
hardwood_radius = 0.0891

# Load all data and geometry into a dictionary
geometry_dict = load_datasets(experiment_plot='Irrigation',
                              tree_number=835,
                              daytime='Morning',
                              date=220510)

# Get all the filepaths for a certain experiment (e.g irrigation tree 835 during the morning of 220510)
filepaths = get_filepaths(experiment_plot='Irrigation',
                          tree_number=835,
                          daytime='Morning',
                          date=220510)
geometry_dict.keys()
one_dataset = geometry_dict['measured_values_835_morning_220510_2']

DataSet = ert.load(one_dataset['data path'])

tree_mesh = create_tree_mesh(geometry_dict['measured_values_835_morning_220510_2'])

pg.show(tree_mesh)


Homogeneous = ert.simulate(tree_mesh, res=180, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
k = 1.0 * Homogeneous('i') / Homogeneous('u')
DataSet.set('k', -k)

DataSet['err'] = np.ones(DataSet['err'].shape[0]) * 0.1

error_calc(filepaths['data'][0],
           filepaths['data'][1],
           filepaths['data'][2])

# Still need to do everything from here downwards

# Set different resistivities within the two cylinders (inner and outer one) and noise to the data
# res = [[1, 180.0],
# [2, 1340.0]]  # This maps markers 1 (outer cylinder) to 10 Ohm m and markers 2 (inner cylinder) to 100 Ohm m
# resistivity values from Bieker and Rust, 2010: Non-Destructive Estimation of Sapwood and Heartwood Width in Scots Pine (Pinus sylvestris L.)
# the error (err) is imported from the error calculation file
err = error_calc('220511_Tree_nr_835_1.txt', '220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')
DataSet.set('err', err)
DataSet.save('Tree_Dataset_GF_morning_220511_1.dat', 'a b m n err rhoa k u i')

# Circular pseudosections:
pg.show(DataSet, circular=True)
#
# ##################################################################################################################################################
#
# 2D Inversion

# load the starting model, and the homogeneous and heterogeneous mesh from the functions file
starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(DataSet, hardwood_radius)

DataSet['u'] = np.abs(DataSet['u'])
DataSet['rhoa'] = np.abs(DataSet['rhoa'])
DataSet['i'] = np.abs(DataSet['i'])

# perform the inverion for a heterogeneous (with hardwood/sapwood boundary) and homogeneous tree
ert_2D_het = ert.ERTManager()
inv_2D_het = ert_2D_het.invert(DataSet,
                               mesh=mesh_2D,
                               startModel=starting_model_2D,
                               verbose=True,
                               calcOnly=True,
                               sr=False)

sim_2D_het = ert_2D_het.simulate(scheme=DataSet,
                                 mesh=mesh_2D,
                                 res=180,
                                 sr=False,
                                 calcOnly=True,
                                 verbose=True)

cov_2D_het = ert_2D_het.coverage()
ert_2D_het.showResult(model=inv_2D_het, coverage=cov_2D_het, cMin=60, cMax=110)  # , logScale=True)
ert_2D_het.showResult(cMin=60, cMax=110, logScale=True)
#
#
#
#
# ert_2D_hom = ert.ERTManager()
# inv_2D_hom = ert_2D_hom.invert(DataSet, mesh=mesh_2D_hom, zWeight=1, startModel=180)
# cov_2D_hom = ert_2D_hom.coverage()
# ert_2D_hom.showResult(model=inv_2D_hom, coverage=cov_2D_hom, cMin=60, cMax=110)#, logScale=True)
# ert_2D_hom.showResult(cMin=60, cMax=110, logScale=True)
#
#
#
# #################################################################################################################################################
#
# # 3D Inversion
# #load the starting model, and the homogeneous and heterogeneous mesh from the functions file
# starting_model_3D,mesh_3D, mesh_3D_hom = create_starting_model_3D(hardwood_radius, height, nel, center_pos, Tree_Geom)
#
# #perform the inverion for a heterogeneous (with hardwood/sapwood boundary) and homogeneous tree
#
# ert_3D_het = ert.ERTManager()
# result_3D_het = ert_3D_het.invert(DataSet, mesh=mesh_3D, startModel=starting_model_3D)
# cov_3D_het = ert_3D_het.coverage()
# #ert_3D_het.saveResult("220511_835_morning_1_het")
# #mesh_3D.exportVTK('Tree_Nr_835_het_220511_morning_1', result_3D_het)
# # Interpolate back to 2D mesh for visualization
# interpolated_het = pg.interpolate(mesh_3D, result_3D_het, mesh_2D.cellCenters())
# ert_2D_het.showResult(interpolated_het, coverage=cov_2D_het, cMin=70, cMax=110)
#
#
#
# ert_3D_hom = ert.ERTManager()
# result_3D_hom = ert_3D_hom.invert(DataSet, mesh=mesh_3D_hom, zWeight=1, startModel=180) #lam=20, Verbose=True)#, startModel=starting_model_3D)
# cov_3D_hom = ert_3D_hom.coverage()
# #ert_3D_hom.saveResult("220511_835_morning_1_hom")
# #mesh_3D_hom.exportVTK('Tree_Nr_835_hom_220511_morning_1', result_3D_hom)
# # Interpolate back to 2D mesh for visualization
# interpolated_hom = pg.interpolate(mesh_3D_hom, result_3D_hom, mesh_2D_hom.cellCenters())
# ert_2D_hom.showResult(interpolated_hom, coverage=cov_2D_hom, cMin=70, cMax=100)
#
#
# plt.plot(ert_2D_het.inv.response, label='response 2D het')
# plt.plot(ert_2D_het.inv.dataVals, label='dataVals 2D het')
# plt.plot(ert_2D_hom.inv.response, label='response 2D hom')
# plt.plot(ert_2D_hom.inv.dataVals, label='dataVals 2D hom')
# plt.plot(ert_3D_het.inv.response, label='response 3D het')
# plt.plot(ert_3D_het.inv.dataVals, label='dataVals 3D het')
# plt.legend()
# plt.show()
#
# chi2_2D_het = ert_2D_het.inv.chi2History
# chi2_2D_hom = ert_2D_hom.inv.chi2History
# chi2_3D_het = ert_3D_het.inv.chi2History
# chi2_3D_hom = ert_3D_hom.inv.chi2History
#
# plt.plot(chi2_2D_het, label='chi^2 2D het')
# plt.plot(chi2_2D_hom, label='chi^2 2D hom')
# plt.plot(chi2_3D_het, label='chi^2 3D het')
# plt.plot(chi2_3D_hom, label='chi^2 3D hom')
# plt.legend()
# plt.show()
#
#
# # #3D visualization by Pyvista
# #cov_3D_het/=cov_3D_het.max()
# pv_mesh = pvpg.pgMesh2pvMesh(mesh_3D)
# pv_mesh['Resistivity [Ohmm]']=result_3D_het
# pv_mesh['opacity']=cov_3D_het
# pv_mesh['opacity']/=pv_mesh['opacity'].max()
# p = pv.Plotter()
# p.add_text('Title')
# p.add_mesh(
#     pv_mesh,
#     clim=[70, 110],
#     below_color='blue',
#     above_color='red',
#     scalars='Resistivity [Ohmm]',
#     opacity='opacity',
#     use_transparency=True,
#     cmap='bwr',
# )
# p.show()
# p.add_scalar_bar('Resistivity [Ohmm]',title_font_size=20,
#     label_font_size=16,
#     shadow=True,
#     n_labels=3,
#     italic=True,
#     fmt="%.1f",
#     font_family="arial",)
# p.show()
#
#
#
#
#     # #Jacobian Matrix: Sensitivtiy:
#     # fop = ert.ERTModelling()
#     # fop.setData(DataSet)
#     # fop.setMesh(mesh_2D)
#     # #model = starting_model
#     # fop.createJacobian(starting_model)
#     # sensitivities = np.sum(fop.jacobian(), axis=0)/mesh_2D.cellSizes() #adding rows: give a weigth to different configs -> no different weight so we can add it
#     # sensitivities /= np.max(np.abs(sensitivities))
#
#
