import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import time
import pygimli.frameworks as pgf

from core.helpers import (create_tree_mesh,
                          ERT_position,
                          error_calc,
                          load_datasets,
                          create_starting_model_3D,
                          get_filepaths)

import pyvista as pv

start = time.time()
tree = 'Tree835_1_220511_morning'
hardwood_radius = 0.0891

geometry_dict = load_datasets(experiment_plot='Control',
                              tree_number=616,
                              daytime='Afternoon',
                              date=220511)

filepaths = get_filepaths(experiment_plot='Control',
                          tree_number=616,
                          daytime='Afternoon',
                          date=220511)

one_dataset = geometry_dict['measured_values_616_afternoon_220511_2']
tree_geom = create_tree_mesh(one_dataset)

DataSet = ert.load(one_dataset['data path'])

Homogeneous = ert.simulate(tree_geom, res=180, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
k = 1.0 * Homogeneous('i') / Homogeneous('u')
DataSet.set('k', -k)



# Still need to do everything from here downwards

# # Set different resistivities within the two cylinders (inner and outer one) and noise to the data
# # res = [[1, 180.0],
# # [2, 1340.0]]  # This maps markers 1 (outer cylinder) to 10 Ohm m and markers 2 (inner cylinder) to 100 Ohm m
# # resistivity values from Bieker and Rust, 2010: Non-Destructive Estimation of Sapwood and Heartwood Width in Scots Pine (Pinus sylvestris L.)
# # the error (err) is imported from the error calculation file
# err = error_calc('220511_Tree_nr_835_1.txt', '220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')
# DataSet.set('err', err)
# DataSet.save('Tree_Dataset_GF_morning_220511_1.dat', 'a b m n err rhoa k u i')

#
# # Circular pseudosections:
# #pg.show(DataSet, circular=True)
#
# ##################################################################################################################################################
#
# # 2D Inversion
#
# #load the starting model, and the homogeneous and heterogeneous mesh from the functions file
# starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(DataSet, hardwood_radius)
#
#
# #perform the inverion for a heterogeneous (with hardwood/sapwood boundary) and homogeneous tree
# ert_2D_het = ert.ERTManager()
# inv_2D_het = ert_2D_het.invert(DataSet, mesh=mesh_2D, startModel=starting_model_2D)
# cov_2D_het = ert_2D_het.coverage()
# ert_2D_het.showResult(model=inv_2D_het, coverage=cov_2D_het, cMin=60, cMax=110)#, logScale=True)
# ert_2D_het.showResult(cMin=60, cMax=110, logScale=True)
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
