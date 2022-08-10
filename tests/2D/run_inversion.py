import os

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import os
from functions import strip_first_col
from functions import Tree_geometry
from functions import ERT_position
from functions import create_starting_model_2D


lambda_start = 20
delta_phi = 1
apply_starting_model = False
hardwood_radius = 0.094884
resistivities = [180, 180]
flip_polarity = True
transform_resistivities = True
negative_k = True

GEOM = np.loadtxt('geometry_835.txt')
nel = GEOM.shape[0]  # n° of electrodes
height = 1.0  # height of the tree
zel = 0.  # electrodes z coordinates

DATA1 = np.loadtxt(strip_first_col('220511_Tree_nr_835_1.txt'))
DATA2 = np.loadtxt(strip_first_col('220511_Tree_nr_835_2.txt'))
DATA3 = np.loadtxt(strip_first_col('220511_Tree_nr_835_3.txt'))

test_arrays = np.sum(np.abs(DATA1[:, 0:4] - DATA2[:, 0:4])) + np.sum(
    np.abs(DATA1[:, 0:4] - DATA3[:, 0:4]))  # should be zero because these are the arrays and should be he same



electrode_arrays = DATA1[:, 0:4]
I_array = np.array([DATA1[:, 4], DATA2[:, 4], DATA3[:, 4]]).T
U_array = np.array([DATA1[:, 5], DATA2[:, 5], DATA3[:, 5]]).T
I_center = np.median(I_array, axis=1)
U_center = np.median(U_array, axis=1)

R = np.abs((U_array / I_array).squeeze())

R_median = np.median(R, axis=1) # this collapses from 3 columns to 1
R_center = np.tile(R_median, (3, 1)).T # create another matrix from the median
R_std = np.std(R - R_center, axis=1)
error = R_std / R_median

datapoint = np.arange(R.shape[0])
fig = plt.figure()
plt.errorbar(datapoint,
             R_median,
             yerr=error * R_median,
             elinewidth=3, linewidth=0)
plt.scatter(datapoint,
            R_median,
            c='r')
plt.xlabel('Datapoint (-)')
plt.ylabel('Measured resistance with  errorbar (Ohm)')
plt.title('Repeatability of measured positive values for 3 runs')
plt.show()
fig.savefig('errorbars.png')
N_data = electrode_arrays.shape[0]  # n° of data points

# remove dat file if already exists and create a new one
try:
    os.remove('pygimli_input_dataset.dat')
except:
    pass

# create Dataset from geometry and measured data:
with open('pygimli_input_dataset.dat', 'w') as f:  #
    f.write('%d' % nel + '\n')  #
    f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
    for i in range(nel):  #
        f.write('%.5f ' % GEOM[i, 1] + '%.5f ' % GEOM[i, 2] + '%.5f ' % zel + '\n')  #
    # write data
    f.write('%d' % N_data + '\n')  #
    f.write('# ' + 'a ' + 'b ' + 'm ' + 'n ' + 'u ' + 'i ' + '\n')  #
    for i in range(N_data):  # ABMNIposUposRposInegUnegRnegStufe
        A = electrode_arrays[i, 0]
        B = electrode_arrays[i, 1]
        M = electrode_arrays[i, 2]
        N = electrode_arrays[i, 3]
        I_pos = I_center[i]
        U_pos = U_center[i]
        if flip_polarity:
            if U_pos > 0:
                f.write('%d ' % A + '%d ' % B + '%d ' % M + '%d ' % N + '%.6f ' % U_pos + '%.6f ' % I_pos + '\n')  #
            else:  # if potential is negative switch polarity and flip M and N !!
                U_pos = -U_pos
                f.write('%d ' % A + '%d ' % B + '%d ' % N + '%d ' % M + '%.6f ' % U_pos + '%.6f ' % I_pos + '\n')  #
        else:
            f.write('%d ' % A + '%d ' % B + '%d ' % M + '%d ' % N + '%.6f ' % U_pos + '%.6f ' % I_pos + '\n')  #

    f.write('%d' % 0)  #
    f.close()


# Plot the data and extracted statistics
current_factor = 1000
fig, axs = plt.subplots(1,2)
axs[0].scatter(DATA1[:, 4]*current_factor, DATA1[:, 5])
axs[0].scatter(DATA2[:, 4]*current_factor, DATA2[:, 5])
axs[0].scatter(DATA3[:, 4]*current_factor, DATA3[:, 5])
axs[1].scatter(DATA1[:, 7]*current_factor, DATA1[:, 8])
axs[1].scatter(DATA2[:, 7]*current_factor, DATA2[:, 8])
axs[1].scatter(DATA3[:, 7]*current_factor, DATA3[:, 8])

axs[0].set_xlim([0, 0.0006*current_factor])
axs[1].set_xlim([0, 0.0006*current_factor])
axs[0].set_xlabel('Current (mA)')
axs[1].set_xlabel('Current (mA)')
axs[0].set_ylabel('Voltage (V)')
#axs[1].set_ylabel('Voltage (V)')
axs[0].set_title('Only positive runs')
axs[1].set_title('Only negative runs')

plt.show()
fig.savefig('raw_data.png')






area = 0.1  # maximum cell size for resulting triangles after mesh generation
quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
center_pos, Tree_Geom = Tree_geometry(GEOM, area, quality, height)
DataSet = ert.load('pygimli_input_dataset.dat')
starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(center_pos, DataSet, hardwood_radius)


# In[10]:





Homogeneous = ert.simulate(mesh_2D_hom, res=1, scheme=DataSet, sr=False, calcOnly=True, verbose=True)

k = 1 * Homogeneous('i') / Homogeneous('u')
DataSet["k"] = k

if negative_k:
    DataSet["rhoa"] = (DataSet['u'] / DataSet['i']) * -k
else:
    DataSet["rhoa"] = (DataSet['u'] / DataSet['i']) * k

DataSet.set('err', error)
DataSet.remove(DataSet["rhoa"] < 20)

DataSet.save('pygimli_input_dataset.dat', 'a b m n err rhoa k u i')
##################################################################################################################################################


# In[11]:


fig, ax = plt.subplots(figsize=(8,8))
pg.show(mesh_2D_hom, ax=ax)
#ax.plot(pg.x(DataSet.sensorPositions()), pg.y(DataSet.sensorPositions()), "ro")
plt.show()

# In[12]:


ert_2D_het = ert.ERTManager(sr=False)
if transform_resistivities:
    ert_2D_het.inv.dataTrans = pg.trans.Trans()


if apply_starting_model:
    inv_2D_het = ert_2D_het.invert(DataSet,
                                   mesh=mesh_2D,
                                   startModel=starting_model_2D,
                                   verbose=True,
                                   lam=lambda_start,
                                   dPhi=delta_phi)
else:
    inv_2D_het = ert_2D_het.invert(DataSet,
                                   mesh=mesh_2D_hom,
                                   verbose=True,
                                   lam=lambda_start,
                                   dPhi=delta_phi)


# In[13]:


fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(15,5))
min_model, max_model = pg.utils.interperc(ert_2D_het.model, islog=True)
pg.show(ert_2D_het.paraDomain, ert_2D_het.model, logScale=True, ax=ax1, label=pg.unit("res"), cMap="Spectral_r", cMin=min_model, cMax=max_model)
ax1.set_title("Inversion result")

min_data, max_data = pg.utils.interperc(DataSet["rhoa"])
ax2.set_title("Measured data\n")
pg.show(DataSet, vals=DataSet["rhoa"], circular=True, ax=ax2, cMin=min_data, cMax=max_data)
pg.show(DataSet, vals=ert_2D_het.inv.response, circular=True, ax=ax3, cMin=min_data, cMax=max_data)
ax3.set_title("Model response\n")
plt.show()

# In[14]:


fig, ax = plt.subplots(figsize=(5,5))
rrms = ert_2D_het.inv.relrms()
pg.show(DataSet, vals=(DataSet["rhoa"]-ert_2D_het.inv.response)/ert_2D_het.inv.response * 100, circular=True, cMap="RdBu_r", cMin=-rrms, cMax=rrms, ax=ax, label="Percentage misfit")



cov_2D_het = ert_2D_het.coverage()
fig_2D_het_cov = ert_2D_het.showResult(model=inv_2D_het,
                                       coverage=cov_2D_het,
                                       cMin=inv_2D_het.array().min(),
                                       cMax=inv_2D_het.array().max(),
                                       logScale=True)
plt.show()

# In[15]:


chi2_2D_het = ert_2D_het.inv.chi2History

plt.figure(3)
plt.plot(chi2_2D_het)
plt.xlabel('Number of Iterations')
plt.ylabel(r'$\chi^2$ error')
plt.title(r'$\chi^2$ misfit')
plt.show()



































#
#
# # 1) Crete the tree geometry (cylinder) from electrodes locations
# area = 0.1  # maximum cell size for resulting triangles after mesh generation
# quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
# # refi_node = 5  # n° of additional nodes for refinement between electrodes
# center_pos, Tree_Geom = Tree_geometry(GEOM, area, quality, height)
# Tree_Geom.translate([0, 0, -height / 2])
#
# # 2) Set the ERT position (along the circumference and equally spaced) and create mesh
# Tree, Tree_Geom = ERT_position(nel, Tree_Geom, GEOM, zel, center_pos, height)
#
# # 3) Simulate response for homogeneous resistivities distributions
# hardwood_cylinder = mt.createCylinder(hardwood_radius,
#                                       height,
#                                       nel, center_pos, marker=2)
# plcMod = Tree_Geom + hardwood_cylinder
# HetTree = mt.createMesh(plcMod)
#
# DataSet = ert.load('pygimli_input_dataset.dat')
#
# # K computation should be the same regardless of test resistivity
# resistance = DataSet['u'] / DataSet['i']
#
# if compute_k:
#     test_resistivity = 1
#     Homogeneous = ert.simulate(Tree, res=test_resistivity, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
#     DataSet["k"] = test_resistivity * Homogeneous['i'] / Homogeneous['u']
#     DataSet["rhoa"] = DataSet["k"] * resistance
#
# else:
#     DataSet["k"] = np.ones(DataSet['i'].array().shape)
#     DataSet["rhoa"] = resistance
#
# DataSet.set('err', error)
# DataSet.save('pygimli_input_dataset.dat', 'a b m n err rhoa k u i')
#
# # Circular pseudosections:
# pg.show(DataSet,
#         cMin=DataSet['rhoa'].array().min(),
#         cMax=DataSet['rhoa'].array().max(),
#         circular=True)
# plt.show()
# ##################################################################################################################################################
#
# # 2D Inversion
# # load the starting model, and the homogeneous and heterogeneous mesh from the functions file
# starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(center_pos,
#                                                                    DataSet,
#                                                                    hardwood_radius,
#                                                                    resistivities=resistivities)
#
# ert_2D_het = ert.ERTManager()
# if transform_resistivities:
#     ert_2D_het.inv.dataTrans = pg.trans.Trans()
# inv_2D_het = ert_2D_het.invert(DataSet,
#                                mesh=mesh_2D,
#                                startModel=starting_model_2D,
#                                verbose=True,
#                                lam=20,
#                                dPhi=1)
#
# # lower lamda for less accuracy, cType = regularization (1=blocky model), dPhi=1: don't abort the inversion that fast
# cov_2D_het = ert_2D_het.coverage()
# fig_2D_het_cov = ert_2D_het.showResult(model=inv_2D_het,
#                                        coverage=cov_2D_het,
#                                        cMin=inv_2D_het.array().min(),
#                                        cMax=inv_2D_het.array().max(),
#                                        logScale=True)
# plt.show()
#
# # data fit
# plt.figure()
# plt.plot(ert_2D_het.inv.response.array(), label='modeled response')
# plt.plot(ert_2D_het.inv.dataVals.array(), label='measured data')
# plt.xlabel('Data point (-)')
# plt.ylabel('Resistivity (Ohm meter)')
# plt.legend()
# plt.show()
#
# # inversion convergence
# chi2_2D_het = ert_2D_het.inv.chi2History
# plt.figure()
# plt.plot(np.arange(len(chi2_2D_het)),
#          chi2_2D_het,
#          label='chi^2 2D het',
#          marker='x')
# plt.legend()
# plt.xlabel('Number of Iterations')
# plt.ylabel('RRMS error')
# plt.title('Chi^2 misfit')
# plt.show()
