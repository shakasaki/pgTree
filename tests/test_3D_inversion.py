import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import os
from core.helpers import load_pickle_obj
from core import OUTPUT_DIR, DATA_DIR

lambda_start = 20
delta_phi = 1
apply_starting_model = False
hardwood_radius = 0.094884
resistivities = [180, 1800]
flip_polarity = False
transform_resistivities = True
negative_k = True

data_dict = load_pickle_obj(directory=OUTPUT_DIR, name='data_dict')
path_dict = load_pickle_obj(directory=OUTPUT_DIR, name='path_dict')
all_electrodes = load_pickle_obj(directory=OUTPUT_DIR, name='electrode_geometry')

data = data_dict['875_Morning_220512']
data = data['/home/alexis/git/pgTree/data/Irrigation Stop/875/Morning/220512/measured_values_875_morning_220512_1.txt']

load_data = pg.load(OUTPUT_DIR + 'Tree_' + str(875) + '_geometry_layout.dat')
electrodes = all_electrodes[875]['electrode positions']
center_pos = all_electrodes[875]['center']

nel = electrodes.shape[0]  # n° of electrodes
height = 1.0  # height of the tree
zel = 0.  # electrodes z coordinates
area = 0.1  # maximum cell size for resulting triangles after mesh generation
quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)

N_data = data.shape[0]

# create Dataset from geometry and measured data:
with open('pg_test.dat', 'w') as f:  #
    f.write('%d' % nel + '\n')  #
    f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
    for i in range(nel):  #
        f.write('%.5f ' % electrodes[i, 0] + '%.5f ' % electrodes[i, 1] + '%.5f ' % zel + '\n')  #
    # write data
    f.write('%d' % N_data + '\n')  #
    f.write('# ' + 'a ' + 'b ' + 'm ' + 'n ' + 'u ' + 'i ' + '\n')  #
    for i in range(N_data):  # ABMNIposUposRposInegUnegRnegStufe
        A = data[i, 1]
        B = data[i, 2]
        M = data[i, 3]
        N = data[i, 4]
        I_pos = data[i, 5]
        U_pos = data[i, 6]
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

ert_data = ert.load('pg_test.dat')
# Plot the data and extracted statistics
current_factor = 1000
fig, axs = plt.subplots(1)
axs.scatter(data[:, 5] * current_factor,
               data[:, 6])
axs.scatter(ert_data['i'] * current_factor,
               ert_data['u'],5,'r')
axs.set_xlabel('Current (mA)')
axs.set_ylabel('Voltage (V)')
plt.show()


if negative_k:
    ert_data['k'] = -load_data['k']
else:
    ert_data['k'] = load_data['k']

ert_data['rhoa'] = ert_data['r']*ert_data['k']
ert_data['err'] = np.ones(load_data['k'].shape)*0.2




El_geom = mt.createPolygon(electrodes,
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
    Tree_Geom.createNode(pg.Pos(electrodes[i, 0],
                                electrodes[i, 1], zel),
                         marker=-99)
    # For sufficient numerical accuracy it is generally a good idea to refine the mesh in the vicinity of the electrodes positions.
    Tree_Geom.createNode([electrodes[i, 0],
                          electrodes[i, 1],
                          zel - 1e-3 / 2])

# Always need dipole current injection since there can be no current flow out of the closed boundaries.
# Define a reference electrode position inside the PLC, with a marker -999, somewhere away from the electrodes (and refine it).
reference_electrode = [electrodes[0, 0],
                       electrodes[0, 1],
                       height / 2]
Tree_Geom.createNode(reference_electrode, marker=-999)
Tree_Geom.createNode([reference_electrode[0],
                      reference_electrode[1],
                      reference_electrode[2] - 1e-3 / 2])

# The second problem for pure Neumann domains is the non-uniqueness of the partial differential equation (there are only partial derivatives of the electric potential so an arbitrary value might be added, i.e. calibrated).
# Add calibration node with marker -1000 where the potential is fixed , somewhere on the boundary and far from the electrodes.

calibration_electrode = [electrodes[0, 0],
                         electrodes[0, 1],
                         -height / 2]
Tree_Geom.createNode(calibration_electrode,
                     marker=-1000)  # height/2
Tree_Geom.createNode([calibration_electrode[0],
                      calibration_electrode[1],
                      calibration_electrode[2] - 1e-3 / 2])

TreeMesh = mt.createMesh(Tree_Geom)

hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)
plcMod = Tree_Geom + hardwood_cylinder
HetTree = mt.createMesh(plcMod)
















#
# # 3) Simulate response for homogeneous resistivities distributions
# hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)
# plcMod = Tree_Geom + hardwood_cylinder
# HetTree = mt.createMesh(plcMod)
#
# # Tree, Tree_Geom = ERT_position(nel, Tree_Geom, GEOM, zel, center_pos, height)
# # hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos,  marker=2)
# # plcMod = Tree_Geom + hardwood_cylinder
# # HetTree = mt.createMesh(plcMod)
# DataSet = ert.load('220510_morning_pygimli_input_dataset.dat')
# # starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(center_pos, DataSet, hardwood_radius)
# starting_model_3D, mesh_3D, mesh_3D_hom = create_starting_model_3D(hardwood_radius, height, nel, center_pos, Tree_Geom)
#
# # In[10]:

#
# Homogeneous = ert.simulate(mesh_3D_hom, res=1, scheme=DataSet, sr=False, calcOnly=True, verbose=True)
#
# # pg.show(Homogeneous, circular=True)
#
# k = 1 * Homogeneous('i') / Homogeneous('u')
# DataSet["k"] = k
#
# if negative_k:
#     DataSet["rhoa"] = (DataSet['u'] / DataSet['i']) * -k
# else:
#     DataSet["rhoa"] = (DataSet['u'] / DataSet['i']) * k
#
# DataSet.set('err', error)
# DataSet.remove(DataSet["rhoa"] < 20)
#
# DataSet.save('220510_morning_pygimli_input_dataset.dat', 'a b m n err rhoa k u i')
##################################################################################################################################################


# In[11]:


# fig, ax = plt.subplots(figsize=(8,8))
# pg.show(mesh_3D, ax=ax)
# #ax.plot(pg.x(DataSet.sensorPositions()), pg.y(DataSet.sensorPositions()), "ro")
# plt.show()
# fig.savefig('220510_morning_mesh_3D_hom.png')

# In[12]:


ert_3D = ert.ERTManager(sr=False)
if transform_resistivities:
    ert_3D.inv.dataTrans = pg.trans.Trans()

if apply_starting_model:
    inv_3D = ert_3D.invert(ert_data,
                           mesh=TreeMesh,
                           startModel=starting_model_3D,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)
else:
    inv_3D = ert_3D.invert(ert_data,
                           mesh=TreeMesh,
                           verbose=True,
                           lam=lambda_start,
                           dPhi=delta_phi)

# In[13]:
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
min_model, max_model = pg.utils.interperc(ert_3D.model, islog=True)
pg.show(ert_3D.paraDomain,
        ert_3D.model,
        logScale=True,
        ax=ax1,
        label=pg.unit("res"), cMap="Spectral_r",
        cMin=min_model,
        cMax=max_model)
ax1.set_title("Inversion result")

min_data, max_data = pg.utils.interperc(ert_data["rhoa"])
ax2.set_title("Measured data\n")
pg.show(ert_data, vals=ert_data["rhoa"], circular=True, cMin=min_data, cMax=max_data)
pg.show(ert_data, vals=ert_3D.inv.response, circular=True, cMin=min_data, cMax=max_data)
ax3.set_title("Model response\n")
plt.show()
fig.savefig('220510_morning_Inversion_result_measured_data_reponse_2D_hom')

# interpolated_het = pg.interpolate(mesh_3D, inv_3D_hom, mesh_2D.cellCenters())
# ert_3D.showResult(interpolated_het, coverage=cov_3D_hom, cMin=20, cMax=1000)
# ert_3D.showResult(interpolated_het, cMin=20, cMax=1000)


# In[13]:
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
min_model, max_model = pg.utils.interperc(ert_3D.model, islog=True)
pg.show(ert_3D.paraDomain, ert_3D.model, logScale=True, ax=ax1, label=pg.unit("res"), cMap="Spectral_r", cMin=min_model,
        cMax=max_model)
ax1.set_title("Inversion result")

min_data, max_data = pg.utils.interperc(ert_data["rhoa"])
ax2.set_title("Measured data\n")
pg.show(ert_data, vals=ert_data["rhoa"], circular=True, ax=ax2, cMin=min_data, cMax=max_data)
pg.show(ert_data, vals=ert_3D.inv.response, circular=True, ax=ax3, cMin=min_data, cMax=max_data)
ax3.set_title("Model response\n")
plt.show()
fig.savefig('Inversion_result_measured_data_reponse_3D_hom')

# In[14]:
plt.figure(1)
fig, ax = plt.subplots(figsize=(5, 5))
rrms = ert_3D.inv.relrms()
pg.show(ert_data, vals=(ert_data["rhoa"] - ert_3D.inv.response) / ert_3D.inv.response * 100, circular=True, cMap="RdBu_r",
        cMin=-rrms, cMax=rrms, ax=ax, label="Percentage misfit")
fig.savefig('220510_morning_relrms.png')

plt.figure(2)
fig_2D_het_cov = ert_3D.showResult(model=inv_3D,
                                   coverage=ert_3D.coverage(),
                                   cMin=inv_3D.array().min(),
                                   cMax=inv_3D.array().max(),
                                   logScale=True)
plt.show()
fig.savefig('220510_morning_coverage.png')

# In[15]:
chi2_2D_het = ert_3D.inv.chi2History
# plt.figure(3)
plt.plot(chi2_2D_het)
plt.xlabel('Number of Iterations')
plt.ylabel(r'$\chi^2$ error')
plt.title(r'$\chi^2$ misfit')
plt.show()
fig.savefig('220510_morning_chi.png')
