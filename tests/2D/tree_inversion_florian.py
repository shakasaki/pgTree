#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
import time

from functions import strip_first_col
from functions import Tree_geometry
from functions import ERT_position
from functions import error_calc
from functions import create_starting_model_2D
from functions import med_u_i


start = time.time()
tree = 'Tree835_1_220511_morning'
hardwood_radius = 0.094884

################################################################################
# 0) Read geometry and data file and Write data for BERT Inversion

# TO DO BEFORE ##################################################################
# geometry.txt: change , with tab/space separation between columns              #
# measured_values.csv: save as .txt + delete first row (the one with letters)   #
################################################################################

GEOM = np.loadtxt('geometry_835.txt')
nel = GEOM.shape[0]  # n° of electrodes
height = 1.0  # height of the tree
zel = 0.  # electrodes z coordinates
DATA = np.loadtxt(strip_first_col('220511_Tree_nr_835_1.txt'))
N = DATA.shape[0]  # n° of data points
med_u, med_i = med_u_i('220511_Tree_nr_835_1.txt', '220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')
err, med_rhoa, sign1, sign2, std1, std2, mean1, mean2, std = error_calc('220511_Tree_nr_835_1.txt','220511_Tree_nr_835_2.txt', '220511_Tree_nr_835_3.txt')

#create Dataset from geometry and measured data:
with open('Tree_Dataset_835_morning_220511_1_test.dat', 'w') as f:  #
    f.write('%d' % nel + '\n')  #
    f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
    for i in range(nel):  #
        f.write('%.5f ' % GEOM[i, 1] + '%.5f ' % GEOM[i, 2] + '%.5f ' % zel + '\n')  #
    # write data                                                                                                                                                                            #
    f.write('%d' % N + '\n')  #
    f.write('# ' + 'a ' + 'b ' + 'm ' + 'n '+ 'u ' + 'i ' + '\n')  #
    for i in range(N):  #
        f.write('%d ' % DATA[i, 0] + '%d ' % DATA[i, 1] + '%d ' % DATA[i, 2] + '%d ' % DATA[i, 3] + '%.6f ' % np.abs(DATA[i, 5]) + '%.6f ' % DATA[i, 4] + '\n')  #
    f.write('%d' % 0)  #
f.close()  #
#########################################################################################################################################################################################


# In[9]:


area = 0.1  # maximum cell size for resulting triangles after mesh generation
quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
center_pos, Tree_Geom = Tree_geometry(GEOM, area, quality, height)
DataSet = ert.load('Tree_Dataset_835_morning_220511_1_test.dat')
starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(center_pos, DataSet, hardwood_radius)


# In[10]:


Homogeneous = ert.simulate(mesh_2D_hom, res=1, scheme=DataSet, sr=False, calcOnly=True, verbose=True)

k = 1 * Homogeneous('i') / Homogeneous('u')
DataSet["k"] = k

DataSet["rhoa"] = (DataSet['u'] / DataSet['i']) * -k

DataSet.set('err', pg.Vector(len(med_u), 0.05))
DataSet.remove(DataSet["rhoa"] < 20)

DataSet.save('Tree_Dataset_GF_morning_220511_1_test.dat', 'a b m n err rhoa k u i')
##################################################################################################################################################


# In[11]:


fig, ax = plt.subplots(figsize=(8,8))
pg.show(mesh_2D_hom, ax=ax)
ax.plot(pg.x(DataSet.sensorPositions()), pg.y(DataSet.sensorPositions()), "ro")


# In[12]:


ert_2D_het = ert.ERTManager(sr=False)
inv_2D_het = ert_2D_het.invert(DataSet, mesh=mesh_2D_hom, lam=20, verbose=True)


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


# In[14]:


fig, ax = plt.subplots(figsize=(5,5))
rrms = ert_2D_het.inv.relrms()
pg.show(DataSet, vals=(DataSet["rhoa"]-ert_2D_het.inv.response)/ert_2D_het.inv.response * 100, circular=True, cMap="RdBu_r", cMin=-rrms, cMax=rrms, ax=ax, label="Percentage misfit")


# In[15]:


chi2_2D_het = ert_2D_het.inv.chi2History

plt.figure(3)
plt.plot(chi2_2D_het)
plt.xlabel('Number of Iterations')
plt.ylabel(r'$\chi^2$ error')
plt.title(r'$\chi^2$ misfit')
plt.show()

