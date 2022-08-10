import matplotlib.pyplot as plt
from pygimli.physics import ert
from functions import *

DataSet = ert.load('data/Tree_Dataset.dat')
DataSet.set('k', -DataSet['k'])
##################################################################################################################################################
hardwood_radius = 0.0891

# 2D Inversion
starting_model_2D, mesh_2D, mesh_2D_hom = create_starting_model_2D(DataSet, hardwood_radius)

# perform the inverion for a heterogeneous (with hardwood/sapwood boundary) and homogeneous tree
ert_2D_het = ert.ERTManager()
inv_2D_het_starting_model = ert_2D_het.invert(DataSet, mesh=mesh_2D, startModel=starting_model_2D, verbose=True)
cov_2D_het = ert_2D_het.coverage()

# During the inversion three seem to be negative values for the resistivities but
# it is not clear how these arise!

# Plot the starting model
import pygimli as pg
pg.show(mesh_2D, starting_model_2D)
plt.show()

# Plot inversion results
ert_2D_het.showResult(model=inv_2D_het_starting_model,
                      coverage=cov_2D_het,
                      cMin=inv_2D_het_starting_model.array().min(),
                      cMax=inv_2D_het_starting_model.array().max())  # , logScale=True)
plt.show()

# Plot chi squared
fig, axs = plt.subplots(2,1)
axs[0].plot(ert_2D_het.inv.chi2History)
axs[0].set_xlabel('Iteration (-)')
axs[0].set_ylabel('chi squared (ohm meter)')

# Compare apparent resistivities
axs[1].plot(DataSet['rhoa'].array(), 'x')
axs[1].plot(ert_2D_het.inv.response.array(),'o')
axs[1].set_xlabel('Datapoint (-)')
axs[1].set_ylabel('Apparent resistivity (Ohm meter) (-)')
axs[1].legend(['Measured data', 'Simulated data'])
fig.show()