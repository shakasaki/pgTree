import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
from pygimli import meshtools as mt
from pygimli.physics.petro import transFwdArchieS as ArchieTrans
from pygimli.frameworks import PetroInversionManager
from pygimli.physics import ert

###############################################################################

# illustration of start: NOT UPDATED 25.04.2022
world3 = mt.createRectangle(start=[0.0, 5.0], end=[50, -10], area=0.9)
tri = mt.createRectangle(start=[0, -5], end=[50, -10],
                         marker=2, area=1.5)
poly = mt.mergePLC([world3, tri])
poly.addRegionMarker([0.2, 0], 1)
poly.addRegionMarker([0.2, -6], 2)

###############################################################################
mMesh = pg.meshtools.createMesh(poly, q=33, smooth=[1, 10])
saturation = pg.solver.parseArgToArray([[1, 0.9], [2, 0.6]], mMesh.cellCount(), mMesh)

###############################################################################

# Two Profiles for EW and NS for Inversion
# Height difference assesed in February through levelling, horizontal distance between electrodes/geophones measured on 21.04.22

# lineNS
# world1 =  mt.createPolygon([[-10,-20],[-10,0],[0,0],[1.1,0],[2.1,0],[3.1,0],[4,0],[5.05,0],[6,0],[7,0.21],[8,0.25],[9,0.27],[9.9,0.2],[11,0.15],[11.9,0],[13,-0.02],[13.95,-0.05],[15,-0.11],[16,-0.3],[17,-0.4],[18,0.0],[19.1,0.14],[20.1,0.2],[21.20,0.34],[22,0.3],[23.30,0.3],[24.20,0.25],[25.30,0.25],[26.10,0.27],[27.10,0.28],[28.10,0.3],[29.10,0.6],[30.1,0.65],[31.25,0.69],[32.30,0.69],[33.30,0.65],[34.10,0.68],[35.05,0.63],[36,0.65],[37.25,0.68],[38.20,0.65],[39.20,0.64],[40.20,0.7],[41.20,0.85],[42.20,1.0],[43.10,1.13],[44.1,1.2],[45.20,1.22],[46,1.25],[47.20,1.27],[48.20,1.3],[49.20,1.32],[60,1.32],[60,-20]], isClosed=True, area= 0.5)

# line EW
line_EW = mt.createPolygon(
    [[-10, -20], [-10, 0], [0, 0], [1, 0], [2.0, 0], [2.9, 0], [3.9, 0], [4.9, 0], [5.9, 0.13], [6.9, 0.2], [7.9, 0.3],
     [8.9, 0.4], [9.9, 0.53], [10.9, 0.7], [11.9, 0.81], [12.9, 0.81], [14, 1.03], [15, 1.05], [16, 1.17],
     [17.05, 1.26], [18.10, 1.25], [19.10, 1.30], [20.10, 1.30], [21.10, 1.35], [22, 1.45], [23.05, 1.50],
     [24.10, 1.64], [25, 1.70], [26.10, 1.85], [27.10, 1.9], [28.10, 2.03], [29.15, 2.05], [30.20, 2.15], [31.20, 2.20],
     [32.20, 2.35], [33.20, 2.4], [34.20, 2.58], [35.20, 2.7], [36.20, 2.7], [37.30, 2.8], [38.25, 2.9], [39.20, 3.0],
     [40.05, 3.1], [41.20, 3.15], [42.15, 3.28], [43.20, 3.30], [44.20, 3.62], [45.20, 3.65], [46.15, 3.72],
     [46.90, 3.89], [47.90, 3.9], [49, 4], [60, 4], [60, -20]], isClosed=True, area=0.5)

pMesh = pg.meshtools.createMesh(line_EW, q=33.0, smooth=[1, 10])

###############################################################################

# line EW to assign resistivity
world = mt.createPolygon(
    [[-10, -20], [-10, 0], [0, 0], [1, 0], [2.0, 0], [2.9, 0], [3.9, 0], [4.9, 0], [5.9, 0.13], [6.9, 0.2], [7.9, 0.3],
     [8.9, 0.4], [9.9, 0.53], [10.9, 0.7], [11.9, 0.81], [12.9, 0.81], [14, 1.03], [15, 1.05], [16, 1.17],
     [17.05, 1.26], [18.10, 1.25], [19.10, 1.30], [20.10, 1.30], [21.10, 1.35], [22, 1.45], [23.05, 1.50],
     [24.10, 1.64], [25, 1.70], [26.10, 1.85], [27.10, 1.9], [28.10, 2.03], [29.15, 2.05], [30.20, 2.15], [31.20, 2.20],
     [32.20, 2.35], [33.20, 2.4], [34.20, 2.58], [35.20, 2.7], [36.20, 2.7], [37.30, 2.8], [38.25, 2.9], [39.20, 3.0],
     [40.05, 3.1], [41.20, 3.15], [42.15, 3.28], [43.20, 3.30], [44.20, 3.62], [45.20, 3.65], [46.15, 3.72],
     [46.90, 3.89], [47.90, 3.9], [49, 4], [60, 4], [60, -20]], isClosed=True, marker=1, area=0.5)
soillayer = mt.createPolygon(
    [[-10, -0.3], [-10, 0], [0, 0], [1, 0], [2.0, 0], [2.9, 0], [3.9, 0], [4.9, 0], [5.9, 0.13], [6.9, 0.2], [7.9, 0.3],
     [8.9, 0.4], [9.9, 0.53], [10.9, 0.7], [11.9, 0.81], [12.9, 0.81], [14, 1.03], [15, 1.05], [16, 1.17],
     [17.05, 1.26], [18.10, 1.25], [19.10, 1.30], [20.10, 1.30], [21.10, 1.35], [22, 1.45], [23.05, 1.50],
     [24.10, 1.64], [25, 1.70], [26.10, 1.85], [27.10, 1.9], [28.10, 2.03], [29.15, 2.05], [30.20, 2.15], [31.20, 2.20],
     [32.20, 2.35], [33.20, 2.4], [34.20, 2.58], [35.20, 2.7], [36.20, 2.7], [37.30, 2.8], [38.25, 2.9], [39.20, 3.0],
     [40.05, 3.1], [41.20, 3.15], [42.15, 3.28], [43.20, 3.30], [44.20, 3.62], [45.20, 3.65], [46.15, 3.72],
     [46.90, 3.89], [47.90, 3.9], [49, 4], [60, 4], [60, 3.55], [49, 3.55], [47.9, 3.45], [46.9, 3.44], [46.15, 3.27],
     [45.20, 3.20], [44.20, 3.17], [43.20, 2.85], [42.15, 2.83], [41.20, 2.70], [40.05, 2.65], [39.20, 2.55],
     [38.25, 2.45], [37.30, 2.35], [36.20, 2.25], [35.20, 2.25], [34.20, 2.13], [33.20, 2.0], [32.20, 1.95],
     [31.20, 1.8], [30.20, 1.75], [29.15, 1.65], [28.10, 1.63], [27.10, 1.5], [25, 1.30], [24.10, 1.24], [23.05, 1.10],
     [22, 1.05], [21.10, 0.95], [20.10, 0.9], [19.10, 0.9], [18.10, 0.85], [17.05, 0.86], [16, 0.87], [15, 0.75],
     [14, 0.73], [12.9, 0.51], [11.9, 0.51], [10.9, 0.4], [9.9, 0.13], [8.9, 0.1], [7.9, 0.0], [6.9, -0.1],
     [5.9, -0.17], [4.9, -0.3], [3.9, -0.3], [2.9, -0.3], [2.0, -0.3], [1.0, -0.3], [0.0, -0.3]], isClosed=True,
    marker=2, area=0.5)

# line NS to assign resistivity
# world       = mt.createPolygon([[-10,-20],[-10,0],[0,0],[1.1,0],[2.1,0],[3.1,0],[4,0],[5.05,0],[6,0],[7,0.21],[8,0.25],[9,0.27],[9.9,0.2],[11,0.15],[11.9,0],[13,-0.02],[13.95,-0.05],[15,-0.11],[16,-0.3],[17,-0.4],[18,0.0],[19.1,0.14],[20.1,0.2],[21.20,0.34],[22,0.3],[23.30,0.3],[24.20,0.25],[25.30,0.25],[26.10,0.27],[27.10,0.28],[28.10,0.3],[29.10,0.6],[30.1,0.65],[31.25,0.69],[32.30,0.69],[33.30,0.65],[34.10,0.68],[35.05,0.63],[36,0.65],[37.25,0.68],[38.20,0.65],[39.20,0.64],[40.20,0.7],[41.20,0.85],[42.20,1.0],[43.10,1.13],[44.1,1.2],[45.20,1.22],[46,1.25],[47.20,1.27],[48.20,1.3],[49.20,1.32],[60,1.32],[60,-20]], isClosed=True, marker = 1, area= 0.5)
# soillayer   = mt.createPolygon([[-10,-0.45],[-10,0],[0,0],[1.1,0],[2.1,0],[3.1,0],[4,0],[5.05,0],[6,0],[7,0.21],[8,0.25],[9,0.27],[9.9,0.2],[11,0.15],[11.9,0],[13,-0.02],[13.95,-0.05],[15,-0.11],[16,-0.3],[17,-0.4],[18,0.0],[19.1,0.14],[20.1,0.2],[21.20,0.34],[22,0.3],[23.30,0.3],[24.20,0.25],[25.30,0.25],[26.10,0.27],[27.10,0.28],[28.10,0.3],[29.10,0.6],[30.1,0.65],[31.25,0.69],[32.30,0.69],[33.30,0.65],[34.10,0.68],[35.05,0.63],[36,0.65],[37.25,0.68],[38.20,0.65],[39.20,0.64],[40.20,0.7],[41.20,0.85],[42.20,1.0],[43.10,1.13],[44.1,1.2],[45.20,1.22],[46,1.25],[47.20,1.27],[48.20,1.3],[49.20,1.32],[60,1.32],[60,1.16],[49.20,1.16],[48.20,1.14],[47.20,1.11],[46,1.09],[45.20,1.06],[44.1,1.04],[43.10,0.98],[42.20,0.85],[41.20,0.7],[40.20,0.55],[39.20,0.49],[38.20,0.50],[37.25,0.53],[36,0.50],[35.05,0.48],[34.10,0.38],[33.30,0.35],[32.30,0.39],[31.25,0.39],[30.1,0.35],[29.10,0.3],[28.10,0.0],[27.10,-0.02],[26.10,-0.03],[25.30,-0.05],[24.20,-0.05],[23.30,0.0],[22,0.0],[21.20,0.04],[20.1,-0.1],[19.1,-0.16],[18,-0.3],[17,-0.7],[16,-0.95],[15,-0.56],[13.95,-0.50],[13,-0.47],[11.9,-0.45],[11,-0.3],[9.9,-0.25],[9,-0.18],[8,-0.2],[7,-0.24],[6,-0.45],[5.05,-0.45],[4,-0.45],[3.1,-0.45],[2.1,-0.45],[1.1,-0.45],[0,-0.45]],isClosed=True, area= 0.5, marker = 2)

# merge to one "two layer world"
world2 = mt.mergePLC([world, soillayer])

# Generate the mesh for the resistivity vector creation
sMesh = pg.meshtools.createMesh(world2, q=33.0, smooth=[1, 10])

##############################################################################
# Start Model
# Set different resistivities for the two layers #2 TopSoil, #1 colluvium present below topsoil

# Resistivities in Ohm*m measured on 21.04.2022
res = [[1, 25.0],
       [2, 10.0]]

# Gives the cells in the topSoil and the colluvium the measured resistivities

n_cells = sMesh.cellCount()  # Number of cells in the mesh
cell_markers = sMesh.cellMarkers()  # Marker assigned to each cell
prior_parameters = np.array(res)
starting_model = np.ones(n_cells)
for cm in np.unique(cell_markers):  # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
    parameter_value = prior_parameters[prior_parameters[:, 0] == cm, 1][0]
    starting_model[cell_markers == cm] = parameter_value
    print('Parameter value ' + str(parameter_value) + ' assigned to cell marker ' + str(cm))

# Show the mesh with the starting model
pg.show(sMesh,
        data=starting_model,
        label=pg.unit('res'),
        showMesh=True,
        xlabel='Distance (m)',
        ylabel='Depth (m)')
plt.show()

###############################################################################


# Introduction of Archies and Wyllies law for callculation of saturation
# parameters like rFluid, phi, cementation exponent, tortuosity etc. need to be defined better: STAND 25.04.2022

ertTrans = ArchieTrans(rFluid=65, phi=0.35)
res = ertTrans(saturation)

#############################################################################

# Load ERT and seismic data, furthermore geometric factor k is added
ERT = ert.ERTManager()
ertData = pg.physics.ert.load('ert_EW_Wen_2107.txt', verbose=False)
ertData['k'] = ert.createGeometricFactors(ertData, numerical=True)

# Invert the data using the starting model
# ERT.fop.setMesh(mMesh, ignoreRegionManager=True)                               #  If not, region 1 is set to background which is not desired
# Topography effect ERT

###############################################################################

# Single inversions

pg.info("ERT Inversion")
resInv = ERT.invert(ertData,
                    mesh=pMesh,
                    zWeight=1,
                    startModel=starting_model,
                    lam=20,
                    verbose=True)
ERT.inv.echoStatus()
ERT.showResultAndFit()
pg.show(ertData, xlabel='m', )

# chi2 = ert.inv.getChi2()
# chi = pg.frameworks.inv.chi2()
# print(chi2)
pg.show(pMesh, ERT.inv.model)

###############################################################################

# Sinle petrogeophysical inversion for saturation
pg.info("ERT Petrogeophysical Inversion")
ERTPetro = PetroInversionManager(petro=ertTrans, mgr=ERT)
satERT = ERTPetro.invert(ertData, mesh=pMesh, limits=[0., 1.], lam=5,
                         verbose=False)
ERTPetro.inv.echoStatus()

# jacobian matrix

fop = ert.ERTModelling()
fop.setData(ertData)
fop.setMesh(pMesh)

model = starting_model  # np.ones(pMesh.cellCount())
fop.createJacobian(model)

sensitivity_matrix = fop.jacobian()

sensitivities = np.sum(sensitivity_matrix, axis=0)

# Log-scaled and normalized sensitivity
normsens = pg.utils.logDropTol(sensitivities / pMesh.cellSizes(), 8e-4)
normsens /= np.max(normsens)

# here you need pandas and openpyxl (latter is to read excel file)
# if you dont have them run
# pip install pandas
# pip install openpyxl
# if you do not manage to install these it is easy to run the code without them.
# but pandas is very useful (creates tables, called dataframes)

import pandas as pd
import numpy as np


def resistivity_from_tree(tree_dataframe,
                          mesh,
                          resistivities,
                          sensitivity,
                          use_sensitivity: bool = True,
                          max_distance: float = 1000,
                          column_name: str = 'test name'):
    '''you give as input a dataframe that contains the tree information (see below)
    and the pygimli mesh and inversion result. the function computes the resistance per tree
    it is resistance because it is resistivity (ohm meter) divided by distance (meter)'''
    centers = mesh.cellCenter().array()
    area = mesh.cellSizes().array()
    resistivity_array = resistivities.array()
    if use_sensitivity:
        sensitivity_array = np.abs(sensitivity.array())
    new_df = tree_dataframe.copy()
    new_df[column_name] = np.nan
    for index, tree in new_df.iterrows():  # loop over tree
        tree_position = np.array([tree['Tree Pos. (m)'], 0, 0])  # Here i assume that Tree Pos. is along x !!
        # sum over all the inversion cells
        total_res = 0
        for res_index, res in enumerate(resistivity_array):  # loop over mesh cells
            distance = np.linalg.norm(tree_position - centers[res_index, :])
            if distance < max_distance:
                if use_sensitivity:
                    total_res += (area[res_index] * sensitivity_array[res_index] * res / np.linalg.norm(
                        tree_position - centers[res_index, :]))
                else:
                    total_res += (area[res_index] * res) / np.linalg.norm(tree_position - centers[res_index, :])
        new_df.at[index, column_name] = total_res
    return new_df


# read the excel file for each line
tree_positions_NS = pd.read_csv('Tree_Positions_NS.csv')
tree_positions_EW = pd.read_csv('Tree_Positions_EW.csv')

# give the tree positions and mesh and inversion results
# function returns the table with added the resistance per tree
tree_positions_NS = resistivity_from_tree(tree_positions_NS,
                                          pMesh,
                                          ERT.inv.model,
                                          sensitivity=normsens,
                                          use_sensitivity=True,
                                          max_distance=10,
                                          column_name='prior to irrigation')

tree_positions_NS = resistivity_from_tree(tree_positions_NS,
                                          pMesh,
                                          ERT.inv.model,
                                          sensitivity=normsens,
                                          use_sensitivity=False,
                                          max_distance=10,
                                          column_name='without sensitivity')

# Plot the results
plt.figure()
plt.plot(tree_positions_NS['Tree Pos. (m)'],
         tree_positions_NS['with sensitivity'],
         'bx')
plt.plot(tree_positions_NS['Tree Pos. (m)'],
         tree_positions_NS['without sensitivity'],
         'ro')
plt.xlabel('Distance along line (m)')
plt.legend(['With sensitivity scaling', 'without'])
plt.ylabel('Resistance index (Ohms * m^2)')
plt.show()

pg.show(pMesh, normsens, cMap="RdGy_r", orientation="vertical",
        label="Normalized\nsensitivity", nLevs=3, cMin=-1, cMax=1)

cov = ERT.coverage()
ERT.showResult(coverage=cov, cMin=50, cMax=20000, xlabel='m', ylabel='m')
plt.legend()
plt.show()

#############################################################################
# plot data with transparency added
# pMesh.addData('res', resInv)
# pMesh.addData('cov', sensitivities)
# pMesh.addData( pMesh, 'pMesh_sens')

##############################################################################
# plot ERT inversion Chi2
chi2 = ERT.inv.chi2History
plt.plot(chi2, label="ERT")
plt.xlabel('Steps')
plt.ylabel('RRMS')
plt.legend()
plt.show()

# other inversion related data
plt.plot(ERT.inv.response)
plt.plot(ERT.inv.dataVals)
plt.show()

##############################################################################

ERT.showData(ertData)


def showModel(ax, model, mesh, petro=1, cMin=None, cMax=None, label=None,
              cMap=None, showMesh=False, logScale=True):
    """Utility function to show and save models for the CG paper."""
    if cMin is None:
        cMin = 0.0
    if cMax is None:
        cMax = 1.0

    if cMap is None:
        cMap = 'viridis'
    if petro:
        ax, _ = pg.show(mesh, model, label=label,
                        logScale=False, cMin=cMin, cMax=cMax, cMap=cMap, ax=ax)
    else:
        ax, _ = pg.show(mesh, model, label=label,
                        logScale=logScale, cMin=cMin, cMax=cMax, cMap=cMap, ax=ax)

    ax.xaxis.set_ticks([-5, 0, 6, 12, 18, 24, 30, 36, 42, 48, 55])
    ax.yaxis.set_ticks([-10, -8, -6, -4, -2, 0, 2, 4, 6])
    ax.set_xlabel("m")
    ax.set_ylabel("m")

    pg.viewer.mpl.drawSensors(ax, ertData.sensorPositions(), diam=0.3, facecolor='white', edgecolor='black')

    # despine(ax=ax, offset=5, trim=True)
    if showMesh:
        pg.viewer.mpl.drawSelectedMeshBoundaries(ax, mesh.boundaries(),
                                                 linewidth=0.1, color="0.2")
    return ax


##############################################################################


axs = [None] * 4

showModel(axs[0], saturation, mMesh, showMesh=True,
          label=r'Saturation (${\tt petro}$)')
showModel(axs[1], res, mMesh, petro=0, cMin=100, cMax=100000, showMesh=1,
          label=pg.unit('res'), logScale=True, cMap=pg.cmap('res'))
showModel(axs[2], resInv, pMesh, 0, cMin=100, cMax=10000,
          label=pg.unit('res'), logScale=True, cMap=pg.cmap('res'))
showModel(axs[3], satERT, pMesh,
          label=r'Saturation (${\tt satERT}$)')
