import numpy as np
import pygimli.meshtools as mt
import pygimli as pg
#Div. functions for Python files
#Error calculation: select rho from different file names

#Strip values:
def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue

#Create Tree Geometry:
def Tree_geometry(GEOM, area, quality, height):
    center_pos = [GEOM[:, 1].min() + (GEOM[:, 1].max() - GEOM[:, 1].min()) / 2,
              GEOM[:, 2].min() + (GEOM[:, 2].max() - GEOM[:, 2].min()) / 2, 0.0]
    El_geom = mt.createPolygon(GEOM[:, 1:3],
                           isClosed=True,
                           area=area,
                           quality=quality,
                           boundaryMarker=1)  # , addNodes=refi_node)
    Tree_Geom = mt.extrude(El_geom, z=height)  # marker=2)
    return center_pos, Tree_Geom

#Set ERT Position and create mesh:
def ERT_position(nel, Tree_Geom, GEOM, zel, center_pos, height):
    for i in range(nel):
    # set the electrodes as nodes with marker -99 to the geometry
    # Index = Tree_Geom.findNearestNode((GEOM[i,1], GEOM[i,2], zel))
    # Tree_Geom.node(Index).setMarker(-99)
        Tree_Geom.createNode(pg.Pos(GEOM[i, 1], GEOM[i, 2], zel), marker=-99)
    # For sufficient numerical accuracy it is generally a good idea to refine the mesh in the vicinity of the electrodes positions.
        Tree_Geom.createNode([GEOM[i, 1], GEOM[i, 2], zel - 1e-3 / 2])
    # Always need dipole current injection since there can be no current flow out of the closed boundaries.
    # Define a reference electrode position inside the PLC, with a marker -999, somewhere away from the electrodes (and refine it).
    Tree_Geom.createNode(center_pos, marker=-999)
    Tree_Geom.createNode([center_pos[0], center_pos[1], center_pos[2] - 1e-3 / 2])
    # The second problem for pure Neumann domains is the non-uniqueness of the partial differential equation (there are only partial derivatives of the electric potential so an arbitrary value might be added, i.e. calibrated).
    # Add calibration node with marker -1000 where the potential is fixed , somewhere on the boundary and far from the electrodes.
    Tree_Geom.createNode([center_pos[0], center_pos[1], height/2], marker=-1000)#height/2
    Tree = mt.createMesh(Tree_Geom)
    return Tree, Tree_Geom


#Error calculation: select rho from different file names
def select_val(fname, line):
    data = np.loadtxt(fname)
    val = data[:,line]
    val = np.abs(val)
    return val

#Error calculation: calculation
def error_calc(file1, file2, file3):
    rho1 = select_val(file1, 7)
    rho2 = select_val(file1, 10)
    rho3 = select_val(file2, 7)
    rho4 = select_val(file2, 10)
    rho5 = select_val(file3, 7)
    rho6 = select_val(file3, 10)
    sign1 = np.sign(rho1) + np.sign(rho3) + np.sign(rho5)
    sign2 = np.sign(rho2) + np.sign(rho4) + np.sign(rho6)

    data_pos = np.array([rho2, rho4, rho6])
    data_neg = np.array([rho1, rho3, rho5])
    #all_data = np.array([rho1, rho2, rho3, rho4, rho5, rho6])

    std1 = np.std(data_neg, axis=0)
    std2 = np.std(data_pos, axis=0)
    mean1 = np.mean(data_neg, axis=0)
    mean2 = np.mean(data_pos, axis=0)

    #std = np.std(all_data, axis=0)
    med_rhoa = np.median(data_neg, axis=0)
    std = np.std(data_neg, axis=0)
    #std_2=((1/6)*rho1)**2+((1/6)*rho2)**2+((1/6)*rho3)**2+((1/6)*rho4)**2+((1/6)*rho5)**2+((1/6)*rho6)**2-(med_rhoa**2)
    #std_1=np.sqrt(abs(std_2))
    err = (med_rhoa-std)/med_rhoa
    #err2 = (med_rhoa-std_1)/med_rhoa
    return err, med_rhoa , sign1, sign2, std1, std2, mean1, mean2, std #, err2, std_2, std_1

#2D starting model:

def create_starting_model_2D(center_pos, DataSet,
                             hardwood_radius,
                             resistivities:list = [180, 180]):
    geom_2D = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)
    geom_2D_start = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)
    geom_2D_hom = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)

    pos = geom_2D.xmax() / 2
    geom_2D += mt.createCircle(pos=[pos, pos], radius=hardwood_radius, nSegments=21, boundaryMarker=1)
    geom_2D_start += mt.createCircle(pos=[pos, pos], radius=hardwood_radius, nSegments=21, marker=2)
    geom_2D_start.createNode([center_pos[0], center_pos[1]], marker=-999)
    geom_2D_start.createNode([center_pos[0], center_pos[1] - 1e-3 / 2])
    geom_2D_start.createNode([center_pos[0]+0.05, center_pos[1]], marker=-1000)

    geom_2D.createNode([center_pos[0], center_pos[1]], marker=-999)
    geom_2D.createNode([center_pos[0], center_pos[1] - 1e-3 / 2])
    geom_2D.createNode([center_pos[0]+0.05, center_pos[1]], marker=-1000)

    geom_2D_hom.createNode([center_pos[0], center_pos[1]], marker=-999)
    geom_2D_hom.createNode([center_pos[0], center_pos[1] - 1e-3 / 2])
    geom_2D_hom.createNode([center_pos[0]+0.05, center_pos[1]], marker=-1000)

    mesh_2D_start = mt.createMesh(geom_2D_start, quality=34, area=1e-4)
    mesh_2D = mt.createMesh(geom_2D, quality=34, area=1e-4)
    mesh_2D_hom = mt.createMesh(geom_2D_hom, quality=34, area=1e-4)

    # Set different resistivities for the two markers
    res = [[1, resistivities[0]],
       [2, resistivities[1]]]

    #Create a starting model
    n_cells = mesh_2D_start.cellCount()  # Number of cells in the mesh
    cell_markers = mesh_2D_start.cellMarkers()  # Marker assigned to each cell
    prior_parameters = np.array(res)
    starting_model_2D = np.ones(n_cells)
    for cm in np.unique(cell_markers): # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
        parameter_value = prior_parameters[prior_parameters[:, 0] == cm, 1][0]
        starting_model_2D[cell_markers == cm] = parameter_value
        print('Parameter value ' + str(parameter_value) + ' assigned to cell marker ' + str(cm))

    return starting_model_2D, mesh_2D, mesh_2D_hom

#3D starting model:

def create_starting_model_3D(hardwood_radius, height, nel, center_pos, Tree_Geom):

    hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, boundaryMarker=1)
    hardwood_cylinder_marker2 = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)

    mesh_3D = mt.createMesh(Tree_Geom + hardwood_cylinder)
    bogus_mesh_3D_start = mt.createMesh(Tree_Geom + hardwood_cylinder_marker2)

    # Set different resistivities for the two markers
    res = [[1, 180.0],
       [2, 1340.0]]

    mesh_3D_hom = mt.createMesh(Tree_Geom)

    # Create a starting model
    n_cells_3D = bogus_mesh_3D_start.cellCount()  # Number of cells in the mesh
    cell_markers_3D = bogus_mesh_3D_start.cellMarkers()  # Marker assigned to each cell
    prior_parameters_3D = np.array(res)
    starting_model_3D = np.ones(n_cells_3D)
    for cm in np.unique(cell_markers_3D):  # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
        parameter_value_3D = prior_parameters_3D[prior_parameters_3D[:, 0] == cm, 1][0]
        starting_model_3D[cell_markers_3D == cm] = parameter_value_3D
        print('Parameter value ' + str(parameter_value_3D) + ' assigned to cell marker ' + str(cm))

    return starting_model_3D, mesh_3D, mesh_3D_hom

#function to calculate the mean of the input data:
#Error calculation: select rho from different file names, also neg values
def select_val_neg(fname, line):
    data = np.loadtxt(fname)
    val_raw = data[:,line]
    #val = np.abs(val)
    return val_raw

def med_u_i_neg(file1, file2, file3):
    u1 = select_val_neg(file1, 6)
    u2 = select_val_neg(file1, 9)
    u3 = select_val_neg(file2, 6)
    u4 = select_val_neg(file2, 9)
    u5 = select_val_neg(file3, 6)
    u6 = select_val_neg(file3, 9)
    i1 = select_val_neg(file1, 5)
    i2 = select_val_neg(file1, 8)
    i3 = select_val_neg(file2, 5)
    i4 = select_val_neg(file2, 8)
    i5 = select_val_neg(file3, 5)
    i6 = select_val_neg(file3, 8)

    #data_pos = np.array([rho2, rho4, rho6])
    #data_neg = np.array([rho1, rho3, rho5])
    all_data_u = np.array([u1, u2, u3, u4, u5, u6])
    all_data_i = np.array([i1, i2, i3, i4, i5, i6])

    #std = np.std(all_data, axis=0)
    med_u_neg = np.median(all_data_u, axis=0)
    med_i_neg = np.median(all_data_i, axis=0)
    return med_u_neg, med_i_neg

def med_u_i(file1, file2, file3):
    u1 = select_val(file1, 6)
    u2 = select_val(file1, 9)
    u3 = select_val(file2, 6)
    u4 = select_val(file2, 9)
    u5 = select_val(file3, 6)
    u6 = select_val(file3, 9)
    i1 = select_val(file1, 5)
    i2 = select_val(file1, 8)
    i3 = select_val(file2, 5)
    i4 = select_val(file2, 8)
    i5 = select_val(file3, 5)
    i6 = select_val(file3, 8)

    #data_pos = np.array([rho2, rho4, rho6])
    #data_neg = np.array([rho1, rho3, rho5])
    all_data_u = np.array([u1, u2, u3, u4, u5, u6])
    all_data_i = np.array([i1, i2, i3, i4, i5, i6])

    #std = np.std(all_data, axis=0)
    med_u = np.median(all_data_u, axis=0)
    med_i = np.median(all_data_i, axis=0)
    return med_u, med_i
