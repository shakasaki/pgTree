import numpy as np
import pygimli.meshtools as mt

# Strip values:
def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue


# Error calculation: select rho from different file names
def select_val(fname, line):
    data = np.loadtxt(fname)
    val = data[:, line]
    val = np.abs(val)
    return val


# Error Calculation:
def error_calc(file1, file2, file3):
    rho1 = select_val(file1, 7)
    rho2 = select_val(file1, 10)
    rho3 = select_val(file2, 7)
    rho4 = select_val(file2, 10)
    rho5 = select_val(file3, 7)
    rho6 = select_val(file3, 10)
    # data_pos = np.array([rho2, rho4, rho6])
    # data_neg = np.array([rho1, rho3, rho5])
    all_data = np.array([rho1, rho2, rho3, rho4, rho5, rho6])
    # std = np.std(all_data, axis=0)
    med_rhoa = np.median(all_data, axis=0)
    std_2 = ((1 / 6) * rho1) ** 2 + ((1 / 6) * rho2) ** 2 + ((1 / 6) * rho3) ** 2 + ((1 / 6) * rho4) ** 2 + (
                (1 / 6) * rho5) ** 2 + ((1 / 6) * rho6) ** 2 - (med_rhoa ** 2)
    std_1 = np.sqrt(abs(std_2))
    err = (med_rhoa - std_1) / med_rhoa
    return err, med_rhoa


# Median calculation:
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
    # data_pos = np.array([rho2, rho4, rho6])
    # data_neg = np.array([rho1, rho3, rho5])
    all_data_u = np.array([u1, u2, u3, u4, u5, u6])
    all_data_i = np.array([i1, i2, i3, i4, i5, i6])
    # std = np.std(all_data, axis=0)
    med_u = np.median(all_data_u, axis=0)
    med_i = np.median(all_data_i, axis=0)
    return med_u, med_i


# load the starting model, and the homogeneous and heterogeneous mesh from the functions file

def create_starting_model_2D(DataSet, hardwood_radius,
                             sapwood_rho: float = 180,
                             hardwood_rho: float = 1340):
    geom_2D = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)
    geom_2D_start = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)
    geom_2D_hom = mt.createPolygon(DataSet.sensorPositions(), isClosed=True)

    pos = geom_2D.xmax() / 2
    geom_2D += mt.createCircle(pos=[pos, pos], radius=hardwood_radius, nSegments=21, boundaryMarker=1)
    geom_2D_start += mt.createCircle(pos=[pos, pos], radius=hardwood_radius, nSegments=21, marker=2)

    mesh_2D_start = mt.createMesh(geom_2D_start, quality=34, area=1e-4)
    mesh_2D = mt.createMesh(geom_2D, quality=34, area=1e-4)
    mesh_2D_hom = mt.createMesh(geom_2D_hom, quality=34, area=1e-4)

    # Set different resistivities for the two markers
    res = [[1, sapwood_rho],
           [2, hardwood_rho]]

    # Create a starting model
    n_cells = mesh_2D_start.cellCount()  # Number of cells in the mesh
    cell_markers = mesh_2D_start.cellMarkers()  # Marker assigned to each cell
    prior_parameters = np.array(res)
    starting_model_2D = np.ones(n_cells)
    for cm in np.unique(cell_markers):  # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
        parameter_value = prior_parameters[prior_parameters[:, 0] == cm, 1][0]
        starting_model_2D[cell_markers == cm] = parameter_value
        print('Parameter value ' + str(parameter_value) + ' assigned to cell marker ' + str(cm))

    return starting_model_2D, mesh_2D, mesh_2D_hom