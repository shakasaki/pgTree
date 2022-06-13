import numpy as np
import pygimli.meshtools as mt
import pygimli as pg
import os
from core import DATA_DIR, OUTPUT_DIR


# Div. functions for Python files
# Error calculation: select rho from different file names

# Strip values:
def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue


# Create Tree Geometry:
def create_tree_geometry(GEOM, area, quality, height):
    center_pos = [GEOM[:, 1].min() + (GEOM[:, 1].max() - GEOM[:, 1].min()) / 2,
                  GEOM[:, 2].min() + (GEOM[:, 2].max() - GEOM[:, 2].min()) / 2, 0.0]
    El_geom = mt.createPolygon(GEOM[:, 1:3],
                               isClosed=True,
                               area=area,
                               quality=quality,
                               boundaryMarker=1)  # , addNodes=refi_node)
    Tree_Geom = mt.extrude(El_geom, z=height)  # marker=2)
    return center_pos, Tree_Geom


# Set ERT Position and create mesh:
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
    Tree_Geom.createNode(center_pos, marker=-999)  # center_pos
    Tree_Geom.createNode([center_pos[0], center_pos[1], center_pos[2] - 1e-3 / 2])
    # The second problem for pure Neumann domains is the non-uniqueness of the partial differential equation (there are only partial derivatives of the electric potential so an arbitrary value might be added, i.e. calibrated).
    # Add calibration node with marker -1000 where the potential is fixed , somewhere on the boundary and far from the electrodes.
    Tree_Geom.createNode([center_pos[0], center_pos[1], height / 2], marker=-1000)
    Tree = mt.createMesh(Tree_Geom)
    return Tree, Tree_Geom


# Error calculation: select rho from different file names
def select_rho(fname, line):
    data = np.loadtxt(fname)
    rho = data[:, line]
    return np.abs(rho)


# Error calculation: calculation
def error_calc(file1, file2, file3):
    rho1 = select_rho(file1, 7)
    rho2 = select_rho(file1, 10)
    rho3 = select_rho(file2, 7)
    rho4 = select_rho(file2, 10)
    rho5 = select_rho(file3, 7)
    rho6 = select_rho(file3, 10)

    # data_pos = np.array([rho2, rho4, rho6])
    # data_neg = np.array([rho1, rho3, rho5])
    all_data = np.array([rho1, rho2, rho3, rho4, rho5, rho6])

    # std = np.std(all_data, axis=0)
    med = np.median(all_data, axis=0)
    std_2 = ((1 / 6) * rho1) ** 2 + ((1 / 6) * rho2) ** 2 + ((1 / 6) * rho3) ** 2 + ((1 / 6) * rho4) ** 2 + (
            (1 / 6) * rho5) ** 2 + ((1 / 6) * rho6) ** 2 - (med ** 2)
    std_1 = np.sqrt(abs(std_2))
    err = (med - std_1) / med
    return err


# 2D starting model:

def create_starting_model_2D(DataSet, hardwood_radius):
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
    res = [[1, 180.0],
           [2, 1340.0]]

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


# 3D starting model:

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
    for cm in np.unique(
            cell_markers_3D):  # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
        parameter_value_3D = prior_parameters_3D[prior_parameters_3D[:, 0] == cm, 1][0]
        starting_model_3D[cell_markers_3D == cm] = parameter_value_3D
        print('Parameter value ' + str(parameter_value_3D) + ' assigned to cell marker ' + str(cm))

    return starting_model_3D, mesh_3D, mesh_3D_hom


def create_directory(directory: str = None):
    """
    Try to make a directory unless it exists
    Args:
        directory: Directory to load the file
    Returns:
        creates the directory if it does not exist
    """
    try:
        os.mkdir(directory)
    except FileExistsError:
        print('Output directory exists')


def get_filepaths(experiment_plot: str = None,
                  tree_number: int = None,
                  daytime: str = None,
                  date: int = None) -> dict:
    """
    Given the experiment details, return paths to datafiles and geometry
    Args:
    :param experiment_plot: can be 'Control', 'Irrigation' or 'Irrigation Stop'
        :param tree_number: can be a valid tree number, e.g. 659
        :param daytime: Can be either 'Morning' or 'Afternoon'
        :param date: Can be a valid date in format YYMMDD, e.g., 220511
    Returns:
        list with data filepaths and geometry path
    """
    main_path = DATA_DIR + experiment_plot + os.sep + str(tree_number) + os.sep
    data_path = main_path + daytime + os.sep + str(date) + os.sep
    geometry_file = main_path + 'geometry_' + str(tree_number) + '.txt'
    return {'data': os.listdir(data_path),
            'geometry': geometry_file}


def load_datasets(experiment_plot: str = None,
                  tree_number: int = None,
                  daytime: str = None,
                  date: int = None,
                  electrode_height: float = 0) -> dict:
    """
    Load a tree dataset and save it in a format for pygimli. Also export data in a python dict
    Args:
        :param experiment_plot: can be 'Control', 'Irrigation' or 'Irrigation Stop'
        :param tree_number: can be a valid tree number, e.g. 659
        :param daytime: Can be either 'Morning' or 'Afternoon'
        :param date: Can be a valid date in format YYMMDD, e.g., 220511
        :param electrode_height: Height of electrodes, usually set to 0 m relative to model
    Returns:
        creates a readable file in the output and returns data as a python dictionary
    """
    main_path = DATA_DIR + experiment_plot + os.sep + str(tree_number) + os.sep
    geometry_file = main_path + 'geometry_' + str(tree_number) + '.txt'
    GEOM = np.loadtxt(geometry_file)
    nel = GEOM.shape[0]  # n° of electrodes
    data_path = main_path + daytime + os.sep + str(date) + os.sep
    files = os.listdir(data_path)
    output_file = OUTPUT_DIR + 'pg_input_data' + os.sep
    create_directory(output_file)
    data_dict = {}
    for index, file in enumerate(files):
        if file.split('.')[-1] == 'txt':
            DATA = np.loadtxt(strip_first_col(data_path + file))
            N = DATA.shape[0]  # n° of data points
            file_out = output_file + file.split('.')[-2] + '.dat'
            # create Dataset from geometry and measured data:
            with open(file_out, 'w') as f:  #
                f.write('%d' % nel + '\n')  #
                f.write('# ' + 'x ' + 'y ' + 'z ' + '\n')  #
                for i in range(nel):  #
                    f.write('%.5f ' % GEOM[i, 1] + '%.5f ' % GEOM[i, 2] + '%.5f ' % electrode_height + '\n')  #
                # write data                                                                                                                                                                            #
                f.write('%d' % N + '\n')  #
                f.write('# ' + 'a ' + 'b ' + 'm ' + 'n ' + 'rhoa ' + 'u ' + 'i ' + 'k ' + '\n')  #
                for i in range(N):  #
                    f.write('%d ' % DATA[i, 0] + '%d ' % DATA[i, 1] + '%d ' % DATA[i, 2] + '%d ' % DATA[
                        i, 3] + '%.6f ' % np.abs(DATA[
                                                     i, 6]) + '%.6f ' % DATA[i, 5] + '%.6f ' % DATA[
                                i, 4] + '%.6f ' % 1 + '\n')  #
                f.write('%d' % 0)  #
            f.close()  #
            file_name = file_out.split('/')[-1].split('.')[-2]
            data_dict[file_name] = {}
            data_dict[file_name]['geometry'] = GEOM
            data_dict[file_name]['data'] = DATA
            data_dict[file_name]['electrode height'] = electrode_height
            data_dict[file_name]['electrodes'] = nel
            data_dict[file_name]['data path'] = file_out
    return data_dict


def create_tree_mesh(geometry_dict,
                     height: float = 0.5,
                     hardwood: bool = False,
                     hardwood_radius: float = 0.1):
    area = 0.1  # maximum cell size for resulting triangles after mesh generation
    quality = 35  # minimum angle of mesh-triangles (increasing this reduces refinement)
    # refi_node = 5  # n° of additional nodes for refinement between electrodes
    geometry = geometry_dict['geometry']
    nel = geometry_dict['electrodes']
    electrode_height = geometry_dict['electrode height']
    center_pos, Tree_Geom = create_tree_geometry(geometry, area, quality, height)
    Tree_Geom.translate([0, 0, -height / 2])

    # 2) Set the ERT position (along the circumference and equally spaced) and create mesh
    Tree, Tree_Geom = ERT_position(nel, Tree_Geom, geometry, electrode_height, center_pos, height)
    if hardwood:
        hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)
        return mt.createMesh(Tree_Geom + hardwood_cylinder)
    else:
        return mt.createMesh(Tree_Geom)
