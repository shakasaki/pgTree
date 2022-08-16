import numpy as np
import pygimli.meshtools as mt
import pygimli as pg
import os
from core import DATA_DIR, OUTPUT_DIR
import pickle
import matplotlib.pyplot as plt


def load_data_for_tree(tree: int = None):
    all_data = load_pickle_obj(directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                               name='data_dict')
    all_keys = list(all_data.keys())
    tree_data = {}
    key_list = list()
    for key in all_keys:
        if key.split('_')[0] == str(tree):
            tree_data[key] = all_data[key]
            key_list.append(key)
    return tree_data, key_list


def create_3D_plc(electrodes: np.array = None,
                  hardwood_radius: float = None,
                  height: float = 1.0,
                  zel: float = 0,
                  area: float = 0.1,
                  quality: int = 35):
    nel = electrodes.shape[0]  # n° of electrodes
    El_geom = mt.createPolygon(electrodes,
                               isClosed=True,
                               area=area,
                               quality=quality,
                               boundaryMarker=1)  # , addNodes=refi_node)
    Tree_Geom = mt.extrude(El_geom, z=height)  # marker=2)
    Tree_Geom.translate([0, 0, -height / 2])

    for i in range(nel):
        # set the electrodes as nodes with marker -99 to the geometry
        Tree_Geom.createNode(pg.Pos(electrodes[i, 0],
                                    electrodes[i, 1], zel), marker=-99)
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

    if hardwood_radius:
        center = [electrodes[:, 0].min() + (electrodes[:, 0].max() - electrodes[:, 0].min()) / 2,
                  electrodes[:, 1].min() + (electrodes[:, 1].max() - electrodes[:, 1].min()) / 2,
                  0.0]
        hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center, marker=2)
        return Tree_Geom + hardwood_cylinder
    else:
        return Tree_Geom


def create_starting_model_3D(electrodes: np.array = None,
                             res_map: list = None,
                             hardwood_radius: float = None,
                             nel: int = None,
                             center_pos: np.array = None,
                             height: float = 1,
                             zel: float = 0):
    """
    Given the experiment details, return paths to datafiles and geometry
    Args:
        :param TreeGeom: pyGimli PLC mesh
        :param res_map: list that assigns a resistivity to a marker, e.g., [1, 100] for 100 Ohm.m to marker 1
        :param hardwood_radius: Radius of hardwood
        :param height: Height of mesh (tree trunk)
        :param nel: Number of electrodes
        :param center_pos: Center position for the mesh (tree axis)
    Returns:
        starting_model_3D: The starting parameters (resistivities)
        mesh_3D: The mesh for the above parameters
        mesh_3D_hom: The homogeneous mesh (without hardwood)
    """
    hardwood_cylinder = mt.createCylinder(hardwood_radius, height, nel, center_pos, boundaryMarker=1)
    hardwood_cylinder_marker2 = mt.createCylinder(hardwood_radius, height, nel, center_pos, marker=2)

    TreeGeom = create_3D_plc(electrodes=electrodes,
                             height=height,
                             zel=zel)

    mesh_3D = mt.createMesh(TreeGeom + hardwood_cylinder)
    bogus_mesh_3D_start = mt.createMesh(TreeGeom + hardwood_cylinder_marker2)

    # Set different resistivities for the two markers
    if res_map == None:
        res_map = [[1, 180.0], [2, 1340.0]]
    if type(res_map) == list:
        res_map = np.array(res_map)
    mesh_3D_hom = mt.createMesh(TreeGeom)
    # Create a starting model
    n_cells_3D = bogus_mesh_3D_start.cellCount()  # Number of cells in the mesh
    cell_markers_3D = bogus_mesh_3D_start.cellMarkers()  # Marker assigned to each cell
    starting_model_3D = np.ones(n_cells_3D)
    for cm in np.unique(
            cell_markers_3D):  # like this you loop only over the unique markers (e.g. 1 and 2 in this example)
        parameter_value_3D = res_map[res_map[:, 0] == cm, 1][0]
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
    return {'data': [data_path + file for file in os.listdir(data_path)],
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


# Error calculation: calculation
def strip_first_col(fname,
                    delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue


# def create_tree_geometry(geometry_xy: np.array = None,
#                          area: float = None,
#                          mesh_quality: float = None,
#                          trunk_height: float = None):
#     center_pos = [geometry_xy[:, 1].min() + (geometry_xy[:, 1].max() - geometry_xy[:, 1].min()) / 2,
#                   geometry_xy[:, 2].min() + (geometry_xy[:, 2].max() - geometry_xy[:, 2].min()) / 2, 0.0]
#     El_geom = mt.createPolygon(geometry_xy[:, 1:3],
#                                isClosed=True,
#                                area=area,
#                                quality=mesh_quality,
#                                boundaryMarker=1)
#     return center_pos, mt.extrude(El_geom, z=trunk_height)


def save_pickle_obj(obj,
                    directory: str = None,
                    name: str = None):
    """
    Saves pickle object in specified directory
    Args:
        obj: The object to save in a pickle file
        directory: Directory to store the file
        name: Name of file
    Returns:
        saved file
    """
    with open(directory + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_pickle_obj(directory: str = None,
                    name: str = None):
    """
    Load pickle object in from directory
    Args:
        directory: Directory to load the file
        name: Name of file
    Returns:
        returns the loaded file
    """
    with open(directory + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def plot_data_fit(ert_inversion,
                  electrode_scheme,
                  name,
                  directory: str = OUTPUT_DIR + 'figures' + os.sep):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, )
    fig.set_figheight(7)
    fig.set_figwidth(15)
    min_data, max_data = pg.utils.interperc(ert_inversion.inv.dataVals)
    ax1.set_title("Simulated data\n")
    pg.show(electrode_scheme,
            vals=ert_inversion.inv.dataVals,
            circular=True,
            cMin=min_data,
            cMax=max_data,
            ax=ax1)
    pg.show(electrode_scheme,
            vals=ert_inversion.inv.response,
            circular=True,
            cMin=min_data,
            cMax=max_data,
            ax=ax2)
    ax2.set_title("Model response\n")
    rho_misfit = ert_inversion.inv.response - ert_inversion.inv.dataVals
    pg.show(electrode_scheme,
            vals=rho_misfit,
            circular=True,
            cMin=np.min(rho_misfit),
            cMax=np.max(rho_misfit),
            ax=ax3)
    ax3.set_title("Misfit\n")
    plt.show()
    fig.savefig(directory + name)
    print('Figure ' + name + ' saved in ' + directory)
    return fig


def plot_chi(ert_3D,
             name,
             directory: str = OUTPUT_DIR + 'figures' + os.sep):
    fig = plt.figure()
    plt.semilogy(ert_3D.inv.chi2History, '-x')
    plt.xlabel('Number of Iterations')
    plt.ylabel(r'$\chi^2$')
    plt.show()
    fig.savefig(directory + name)
    print('Figure ' + name + ' saved in ' + directory)
    return fig


def plot_field_data(field_data,
                    name: str = 'choose_a_name',
                    only_positive: bool = True):
    data_entries = list(field_data.keys())
    # [  0     1  2  3  4    5   6      7     8     9     10   11
    # [ index, A, B, M, N, Ipos, Upos, Rpos, Ineg, Uneg, Rneg, Level]
    if only_positive:
        fig, axs = plt.subplots(1, 1)
        fig.set_figheight(15)
        fig.set_figwidth(15)
        for entry in data_entries:
            axs.plot(field_data[entry][:, 7], 'x')
        axs.set_xlabel('Measurement (-)')
        axs.set_ylabel('Apparent resistivity (Ohm meters)')
        axs.set_title('Positive run')
        fig.savefig(OUTPUT_DIR + 'figures' + os.sep + name)
        return fig
    else:
        fig, axs = plt.subplots(2, 1)
        fig.set_figheight(15)
        fig.set_figwidth(15)
        for entry in data_entries:
            axs[0].plot(field_data[entry][:, 7], 'x')
        axs[0].set_xlabel('Measurement (-)')
        axs[0].set_ylabel('Apparent resistivity (Ohm meters)')
        axs[0].set_title('Positive run')
        for entry in data_entries:
            axs[1].plot(field_data[entry][:, 10], 'x')
        axs[1].set_xlabel('Measurement (-)')
        axs[1].set_ylabel('Apparent resistivity (Ohm meters)')
        axs[1].set_title('Negative run')
        fig.savefig(OUTPUT_DIR + 'figures' + os.sep + name)
        return fig


def get_data_error(field_data,
                   use_both: bool = False):
    data_entries = list(field_data.keys())
    # [  0     1  2  3  4    5   6      7     8     9     10   11
    # [ index, A, B, M, N, Ipos, Upos, Rpos, Ineg, Uneg, Rneg, Level]
    current = list()
    voltage = list()
    resistance = list()
    if use_both:
        for entry in data_entries:
            current.append(field_data[entry][:, 5])
            current.append(field_data[entry][:, 8])
            voltage.append(field_data[entry][:, 6])
            voltage.append(field_data[entry][:, 9])
            resistance.append(field_data[entry][:, 7])
            resistance.append(field_data[entry][:, 10])
            cols = len(data_entries) * 2
    else:
        for entry in data_entries:
            current.append(field_data[entry][:, 5])
            voltage.append(field_data[entry][:, 6])
            resistance.append(field_data[entry][:, 7])
            cols = len(data_entries)

    current = np.array(current).T
    voltage = np.array(voltage).T
    resistance = np.abs(np.array(resistance)).T

    md_current = np.median(current, axis=1)
    md_voltage = np.median(voltage, axis=1)
    md_resistance = np.median(resistance, axis=1)

    mdc_tiled = np.tile(md_current, (cols, 1)).T
    mdv_tiled = np.tile(md_voltage, (cols, 1)).T
    mdr_tiled = np.tile(md_resistance, (cols, 1)).T

    err_current = np.std(current - mdc_tiled, axis=1) / md_current
    err_voltage = np.std(voltage - mdv_tiled, axis=1) / md_voltage
    err_resistance = np.std(resistance - mdr_tiled, axis=1) / md_resistance

    medians = np.array([md_current,
                        md_voltage,
                        md_resistance]).T

    errors = np.array([err_current,
                       err_voltage,
                       err_resistance]).T

    return medians, errors


def get_pg_dataset(field_data,
                   tree: int = None):
    from pygimli.physics import ert
    DataSet = ert.load(OUTPUT_DIR + 'geometry' + os.sep + 'Tree_' + str(tree) + '_geometry_layout.dat')
    medians, errors = get_data_error(field_data)
    DataSet['i'] = medians[:, 0]
    DataSet['u'] = medians[:, 1]
    DataSet['err'] = errors
    return DataSet
