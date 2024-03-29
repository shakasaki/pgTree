import matplotlib.pyplot as plt
from core import OUTPUT_DIR, DATA_DIR
import numpy as np
from core.helpers import save_pickle_obj
import os

try:
    os.mkdir(OUTPUT_DIR + 'data_plots')
except FileExistsError:
    print('directory exists')

output_dir = OUTPUT_DIR + 'data_plots' + os.sep
def return_directories(file_to_search):
    dirlist = []
    for filename in os.listdir(file_to_search):
        if os.path.isdir(os.path.join(file_to_search, filename)):
            dirlist.append(filename)
    return dirlist

plots = return_directories(DATA_DIR)
path_dict = {}
for plot in plots:
    trees = return_directories(DATA_DIR + plot + os.sep)
    for tree in trees:
        times = return_directories(DATA_DIR + plot + os.sep + tree + os.sep)
        for time in times:
            dates = return_directories(DATA_DIR + plot + os.sep + tree + os.sep + time + os.sep)
            for date in dates:
                full_path = DATA_DIR + plot + os.sep + tree + os.sep + time + os.sep + date + os.sep
                path_dict[tree+'_'+time+'_'+date] = [full_path + filename for filename in os.listdir(full_path)]

data_dict = {}
for survey in path_dict.keys():
    entries = path_dict[survey]
    data_dict[survey] = {}
    print('Data path : ' + survey)
    if len(entries) > 1:
        fig, axs = plt.subplots(len(entries), 1)
        for index, entry in enumerate(entries):
            data = np.loadtxt(entry)
            print('data dimensions ' + str(data.shape))
            data_dict[survey][entry] = data
            axs[index].plot(data[:, 7],'--x')
            axs[index].set_xlabel('Measurement (-)')
            axs[index].set_ylabel('Resistance (Ohms)')
        fig.savefig(output_dir + survey + '.png')
        plt.close()
    else:
        fig, axs = plt.subplots(1,1)
        data = np.loadtxt(entries[0])
        print('data dimensions ' + str(data.shape))
        data_dict[survey][entries[0]] = data
        axs.plot(data[:, 7],'--x')
        axs.set_xlabel('Measurement (-)')
        axs.set_ylabel('Resistance (Ohms)')
        fig.savefig(output_dir + survey + '.png')
        plt.close()

save_pickle_obj(data_dict, directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                name='data_dict')
save_pickle_obj(path_dict, directory=OUTPUT_DIR + 'dictionaries' + os.sep,
                name='path_dict')