__author__ = 'ruben'
__doc__ = 'Takes the .spk .clu .whl files of the HC databases and convert them to a unified file .mat, npy'

import os
import numpy as np
import matplotlib.pyplot as plt


def update_axis(fig):
    """Change the geometry of the existing axis in a figure to add one more subplot"""
    n = len(fig.axes)
    for i in range(n):
        fig.axes[i].change_geometry(n + 1, 1, i + 1)
    plt.subplots_adjust(hspace=0.01)
    ax = fig.add_subplot(n + 1, 1, n + 1)
    return ax


def raster(cells, title='', fig=None):
    """
    Creates a raster plot of spikes events

    :param cells: list (N x 1) containing spk times

    """
    space = 0
    if fig is None:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
    else:
        ax = update_axis(fig)

    for n in range(len(cells)):
        ax.plot(cells[n], np.ones(len(cells[n])) + space, '|')
        space += 1
    plt.xlabel('Samples')
    plt.ylabel('Cell Num.')
    plt.title(title)

    return fig


def get_files(folder_base):
    """
    Reads the .spk .clu .whl .res .eeg from hc databases from Buszaki Lab

    :param folder_base: the main folder containing the hC database files
    :return: spk_cells: [n_cells, :] spike times for each cell
     """

    pattern_folder = ['.clu', '.res']
    hc_data = dict()
    files = sorted(os.listdir(folder_base), key=lambda x: x[-1])

    for reg in pattern_folder:
        values = list()
        max_clu = 0
        for case in files:
            if case.__contains__(reg):
                data = np.loadtxt(os.path.join(folder_base, case))
                if reg == '.clu':
                    values.extend(data[1:] + max_clu)
                    max_clu = max(values)
                    print max(values), data[0], case

                else:
                    values.extend(data)
        hc_data[reg] = values

    assert len(hc_data['.clu']) == len(hc_data['.res']), "Number of .clu {} and .res {} files should be equal.".format(
        len(hc_data['.clu']), len(hc_data['.res']))
    # organize according to cluster number in ascended order
    data = [(t, c) for t, c in zip(hc_data['.res'], hc_data['.clu'])]
    data = sorted(data, key=lambda x: x[1])

    n_cells = int(max(hc_data['.clu']))
    assert n_cells != 0, 'No clusters found including clusters artifact and noise clusters (0 & 1)'
    print '{} cells found'.format(n_cells)
    spk_cells = list()
    events = 0
    # Removes clusters artifact and noise clusters (0 & 1)
    for c in range(2, n_cells + 1):
        spks = [x[0] for x in data if x[1] == c]
        spk_cells.append(spks)
        events += len(spks)

    # TODO: to extract EEG has to read the XML file to know number of channels, int16
    found = False
    eeg = list()
    for case in os.listdir(folder_base):
        if case.__contains__('.eeg'):
            found = True
            eeg.append(np.fromfile(os.path.join(folder_base, case), dtype=np.int16))
            print 'EEG file with {} samples'.format(len(eeg[0]))
            break

    if not found:
        print 'eeg file not found'

    return spk_cells, eeg


if __name__ == '__main__':
    folder_base = '/media/bigdata/hc-3/ec013.15/ec013.157'
    spikes, eeg = get_files(folder_base)
    # fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    # ax = fig.add_subplot(111)
    # ax.plot(eeg)
    fig = raster(spikes, fig=None)
    plt.show()
