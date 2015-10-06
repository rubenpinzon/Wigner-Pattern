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
    ax = fig.add_subplot(n + 1, 1, n + 1, sharex=fig.axes[n - 1])
    return ax


def raster(cells, title='', fig=None):
    """
    Creates a raster plot of spikes events

    :param cells: list (N x 1) containing spk times

    """
    if fig is None:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
    else:
        ax = update_axis(fig)

    for c in np.unique(cells[:, 1]):
        spks = [x[0] for x in cells if x[1] == c]
        ax.plot(spks, np.ones(len(spks)) * c, '|')

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
    found = False
    spk_cells = []
    for case in os.listdir(folder_base):
        if case.__contains__('_spikes'):
            found = True
            spk_cells = np.loadtxt(os.path.join(folder_base, case))
            print 'Spikes file with {} samples'.format(len(spk_cells))
            break

    if not found:
        print 'Spike file not found'

    found = False
    eeg = []
    for case in os.listdir(folder_base):
        if case.__contains__('_eeg'):
            found = True
            eeg = np.loadtxt(os.path.join(folder_base, case))
            print 'EEG file with {} samples'.format(len(eeg))
            break

    if not found:
        print 'eeg file not found'

    found = False
    pos = []
    for case in os.listdir(folder_base):
        if case.__contains__('_posX'):
            found = True
            pos = np.loadtxt(os.path.join(folder_base, case))
            print 'Pos file with {} samples'.format(len(eeg))
            break

    if not found:
        print 'pos file not found'

    return spk_cells, eeg, pos


if __name__ == '__main__':
    folder_base = '/media/bigdata/hc-3/ec013.15/ec013.156'
    spikes, eeg, pos = get_files(folder_base)
    fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    time = np.linspace(0, len(eeg) / 1.25, len(eeg))
    ax.plot(time, eeg / 50)
    ax.plot(time, pos)
    fig = raster(spikes, fig=fig)
    plt.show()

# TODO: get eh inhibitory neurons out using the Map variable from LoadCluRes.m