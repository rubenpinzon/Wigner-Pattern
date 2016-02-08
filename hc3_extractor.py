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


def raster(cells, title='', fig=None, exclude=None):
    """
    Creates a raster plot of spikes events

    :param cells: list (N x 1) containing spk times

    """
    if exclude is None:
        exclude = [0]
    if fig is None:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
    else:
        ax = update_axis(fig)

    for c in np.unique(cells[:, 1]):
        if c not in exclude:
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
    exts = ['_spikes', '_eeg', '_posX', '_mapElec']
    data = dict()
    for e in exts:
        found = False
        file = []
        for case in os.listdir(folder_base):
            if case.__contains__(e):
                found = True
                file = np.loadtxt(os.path.join(folder_base, case))
                print '{} file with {} samples'.format(e, len(file))
                break
        data[e] = file
        if not found:
            print '{} file not found'.format(e)

    return data['_spikes'], data['_eeg'], data['_posX'], data['_mapElec']


def interneurons(folder_base):
    """ Reads the isIntern.txt file with the description of the type and region of neurons in the hc-3 db

    :param folder_base: to search for the isIntern.txt file
    :return: chn: is the channel information from the hc-3 db extracted from table cell in hc3-tables.db
    """

    folder = os.path.split(folder_base)[0]
    file = []
    for case in os.listdir(folder):
        if case.__contains__('isIntern'):
            file = np.loadtxt(os.path.join(folder, case), dtype=str)
            chn = [(int(i[0]), int(i[1]), int(i[2]), int(i[3])) for i in file if i[2] != 'ele']
        return chn

    print 'IsIntern.txt file not found at'.format(folder)
    return


def match_cell_type(map, chn_info):
    """

    :param map: map file from is a matrix displaying the correspondance between new cluster numbers
                (first column) and inital shank # (second column) and cluster # (third column) from LoadClusRes.m
    :param chn_info: is the channel information from the hc-3 db extracted from table cell in hc3-tables.db
    :return: is_inh : a boolean vector
    """
    is_inh = list()
    for i in range(len(map)):
        shank = int(map[i, 1])
        cluster = int(map[i, 2])
        cell_id = False
        cell_id = [True for c in chn_info if c[2] == shank and c[3] == cluster and c[1] == 1]
        if cell_id:
            is_inh.extend([int(map[i, 0])])
    return is_inh


def sfft(data, fs=1250., fft_size=64, overlap_fac=0.5):
    # data = a numpy array containing the signal to be processed
    # fs = a scalar which is the sampling frequency of the data

    hop_size = np.int32(np.floor(fft_size * (1 - overlap_fac)))
    pad_end_size = fft_size  # the last segment can overlap the end of the data array by no more than one window size
    n_segments = np.int32(np.ceil(len(data) / np.float32(hop_size)))
    t_max = len(data) / fs

    window = np.hanning(fft_size)  # our half cosine window
    inner_pad = np.zeros(fft_size)  # the zeros which will be used to double each segment size

    proc = np.concatenate((data, np.zeros(pad_end_size)))  # the data to process
    power_db = np.empty((n_segments, fft_size), dtype=np.float32)  # space to hold the result

    for i in xrange(n_segments):  # for each segment
        current_hop = hop_size * i  # figure out the current segment offset
        segment = proc[current_hop:current_hop + fft_size]  # get the current segment
        windowed = segment * window  # multiply by the half cosine function
        padded = np.append(windowed, inner_pad)  # add 0s to double the length of the data
        spectrum = np.fft.fft(padded) / fft_size  # take the Fourier Transform and scale by the number of samples
        autopower = np.abs(spectrum * np.conj(spectrum))  # find the autopower spectrum
        power_db[i, :] = autopower[:fft_size]  # append to the results array

    power_db = 20 * np.log10(power_db)  # scale to db
    power_db = np.clip(power_db, -40, 200)  # clip values

    time = np.linspace(0, len(eeg) / 1.25, n_segments)
    return power_db.T, time


def psth(spikes, width=100, fs=1250.):
    """

    :param spikes: vector or list with spike times in ms
    :param width: in ms
    :return: smoothed psth
    """

    import scipy.ndimage.filters as sfil
    range = max(spikes) - min(spikes)
    print range, np.ceil(range / width)
    hist, _ = np.histogram(a=spikes, bins=np.ceil(range / width), range=(min(spikes), max(spikes)))
    y = sfil.gaussian_filter1d(hist, sigma=.1)
    t = np.linspace(0, max(spikes), len(y))  # ms
    return y, t


if __name__ == '__main__':
    folder_base = '/media/bigdata/hc-3/ec013.15/ec013.156'
    chn_info = interneurons(folder_base)
    spikes, eeg, pos, map = get_files(folder_base)
    is_inh = match_cell_type(map, chn_info)

    power, time_sfft = sfft(eeg, fft_size=32)

    psth, t_psth = psth(np.array(spikes)[:, 0])

    fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    time = np.linspace(0, len(eeg) / 1.25, len(eeg))
    ax.plot(time, eeg / 50)
    ax.plot(time, pos / 4, color='r')
    fig = raster(spikes, fig=fig, exclude=is_inh)
    ax = update_axis(fig)
    # img = ax.imshow(power, cmap='jet', interpolation='nearest', aspect='auto')
    # img.set_extent([0, time_sfft[-1], 0, 256])
    ax.plot(t_psth, psth)
    plt.show()
