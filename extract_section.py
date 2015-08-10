__author__ = 'ruben'
__doc__ = 'Extract and creates a npy and txt files with the spikes events in the section of interest'

import scipy.io as sio
import numpy as np
import re

files = {'rat005': '/media/bigdata/i01_maze05.005/i01_maze05_MS.005_BehavElectrData.mat',
         'rat006': '/media/bigdata/i01_maze06.002/i01_maze06_MS.002_BehavElectrData.mat'}


def get_cells(path, section=None, only_pyr=None):
    """
    Extract the spikes events from the MAT file of the HC-5 DB.
    if section is provided, then spikes are split accordingly

    :param path: Path to the processed MAT file
    :param section: Name of the section to extract (default: None)
                    Run, Wheel, or Other
    :param only_pyr: return only pyramidal cells
    :return neuron: list with the spikes for each cell and lap is section
                    provided

    """
    # TODO: implement section wheel and other
    data = sio.loadmat(path)
    clusters = np.squeeze(data['Spike']['totclu'][0, 0])
    spikes = np.squeeze(data['Spike']['res'][0, 0])
    isIntern = np.squeeze(data['Clu']['isIntern'][0, 0]) == 1
    sections = np.squeeze(data['Par']['MazeSectEnterLeft'][0, 0])
    # Separate spikes by neuron number
    neuron = list()
    for n in range(0, max(clusters)):

        if only_pyr and isIntern[n]:
            continue
        spk = spikes[clusters == n]

        if section == 'Run':
            # Get intervals of interest: Section 2 to 6 that correspond to the running sections.
            # Section 1 seems to be between the running wheel and the central arm of the Maze.
            #  It is not clear the boundary thus it is avoided
            num_laps = len(sections)
            laps = list()
            for i in range(num_laps):
                start_run, end_run = sections[i][1, 0], sum(sections[i][6:8, 0]) - 1
                idx = np.where(np.logical_and(spk >= start_run, spk < end_run))
                # save spike events aligned to the entering to sect 2.
                laps.append(spk[idx] - start_run)
            neuron.append(laps)
        else:
            neuron.append(spk)
    print '{} cells extracted'.format(len(neuron))
    print 'Loading completed'
    return neuron


def raster(cells, title=''):
    """
    Creates a raster plot of spikes events

    :param cells: list (N x 1) containing spk times

    """
    import matplotlib.pyplot as plt
    space = 0
    fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    for n in range(len(cells)):
        plt.plot(cells[n], np.ones(len(cells[n])) + space, '|')
        space += 1
    plt.xlabel('Samples')
    plt.ylabel('Cell Num.')
    plt.title(title)


def spike_train(spks, length=1000):
    """
    Conver spk events to matrix of spike trains (1's and 0's)
    :param spks: spike events
    :param length: maximum length to extract
    :return: trains: matrix of spike trains (N x length) or (N x lenght x laps)

    """
    n, l = np.shape(spks)
    trains = np.zeros([n, length, l])
    c = 0
    for cell in spks:
        lp = 0
        for lap in cell:
            inside = lap[lap < length]
            trains[c, inside, lp] = 1.
            lp += 1
        c += 1
    return trains


def smooth_spk(train, width=0.1, plot=False, normalize=False):
    """
    Gaussian convolution of spike trains, averaged across trials
    :param train: spikes trains (N x Length x Trials)
    :param width: width of the gaussian kernel
    :return: smo: smoothed trains
    """
    import scipy.ndimage.filters as fil

    ave = np.mean(train, axis=2)
    smo = list()
    for n in range(len(train)):
        y = fil.gaussian_filter1d(ave[n, :], sigma=width)
        if normalize:
            den = np.max(y) - np.min(y)
            y = (y - np.min(y)) / den if den != 0. else y
        smo.append(y)

    if plot:
        import matplotlib.pyplot as plt
        space = 0.
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        for n in range(len(train)):
            plt.plot(smo[n] + space)
            space += 1.
        plt.xlabel('Samples')
        plt.ylabel('Cell Num.')

    return np.array(smo)


y = {}
for k, source in files.iteritems():
    data = get_cells(source, only_pyr=True, section='Run')
    # sanity check: raster one lap
    raster([x[0] for x in data], title='{} Lap 0'.format(k))
    y = spike_train(data, length=int(3 * 1250))
    y_smo = smooth_spk(y, width=int(0.05 * 1250), plot=True, normalize=True)
    name_file = re.findall(r'(i0\w+.\d+)', source)[0] + '_firings.npy'
    np.save(name_file, y_smo)
