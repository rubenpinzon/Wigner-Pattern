__author__ = 'ruben'

import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import ndimage
import argparse
import numpy as np
from neuronpy.graphics import spikeplot
from neuronpy.math import kernel

parser = argparse.ArgumentParser(description='Function to plot a matlab processed hc-5 database file')
parser.add_argument('PATH', type=str, nargs='+',
                    help='Path to the mat file ##.BehavElectrData.mat')
parser.add_argument("--verbosity", help="increase output verbosity")
color_ascii = {'green': '\033[1;32m{t}\033[1;m', 'red': '\033[1;31m{t}\033[1;m'}

#  Manage mouse events in the raster plot and print neurons name
class LabelPrinter:
    def __init__(self, labels, canvas):
        self.labels = labels
        self.canvas = canvas
        self.cid = canvas.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes != self.canvas.axes:
            return
        print self.labels[int(event.ydata / 5)]


def get_psth(data, bin_size, **kwargs):
    if 'keys' in kwargs:
        keys = kwargs['keys']
    else:
        keys = data.keys()

    if len(keys) != 0:
        bin_min, bin_max = (0, np.ceil(max(data[max(data)])))
        psth = np.zeros(np.ceil(bin_max/bin_size), dtype=int)
        for k in keys:
            spike_bin = data[k]
            hist, _ = np.histogram(a=spike_bin, bins=np.ceil(bin_max/bin_size), range=(bin_min, bin_max))
            psth += hist
        return psth/(len(keys)*bin_size)
    else:
        return None


if __name__ == '__main__':
    args = parser.parse_args()
    data = sio.loadmat(args.PATH[0])

    if args.verbosity:
        for key, value in data.iteritems():
            print 'Field {v} of type {p} found'.format(v=color_ascii['green'].format(t=key), p=type(value))
            if type(value).__module__ == np.__name__:
                variables = value[0].dtype.names
                print 'Number of variables {} : {}'.format(len(variables), variables)

                # TODO: add test for the presence of the fields in the MAT file
    # Unpack data of interest from the mat file
    clusters = data['Spike'][0]['totclu'][0]
    laps = data['Laps'][0]['StartLaps'][0]
    events = data['Spike'][0]['res'][0]
    fs = 1250.0

    #  extract the spikes times for each neuron
    spikes = {}
    for neuron in range(1, max(clusters)):
        spikes['neuron {}'.format(neuron)] = events[clusters == neuron] / fs
    # max time based on the time of the last spike
    time_max = max(spikes[max(spikes)])
    #  Configure the figure's window
    fig = plt.figure(frameon=False, figsize=(17, 10), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    #  Plot the spike events using vertical marker. s is to space out the spikes vertically
    c = 0
    for key, times in spikes.iteritems():
        plt.plot(times, c * np.ones(np.shape(times)), linestyle='None', marker='|', label=key)
        c += 5

    plt.ylim((-10, c + 5))
    psth = get_psth(data=spikes, bin_size=.5)
    fig2 = plt.figure(frameon=False, figsize=(17, 2), dpi=80, facecolor='w', edgecolor='k')
    ax2 = fig2.add_axes([0, 0, 1, 1])
    ax2.axis('off')
    ax2.plot(psth, marker='o')

    canvas, = ax.plot([0], [0])  # empty line to get the axis
    LabelPrinter = LabelPrinter(spikes.keys(), canvas)

    plt.show()


    # TODO: add names of the neurons in the plot
