__author__ = 'ruben'

import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import signal
import argparse
import numpy as np

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
        return psth*1000./(len(keys)*bin_size)
    else:
        return None

def get_scalogram(data, **kwargs):
    if 'wavelet' in kwargs:
        wavelet = kwargs['wavelet']
    else:
        wavelet = signal.ricker()
    if 'levels' in kwargs:
        levels = kwargs['levels']
    else:
        levels = np.arange(1, 11)
    return signal.cwt(data, wavelet, levels)

# First-order statistics
def firing_rate(spikes):
    '''
    Rate of the spike train.
    '''
    return (len(spikes)-1)/(spikes[-1]-spikes[0])

def CV(spikes):
    '''
    Coefficient of variation.
    '''
    ISI = np.diff(spikes) # interspike intervals
    return np.std(ISI)/np.mean(ISI)


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
    eeg = data['Track'][0]['eeg'][0]
    mazeId = data['Laps'][0]['MazeSection'][0]
    fs = 1250.0

    #  extract the spikes times for each neuron
    spikes = {}
    average_firing = list()
    coeff_var = list()
    for neuron in range(1, max(clusters)):
        spikes['neuron {}'.format(neuron)] = events[clusters == neuron] / fs
        average_firing.append(firing_rate(spikes['neuron {}'.format(neuron)]))
        coeff_var.append(CV(spikes['neuron {}'.format(neuron)]))
        print 'Firing rate of neuron {} = {:03.2f} Hz, and coeff. variation {:03.2f}'.format(neuron,average_firing[-1], coeff_var[-1])

    # distribution of firing rates
    time_max = max(spikes[max(spikes)])

    n, bins = np.histogram(average_firing, bins=10)
    his_fr = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    his_ax = his_fr.add_axes([0.1, 0.1, 0.8, 0.8])
    pos = np.arange(len(n))+0.5
    his_ax.barh(pos, n, align='center', height=.5, color='m')
    plt.title('Firing rate distribution')
    plt.yticks(pos-.25, ['{:03.2f} Hz'.format(b) for b in bins])

    # distribution of ISI
    n, bins = np.histogram(coeff_var, bins=10)
    his_cv = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    his_cv_ax = his_cv.add_axes([0.1, 0.1, 0.8, 0.8])
    pos = np.arange(len(n))+0.5
    his_cv_ax.barh(pos, n, align='center', height=.5, color='m')
    plt.title('Distribution of CV')
    plt.yticks(pos-.25, ['{:03.2f}'.format(b) for b in bins])

    #  Configure the figure's window
    fig = plt.figure(frameon=False, figsize=(17, 10), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_axes([0, 0.3, 1, 0.7])
    ax.axis('off')

    #  Plot the spike events using vertical marker. s is to space out the spikes vertically
    c = 0
    for key, times in spikes.iteritems():
        plt.plot(times, c * np.ones(np.shape(times)), linestyle='None', marker='|', label=key)
        c += 5

    plt.ylim((-10, c + 5))
    psth = get_psth(data=spikes, bin_size=.05)
    ax2 = fig.add_axes([0, 0.2, 1, 0.1], sharex=ax)
    ax2.axis('off')
    ax2.plot(np.arange(start=0, stop=time_max, step=time_max/len(psth)), psth)

    ax3 = fig.add_axes([0, 0, 1, 0.2], sharex=ax)
    ax3.axis('off')
    ax3.plot(np.arange(start=0, stop=time_max, step=time_max/len(eeg)), eeg, color='r')
    ax3.plot(np.arange(start=0, stop=time_max, step=time_max/len(mazeId)), (mazeId*2000)-10000., color='k')

    canvas, = ax.plot([0], [0])  # empty line to get the axis
    LabelPrinter = LabelPrinter(spikes.keys(), canvas)

    plt.show()


    # TODO: add names of the neurons in the plot
