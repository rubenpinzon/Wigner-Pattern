__author__ = 'ruben'

import argparse
import scipy.io as sio
import brian
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Function to calculate correlation of 2 trains of spikes')
parser.add_argument("-p", '--PATH', type=str, nargs='*',
                    help='Path to the mat file ##.BehavElectrData.mat')
parser.add_argument("-w", "--window", type=float, nargs=1,
                    help="length of the correlation window. Default 10 ms",
                    default=0.01)
parser.add_argument("-n", "--neurons", type=int, nargs=1,
                    help="Neurons to compute. Default empty: i.e., autocorrelogram",
                    default=tuple())

color_ascii = {'green': '\033[1;32m{t}\033[1;m', 'red': '\033[1;31m{t}\033[1;m'}

if __name__ == '__main__':
    args = parser.parse_args()
    data = sio.loadmat(args.PATH[0])
    window = args.window[0]
    #neu_num = args.neurons[0]
    events = data['Spike'][0]['res'][0]
    clusters = data['Spike'][0]['totclu'][0]
    fs = 1250.0
    dt = 1/fs

    #  extract the spikes times for each neuron
    fig = plt.figure(frameon=False, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_axes([0, 0.2, 1, 0.8])
    ax2 = fig.add_axes([0, 0, 1, 0.1])
    ax1.axis('off')
    ax2.axis('off')

    spikes = {}
    for neuron in range(1, max(clusters)):
        spikes['neuron {}'.format(neuron)] = events[clusters == neuron] / fs

    # max time based on the time of the last spike
    fig = plt.figure(frameon=False, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_axes([0, 0.2, 1, 0.8])
    ax1.axis('off')
    corr = np.zeros((max(clusters)-1, int(2*np.ceil(window/dt))))
    for n in range(len(spikes)):
        corr[n, :] = brian.statistics.autocorrelogram(T0=spikes['neuron {}'.format(n+1)], width=window, bin=dt)
        #  ax1.bar(left=range(int(2*np.ceil(window/dt))), height=corr[n, :])
        ax1.plot(corr[n, :])

    np.savetxt('autocorr1000ms.txt',corr, fmt='%3.4f\t')

    plt.ylim((0, 0.1))
    plt.show()