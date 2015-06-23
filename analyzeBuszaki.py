__author__ = 'ruben'

import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.ndimage.filters as sfil
import argparse
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

parser = argparse.ArgumentParser(description='Function to plot a matlab processed hc-5 database file')
parser.add_argument('PATH', type=str, nargs='+',
                    help='Path to the mat file ##.BehavElectrData.mat')
parser.add_argument("--verbosity", help="increase output verbosity")

col_ascii = {'green': '\033[1;32m{t}\033[1;m',
             'red': '\033[1;31m{t}\033[1;m',
             'blue': '\033[1;34m{t}\033[1;m'}


class LoadBuszaki(object):
    """
    Class to extract and store information from Buszaki's databases

    2015 RP
    """

    def __init__(self, path=None):
        #     Initialization various unitary variables and loading the MAT file
        self.verbose = False
        self.path = path
        assert len(self.path) > 0, 'File was not loaded'
        self.data = sio.loadmat(args.PATH[0])
        self.vars = self.get_varnames()
        if self.verbose:
            print self.vars

        # Loading data
        self.clusters = self.get_var('totclu')
        self.spikes = self.get_var('res')
        self.fs = float(self.get_var('SamplingFrequency'))
        self.spk_per_neuron = self.get_neurons()
        self.laps = np.trim_zeros(self.get_var('StartLaps'))
        self.eeg = self.get_var('eeg')
        self.time = np.linspace(start=0., stop=len(self.eeg) / self.fs, num=len(self.eeg))
        self.numN = max(self.clusters)

        self.x = self.get_var('X')
        self.y = self.get_var('Y')
        self.speed = self.get_var('speed')
        self.whspeed = self.get_var('WhlSpeedCW')
        self.section = self.get_var('MazeSection')
        self.zone = self.set_zones('whspeed')
        self.split = False
        #   gui and visuals
        self.app = QtGui.QApplication([])
        self.win = QtGui.QWidget()
        self.win.resize(1300, 900)
        self.layout = QtGui.QGridLayout()
        self.layout.setVerticalSpacing(0)
        self.win.setLayout(self.layout)
        self.rasterWidget = None
        self.lastWidget = None

    def get_varnames(self):
        var_names = []
        for key, value in self.data.iteritems():
            if type(value).__module__ == np.__name__:
                variables = value[0].dtype.names
                for v in variables:
                    var_names.append((key, v))
        return var_names

    def get_var(self, field):
        """
        Extract an existing variable from the Buszaki's MAT file

        :param field: the name of the variable
        :return: variable
        """
        for tup in self.vars:
            if field == tup[1]:
                return np.squeeze(self.data[tup[0]][0][tup[1]][0])

        print 'Variable {} not found in {}'.format(field, self.path)
        return None

    def get_neurons(self):
        spikes = {}
        for neuron in range(0, max(self.clusters)):
            spikes['neuron {}'.format(neuron)] = self.spikes[self.clusters == neuron]
        return spikes

    def raster(self, lap=0):
        """
        Plot spike raster
        :param trigger:
        :return:
        """
        pw = pg.PlotWidget()
        space = 0
        for key, times in self.spk_per_neuron.iteritems():
            if lap == 0:
                spk = times
            else:
                idx = np.where(np.logical_and(times > self.laps[lap], times <= self.laps[lap + 1]))
                spk = times[idx]
            pw.plot(spk / self.fs, np.ones(np.shape(spk)) + space, symbol='s', symbolPen=(space, self.numN), pen=None,
                    symbolSize=1)
            space += 1
        self.rasterWidget = pw
        self.layout.addWidget(pw, 0, 0)

    def get_psth(self, bin_size=0.05, lap=0):

        if lap != 0:
            bin_min, bin_max = (self.laps[lap], self.laps[lap + 1])
            psth = np.zeros(np.ceil((bin_max - bin_min) / bin_size), dtype=int)
            for key, times in self.spk_per_neuron.iteritems():
                hist, _ = np.histogram(a=times, bins=np.ceil((bin_max - bin_min) / bin_size), range=(bin_min, bin_max))
                psth += hist
            y = sfil.gaussian_filter1d(psth * self.fs / (bin_size * self.numN), sigma=100.)
            time = np.linspace(self.laps[lap] / self.fs, self.laps[lap + 1] / self.fs, num=len(y))
            self.addExtraWidget(time, y, 'b')
        else:
            return None

    def plot(self, field, lap=None, color='g', **kwargs):
        """
        Create a new figure with the specified variable If Split is active,
        data oer lap is shown on top of each other

        :param field:
        :return:
        """
        if lap:
            idx = range(self.laps[lap], self.laps[lap + 1])
            y = self.__getattribute__(field)[idx]
            time = np.linspace(self.laps[lap] / self.fs, self.laps[lap + 1] / self.fs, num=len(y))
            self.addExtraWidget(time, y, color)
        elif self.split:
            for x in range(len(self.laps) - 1):
                idx = range(self.laps[x], self.laps[x + 1])
                y = self.__getattribute__(field)[idx]
                time = np.linspace(start=0., stop=len(y) / self.fs, num=len(y))
                pw = pg.PlotWidget()
                pw.plot(time, y, pen=(x, 10))
                self.layout.addWidget(pw, x, 0)
            self.lastWidget = pw
        elif 'type' in kwargs:
            pw = self.lastWidget
            for s in self.__getattribute__(field):
                pw.plot(s / self.fs * np.array([1, 1]), np.array([0, 1]), pen=color, symbol='o')
        else:
            pw = pg.PlotWidget()
            y = self.__getattribute__(field)
            pw.plot(self.time, y, pen=(255, 0, 125))
            self.layout.addWidget(pw, self.layout.count() + 1, 0)
            self.lastWidget = pw

    def addExtraWidget(self, time, y, color):
        pw = pg.PlotWidget()

        pw.plot(time, y, pen=color)
        pw.setFixedHeight(100)
        pw.showGrid(x=True, y=True)
        if self.rasterWidget:
            pw.setXLink(self.rasterWidget)
        self.layout.addWidget(pw, self.layout.count() + 1, 0)
        self.lastWidget = pw

    def show(self):
        return self.win.show()

    def set_zones(self, field, lo=0, hi=1):
        """
        Sets the intervals to divide the laps into zones according to mazeids, whspeed, and speed
        :param intervals: tuples with sections belonging to mazeids
        :return: self.zone
        """

        return self.smith_trigger(field, hi, lo)

    def smith_trigger(self, field, hi, lo=0.):
        time_events = list()
        idx = 0
        value_old = self.__getattribute__(field)[0]

        for value in self.__getattribute__(field):
            if value >= hi > value_old:
                #   crossed thr to hi
                time_events.append(idx)
            elif value <= lo < value_old:
                #   crossed thr to lo
                time_events.append(idx)
            value_old = value
            idx += 1
        return time_events


if __name__ == '__main__':
    args = parser.parse_args()
    data = LoadBuszaki(args.PATH[0])

    #   Analyze by laps
    # data.split = True
    lap = None
    # data.raster(lap=lap)
    #data.plot('section', lap=lap, color='m')
    data.plot('speed', lap=lap)
    data.plot('whspeed', lap=lap)
    data.plot('zone', type='event')


    # data.plot('eeg', lap=lap, color='r')
    # data.get_psth(bin_size=0.1, lap=lap)

    #   Analyze by sections
    data.show()

    QtGui.QApplication.instance().exec_()
