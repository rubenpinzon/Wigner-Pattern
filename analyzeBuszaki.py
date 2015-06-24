__author__ = 'ruben'

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
        self.sectionIO = self.get_var('MazeSectEnterLeft')
        self.tmax = self.get_var('SyncOff')
        self.panels = {'run': 1, 'whl': 0}
        self.panels_count = {'run': 0, 'whl': 0}


        self.x = self.get_var('X')
        self.y = self.get_var('Y')
        self.speed = self.get_var('speed')
        self.whspeed = self.get_var('WhlSpeedCW')
        self.section = self.get_var('MazeSection')
        self.split = True
        #   gui and visuals
        self.app = QtGui.QApplication([])
        self.win = QtGui.QWidget()
        self.win.resize(1300, 900)
        self.layout = QtGui.QGridLayout()
        self.layout.setVerticalSpacing(0)
        self.win.setLayout(self.layout)
        self.rasterWidget = None

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
        Need simplification
        """
        widget = {'run': pg.PlotWidget(), 'whl': pg.PlotWidget()}
        space = 0
        for key, spikes in self.spk_per_neuron.iteritems():
            if lap == 0:
                spk = spikes
            else:
                win = self.extract_setions(lap)
                for name, w in win.iteritems():
                    idx = np.where(np.logical_and(spikes > w[0], spikes <= w[1]))
                    spk = spikes[idx] - w[0]
                    widget[name].plot(spk / self.fs, np.ones(np.shape(spk)) + space, symbol='s',
                                      symbolPen=(space, self.numN), pen=None, symbolSize=1)
            space += 1
        self.rasterWidget = widget
        for key, pw in widget.iteritems():
            self.layout.addWidget(pw, 0, self.panels[key])
            self.panels_count[key] += 1

    def get_psth(self, bin_size=0.05, lap=0):

        if lap != 0:
            win = self.extract_setions(lap)
            for name, w in win.iteritems():
                bin_min, bin_max = (w[0], w[1])
                psth = np.zeros(np.ceil((bin_max - bin_min) / bin_size), dtype=int)
                for key, times in self.spk_per_neuron.iteritems():
                    hist, _ = np.histogram(a=times, bins=np.ceil((bin_max - bin_min) / bin_size), range=(bin_min, bin_max))
                    psth += hist
                y = sfil.gaussian_filter1d(psth * self.fs / (bin_size * self.numN), sigma=100.)
                time = np.linspace(0., (w[1] - w[0]) / self.fs, num=len(y))
                self.addWidget(time, y, 'c', name)
        else:
            return None

    def plot(self, field, lap=0, color='g'):
        """
        Create a new figure with the specified variable If Split is active,
        data oer lap is shown on top of each other

        :param field:
        :return:
        """
        if lap != 0:
            win = self.extract_setions(lap)
            for panel, w in win.iteritems():
                y = self.__getattribute__(field)[w[0]:w[1]]
                time = np.linspace(0, (w[1] - w[0]) / self.fs, num=len(y))
                self.addWidget(time, y, color, panel)
        else:
            pw = pg.PlotWidget()
            y = self.__getattribute__(field)
            pw.plot(self.time, y, pen=(255, 0, 125))
            self.layout.addWidget(pw, 0, 0)

    def addWidget(self, time, y, color, parent):
        pw = pg.PlotWidget()

        pw.plot(time, y, pen=color)
        pw.setFixedHeight(100)
        pw.showGrid(x=True, y=True)
        if self.rasterWidget:
            pw.setXLink(self.rasterWidget[parent])
        self.layout.addWidget(pw, self.panels_count[parent], self.panels[parent])
        self.panels_count[parent] += 1


    def show(self):
        return self.win.show()

    def extract_setions(self, lap):
        """splits the laps into two regions, maze and running
            running 1-6
            water 8-10
            wheel 11-13
        """
        enter = self.sectionIO[lap]

        intervals = {'run': (enter[0, 0], enter[5, 1] if enter[4, 1] != 0 else enter[5, 1]),
                     'whl': (enter[11, 0], enter[12, 1])}

        return intervals


if __name__ == '__main__':
    args = parser.parse_args()
    data = LoadBuszaki(args.PATH[0])

    #   Analyze by laps
    data.split = True
    lap = 5
    data.raster(lap=lap)
    data.plot('section', lap=lap, color='m')
    data.plot('whspeed', lap=lap)
    data.plot('speed', lap=lap)
    #
    # data.plot('eeg', lap=lap, color='r')
    data.get_psth(bin_size=0.05, lap=lap)

    #   Analyze by sections
    data.show()
    QtGui.QApplication.instance().exec_()

