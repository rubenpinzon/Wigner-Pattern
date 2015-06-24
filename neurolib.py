import scipy.io as sio
import scipy.ndimage.filters as sfil
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui


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
        self.data = sio.loadmat(self.path)
        self.vars = self.get_varnames()
        if self.verbose:
            print self.vars

        # Loading data
        self.clusters = self.get_var('totclu')
        self.spikes = self.get_var('res')
        self.fs = float(self.get_var('SamplingFrequency'))
        self.spk_per_neuron = self.get_neurons()
        self.ave_firing = self.firing_rate()
        self.laps = np.trim_zeros(self.get_var('StartLaps'))
        self.eeg = self.get_var('eeg')
        self.time = np.linspace(start=0., stop=len(self.eeg) / self.fs, num=len(self.eeg))
        self.numN = max(self.clusters)
        self.sectionIO = self.get_var('MazeSectEnterLeft')
        self.tmax = self.get_var('SyncOff')
        self.panels = {'stats': 2, 'run': 1, 'whl': 0}
        self.panels_count = {'run': 0, 'whl': 0, 'stats': 0}
        self.autocorr = None

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
                    hist, _ = np.histogram(a=times, bins=np.ceil((bin_max - bin_min) / bin_size),
                                           range=(bin_min, bin_max))
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
            pw.plot(y, pen=(1, self.numN))
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
        self.win.show()
        QtGui.QApplication.instance().exec_()

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

    def firing_rate(self):
        firing = {}
        for key, spikes in self.spk_per_neuron.iteritems():
            firing[key] = float(len(spikes)) * self.fs / (spikes[-1] - spikes[0]) if len(spikes) != 0. else 0.
        return firing

        # Second-order statistics from Brian

    def correlogram(self, T1, T2, width=0.020, bin=0.001, T=None):
        """
            Returns a cross-correlogram with lag in [-width,width] and given bin size.
            T is the total duration (optional) and should be greater than the duration of T1 and T2.
            The result is in Hz (rate of coincidences in each bin).
        """
        if (T1 == []) or (T2 == []):  # empty spike train
            return None
        # Remove units
        width = float(width)
        T1 = np.array(T1)
        T2 = np.array(T2)
        i = 0
        j = 0
        n = int(np.ceil(width / bin))  # Histogram length
        l = []
        for t in T1:
            while i < len(T2) and T2[i] < t - width:  # other possibility use searchsorted
                i += 1
            while j < len(T2) and T2[j] < t + width:
                j += 1
            l.extend(T2[i:j] - t)
        H, _ = np.histogram(l, bins=np.arange(2 * n + 1) * bin - n * bin)  # , new = True)

        # Divide by time to get rate
        if T is None:
            T = max(T1[-1], T2[-1]) - min(T1[0], T2[0])
        # Windowing function (triangle)
        W = np.zeros(2 * n)
        W[:n] = T - bin * np.arange(n - 1, -1, -1)
        W[n:] = T - bin * np.arange(n)

        return H / W

    def autocorrelogram(self, width=0.002, bin=0.001, plot=True):
        """
        Returns an autocorrelogram with lag in [-width,width] and given bin size.
        T is the total duration (optional) and should be greater than the duration of T1 and T2.
        The result is in Hz (rate of coincidences in each bin).
        """
        aucorr = {}
        for key, spikes in self.spk_per_neuron.iteritems():
            aucorr[key] = self.correlogram(T1=spikes/self.fs, T2=spikes/self.fs, width=width, bin=bin) if len(spikes) != 0 else np.array([0,0])
        if plot:
            pw = pg.PlotWidget()
            cnt = 0
            for key, value in aucorr.iteritems():
                pw.plot(value, pen=(cnt, self.numN))
                cnt += 1
            self.layout.addWidget(pw, 0, 2)
            self.panels_count['stats'] += 1

        return aucorr
