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
        self.numN = max(self.clusters)
        self.spikes = self.get_var('res')
        self.fs = float(self.get_var('SamplingFrequency'))
        self.celltype = self.get_var('isIntern') #TODO: separate spks based on neuron type
        self.spk_per_neuron = self.get_neurons()
        self.ave_firing = self.firing_rate()
        self.laps = np.trim_zeros(self.get_var('StartLaps'))
        self.eeg = self.get_var('eeg')
        self.time = np.linspace(start=0., stop=len(self.eeg) / self.fs, num=len(self.eeg))
        self.sectionIO = self.get_var('MazeSectEnterLeft')
        self.tmax = self.get_var('SyncOff')
        self.panels = {'stats': 2, 'run': 1, 'whl': 0}
        self.panels_count = {'run': 0, 'whl': 0, 'stats': 0}
        self.autocorr = None
        self.autocorr_mean = None  # mean of the autocorrelogram
        self.fr_smo = None  # smooth firing rates

        self.x = self.get_var('X')
        self.y = self.get_var('Y')
        self.speed = self.get_var('speed')
        self.whspeed = self.get_var('WhlSpeedCW')
        self.section = self.get_var('MazeSection')
        self.theta = self.get_var('thetaPh')
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
        for neuron in range(self.numN):
            if 'neuron {}'.format(neuron) in self.spk_per_neuron:
                spikes = self.spk_per_neuron['neuron {}'.format(neuron)]
                if lap == 0:
                    spk = spikes
                    widget['run'].plot(spk / self.fs, np.ones(np.shape(spk)) + space, symbol='s',
                                          symbolPen=(space, self.numN), pen=None, symbolSize=1)
                    # TODO: has to add a new case for the raster of the whole record
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
            pw.showGrid(x=True, y=True)

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
                y = sfil.gaussian_filter1d(psth * self.fs / (bin_size * self.numN), sigma=self.fs * .1)
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
        intervals = {'run': (enter[0, 0], enter[4, 1]-1 if enter[4, 1] != 0 else enter[5, 1]-1),
                     'whl': (enter[11, 0], enter[12, 1]-1)}
        return intervals

    def firing_rate(self):
        firing = list()
        for neuron in range(self.numN):
            spk = self.spk_per_neuron['neuron {}'.format(neuron)]
            firing.append(float(len(spk)) * self.fs / (spk[-1] - spk[0]) if len(spk) != 0. else 0.)
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
        full = H / W

        return full[0:np.ceil(width / bin)]

    def autocorrelogram(self, width=0.002, bin=0.001, plot=True):
        """
        Returns an autocorrelogram with lag in [-width,width] and given bin size.
        T is the total duration (optional) and should be greater than the duration of T1 and T2.
        The result is in Hz (rate of coincidences in each bin).
        """
        aucorr = list()
        mean_aucorr = list()
        for neuron in range(self.numN):
            spk = self.spk_per_neuron['neuron {}'.format(neuron)] / self.fs
            aucorr.append(self.correlogram(T1=spk, T2=spk, width=width, bin=bin) if len(spk) != 0 else np.array([0, 0]))
            idx = np.abs(aucorr[neuron] - np.mean(aucorr[neuron])).argmin()
            mean_aucorr.append(idx * bin)
        if plot:
            pw = pg.PlotWidget()
            cnt = 0
            for value in aucorr:
                pw.plot(value, pen=(cnt, self.numN))
                cnt += 1
            self.layout.addWidget(pw, self.panels_count['stats'], 2)
            self.panels_count['stats'] += 1
        self.autocorr = aucorr
        self.autocorr_mean = mean_aucorr

    def plotxy(self, varx, vary):

        pw = pg.PlotWidget()
        y = self.__getattribute__(vary)
        x = self.__getattribute__(varx)
        pw.plot(x, y, pen=None, symbol='o', symbolpen='b')

        self.layout.addWidget(pw, self.panels_count['stats'], 2)
        self.panels_count['stats'] += 1

    def smooth_spikes(self, lap, plot=True):

        if lap != 0:
            widget = {'run': pg.PlotWidget(), 'whl': pg.PlotWidget()}
            space = 0
            win = self.extract_setions(lap)
            #   firing rates smoothed per panel
            fr_smo = {'run': list(), 'whl': list()}
            for neuron in range(self.numN):
                if 'neuron {}'.format(neuron) in self.spk_per_neuron:
                    spk = self.spk_per_neuron['neuron {}'.format(neuron)]
                    for panel, w in win.iteritems():
                        idx = np.where(np.logical_and(spk > w[0], spk <= w[1]))
                        y = np.zeros(w[1] - w[0] + 1)
                        y[spk[idx] - w[0]] = 1.
                        fr = sfil.gaussian_filter1d(y, sigma=self.fs * .1)
                        fr_smo[panel].append(fr)
                        if plot:
                            time = np.linspace(0, (w[1] - w[0]) / self.fs, num=len(fr))
                            if np.mean(fr) > (.5 / self.fs):
                                widget[panel].plot(time, fr/max(fr) + space, pen=(space, self.numN))

                    space += 1

            if plot:
                for key, pw in widget.iteritems():
                    self.layout.addWidget(pw, self.panels_count[key], self.panels[key])
                    self.panels_count[key] += 1
                    if self.rasterWidget:
                        pw.setXLink(self.rasterWidget[key])
            self.fr_smo = fr_smo

    def sort_firing(self, plot=True):

        fr_sort = {'run': list(), 'whl': list()}
        widget = {'run': pg.PlotWidget(), 'whl': pg.PlotWidget()}

        for panel in self.panels:
            idx = list()
            if panel in self.fr_smo:
                firing = self.fr_smo[panel]
                for f in firing:
                    idx.append(np.argmax(f))
                fr_sort[panel] = np.argsort(idx)
                if plot:
                    space = 0
                    for n in fr_sort[panel]:
                        fr = self.fr_smo[panel][n]
                        if np.mean(fr) > (.5 / self.fs):
                            time = np.linspace(0., len(fr)/self.fs, num=len(fr))
                            widget[panel].plot(time, fr/max(fr) + space, pen=(space, self.numN))
                        space += 1
                    for key, pw in widget.iteritems():
                        self.layout.addWidget(pw, self.panels_count[key], self.panels[key])
                        self.panels_count[key] += 1
                        if self.rasterWidget:
                            pw.setXLink(self.rasterWidget[key])

    def plot_maze(self, neuron):

        pw = pg.PlotWidget()
        y = self.x
        x = self.y
        pw.plot(x, y, pen=None, symbol='o', symbolpen='b')

        self.layout.addWidget(pw, self.panels_count['stats'], 2)
        self.panels_count['stats'] += 1

    def no_inh(self):
        """remove inhibitory neurons from spikes"""
        for i in range(len(self.celltype)):
            if self.celltype[i] == 1:
                print 'Neuron {} is inhibitory, removing'.format(i)
                del self.spk_per_neuron['neuron {}'.format(i)]
        self.numN = len(self.spk_per_neuron)