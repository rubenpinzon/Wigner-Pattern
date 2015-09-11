import scipy.io as sio
import scipy.ndimage.filters as sfil
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import os


# Old class to extract the data from the Matlab processed files of the hc-5 database
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
        self.celltype = self.get_var('isIntern')  # TODO: separate spks based on neuron type
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
        intervals = {'run': (enter[0, 0], enter[4, 1] - 1 if enter[4, 1] != 0 else enter[5, 1] - 1),
                     'whl': (enter[11, 0], enter[12, 1] - 1)}
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
                                widget[panel].plot(time, fr / max(fr) + space, pen=(space, self.numN))

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
                        if np.mean(fr) > (.2 / self.fs):
                            fr = fr - fr.min()
                            time = np.linspace(0., len(fr) / self.fs, num=len(fr))
                            widget[panel].plot(time, fr / max(fr) + space, pen=(space, self.numN))
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


def find_files(folder_base):
    """
    Finds the matlab files in the hc-5 database containing the experimental data.
    see hc-5 description for further details

    :param folder_base: the main folder containing the hC-5 database files
    :return: names: (name, paths) to the matlab files with the processed data
    """
    pattern_folder = 'i01_maze'
    patters_file = 'BehavElectrData'
    names = list()
    for case in os.listdir(folder_base):
        if case.__contains__(pattern_folder):
            for f in os.listdir(os.path.join(folder_base, case)):
                if f.__contains__(patters_file):
                    print os.path.join(folder_base, case, f)
                    names.append((case, os.path.join(folder_base, case, f)))
    return names


def get_cells(path, section=None, only_pyr=None, verbose=False):
    """
    Extract the spikes events from the MAT file of the HC-5 DB.
    if section is provided, then spikes are split accordingly

    :param path: Path to the processed MAT file
    :param section: Name of the section to extract (default: None)
                    Run, Wheel, or Other
    :param only_pyr: return only pyramidal cells
    :return neuron: (N, 2, Laps), where the second dimension of the list
                    contains the x-y coordinates of the N-spike
    :return trajectory_dict: dict with X, Y coordinates of animal in the running section
                    with 'left' and 'right' keys.

    """
    # TODO: implement section wheel and other sections if needed
    from scipy.ndimage import filters as fil

    data = sio.loadmat(path)
    clusters = np.squeeze(data['Spike']['totclu'][0, 0])
    spikes = np.squeeze(data['Spike']['res'][0, 0])
    isIntern = np.squeeze(data['Clu']['isIntern'][0, 0]) == 1
    sections = np.squeeze(data['Par']['MazeSectEnterLeft'][0, 0])
    X = np.squeeze(data['Track']['X'][0, 0])
    Y = np.squeeze(data['Track']['Y'][0, 0])
    x_spk = np.squeeze(data['Spike']['X'][0, 0])
    y_spk = np.squeeze(data['Spike']['Y'][0, 0])
    direction = np.squeeze(data['Laps']['TrialType'][0, 0])
    hit = np.squeeze(data['Par']['BehavType'][0, 0]) == 1
    # Separate spikes by neuron number
    neuron_xy = list()
    neuron_spk = list()
    trajectory_left = list()
    trajectory_right = list()
    trajectory_all = list()
    num_laps = len(sections)

    for n in range(1, max(clusters) + 1):

        if only_pyr and isIntern[n - 1]:
            continue
        spk = spikes[clusters == n]
        xy_cell = [x_spk[clusters == n], y_spk[clusters == n]]
        if verbose:
            print 'neuron {}-th with {} spks'.format(n, len(spk))

        if section == 'Run':
            # Get intervals of interest: Section 2 to 6 that correspond to the running sections.
            # Section 1 seems to be between the running wheel and the central arm of the Maze.
            #  It is not clear the boundary thus it is avoided
            if verbose:
                print 'neuron {}-th is pyramidal with {} spks'.format(n, len(spk))

            laps = list()
            xy_laps = list()
            for i in range(num_laps):
                start_run, end_run = sections[i][1, 0], sum(sections[i][4:6, 1])
                idx = np.where(np.logical_and(spk >= start_run, spk <= end_run))
                # save spike events aligned to the entering to sect 2.
                laps.append(spk[idx] - start_run)
                xy_laps.append((xy_cell[0][idx], xy_cell[1][idx]))
            neuron_spk.append(laps)
            neuron_xy.append(xy_laps)
        else:
            neuron_spk.append(spk)
            neuron_xy.append(xy_cell)

    # duration of each lap
    duration = list()
    arm = list()
    for i in range(num_laps):
        start_run, end_run = sections[i][1, 0], sum(sections[i][4:6, 1])
        duration.append(end_run - start_run)
        arm.append(direction[start_run])

    neuron = {'spikes': neuron_spk, 'xy': neuron_xy, 'direction': arm, 'hit': hit}
    # extract position per lap
    for i in range(num_laps):
        start_run, end_run = sections[i][1, 0], sum(sections[i][4:6, 1])
        # split the trajectories based on the Y position
        trajectory_all.append((X[start_run:end_run], Y[start_run:end_run]))
        if hit[i]:
            if Y[end_run] > Y[start_run]:
                trajectory_left.append((X[start_run:end_run], Y[start_run:end_run]))
                assert Y[end_run] > Y[start_run], 'lap {} is not to the left as supposed'.format(i)
            else:
                trajectory_right.append((X[start_run:end_run], Y[start_run:end_run]))
                assert Y[end_run] < Y[start_run], 'lap {} is not to the right as supposed'.format(i)
        else:
            print 'Skipping lap {} because animal failed it'.format(i)

    trajectory_dict = {'left': trajectory_left, 'right': trajectory_right}
    # compute the average trajectories by interpolating to the same length all
    #  the trajectories in the same direction

    trajectory_dict_interp = dict()
    for dir, trajectory in trajectory_dict.iteritems():
        _, max_time = np.shape(max(trajectory, key=lambda x: np.shape(x)[1]))
        trajectory_inter = list()
        for tra in trajectory:
            x, y = tra[0], tra[1]
            xp = np.interp(range(max_time), range(len(x)), x)
            yp = np.interp(range(max_time), range(len(y)), y)
            trajectory_inter.append([xp, yp])
        trajectory_dict_interp[dir + '_interp'] = trajectory_inter
        trajectory_median = np.median(trajectory_inter, axis=0)
        trajectory_dict_interp[dir + '_median'] = (
            fil.gaussian_filter1d(trajectory_median[0], 50.), fil.gaussian_filter1d(trajectory_median[1], 50.))

    trajectory_dict.update(trajectory_dict_interp)
    trajectory_dict.update({'all_traj': trajectory_all})

    print '{} cells extracted'.format(len(neuron['spikes']))
    print '{} Loading completed'.format(path)
    return neuron, trajectory_dict, duration


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


def plot_position(space, title='', **kwargs):
    """

    :param space: list of X,Y coordinates per lap

    """
    import matplotlib.pyplot as plt

    fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    for position in space:
        plt.plot(position[0], position[1], **kwargs)
    plt.xlabel('Position (mm)')
    plt.ylabel('Position (mm)')
    plt.title(title)
    plt.show()


def spike_train(spikes, length=1000, threshold=0.):
    """
    Conver spk events to matrix of spike trains (1's and 0's)
    :param spikes: spike events
    :param length: maximum length to extract
    :return: trains: matrix of spike trains (N x length) or (N x lenght x laps)

    """
    n, l = np.shape(spikes)
    trains = np.zeros([n, length, l])
    for icell, cell in enumerate(spikes):
        for ilap, lap in enumerate(cell):
            inside = lap[lap < length]
            trains[icell, inside, ilap] = 1.

    if threshold != 0.:
        m = np.mean(trains.reshape([n, -1]), axis=1) * 1250.
        keep = m >= threshold
        # print '{} Neurons removed with firing rate below {}'.format(sum(~keep), threshold)
        return trains[keep]
    return trains


def binned(train, bin_size=0.1):
    """
    Binned and square root transformed spike trains

    :param train:  spike trains
    :param bin_size: in ms
    :return: y : spike trains binned
    """
    # q x (Bins * Laps)
    q, t, laps = np.shape(train)
    bin_width = int(bin_size * 1250)
    T = int(t / bin_width)
    y = np.zeros([q, T, laps])
    for lap in range(laps):
        for iter_bin in range(T):
            bin_start, bin_end = iter_bin * bin_width, (iter_bin + 1) * bin_width - 1
            y[:, iter_bin, lap] = np.sum(train[:, bin_start:bin_end, lap], axis=1)
    return y


def rect_roi(centre, size, angle):
    ROI = ([size[0] / 2, size[1] / 2],
           [size[0] / 2, -size[1] / 2],
           [-size[0] / 2, -size[1] / 2],
           [-size[0] / 2, size[1] / 2],
           [size[0] / 2, size[1] / 2])
    rotation = ([np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)])
    return np.dot(ROI, rotation) + centre


def point_in_poly(x, y, poly):
    """ Determine if a point is inside a given polygon or not
        Polygon is a list of (x,y) pairs. This function
        returns True or False.  The algorithm is called
        the "Ray Casting Method.

    :param x: list
    :param y: list
    :param poly: list of (x,y) pairs
    :return: true false
    """
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def multi_point_in_poly(x, y, rois):
    count = list()
    for poly in rois:
        inside = 0
        for cx, cy in zip(x, y):
            if point_in_poly(cx, cy, poly.T):
                inside += 1
        count.append(inside)
    return count


def count_spikes(cells, arena_shape, grid_shape, verbose=-1):
    """counts the spike inside the given grid

    """
    if verbose != -1:
        import matplotlib.pyplot as plt
        print 'verbose activated, showing cell {}'.format(verbose)
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    num_bins = arena_shape / grid_shape
    rois = list()
    centers = list()
    for x in range(0, arena_shape[0], grid_shape[0]):
        for y in range(0, arena_shape[1], grid_shape[1]):
            roi = rect_roi([x, y], grid_shape, 0).T
            rois.append(roi)
            centers.append((x, y))
            if verbose != -1:
                plt.plot(roi[0], roi[1], color=[0.5, 0.9, 0.5])

    spike_count = list()
    cell_idx = 0
    for c in cells:
        cx = [x for xy in c for x in xy[0]]
        cy = [y for xy in c for y in xy[1]]
        if verbose != -1 and cell_idx == verbose:
            plt.plot(cx, cy, 'x', color='b')

        count = list()
        img = np.zeros(num_bins)
        for poly, center in zip(rois, centers):
            inside = 0
            for x, y in zip(cx, cy):
                if point_in_poly(x, y, poly.T):
                    inside += 1
            count.append(inside)
            idx = center / grid_shape
            img[idx[0], idx[1]] = inside
            if verbose != -1 and cell_idx == verbose:
                plt.text(center[0], center[1], '{}'.format(inside))
                plt.title('Cell {}'.format(cell_idx))

        if verbose != -1:
            print 'Spike counted for Cell {}'.format(cell_idx)
        spike_count.append(img)
        cell_idx += 1

    return spike_count


def covFun(time, tau, sigma):
    """ Gaussian process covariance function square exp"""
    sn = sigma
    sf = 1 - sn
    kernel = sf * np.exp(-0.5 * np.exp(tau) * time ** 2) + sn * (time == 0)
    return kernel


def estimate_hidden(model):
    """
    Estimate latent factors (and log-likelihood).
       EX, VarX, logLike = estX(Y) returns the expected
       value (EX) of the latent state X, its variance (VarX), and
       the log-likelihood (logLike) of the data Y.

    :param model:
    :return: Expentation hidden, variance and lok-likelihood
    """
    R = model['R']
    C = model['C']
    T = model['T']
    p = model['p']
    q = model['q']
    N = model['N']
    gamma_gp = model['gamma']
    sigma_noise = model['sigmaN']
    Y = model['y']

    from scipy.linalg import toeplitz

    # compute GP covariance and its inverse T*p x T*p
    Kb = np.zeros([T * p, T * p])
    Kbi = np.zeros([T * p, T * p])
    logdetKb = 0.

    for i in range(p):
        K = toeplitz(covFun(np.arange(0, T), gamma_gp[i], sigma_noise))
        # TODO: iplement a fast inverse method for symmetric Toeplitz matrices
        ndx = np.arange(i, T * p, p) + (T * p) * np.arange(i, T * p, p)[:, np.newaxis]
        Kbi.flat[ndx] = np.linalg.inv(K)
        Kb.flat[ndx] = K
        sig, det = np.linalg.slogdet(K)
        logdetKb += sig * det


    # Perform E tep (1) C'inv(R)*C
    RiC = np.divide(1., R)[:, np.newaxis] * C
    CRiC = np.dot(C.T, RiC)
    VarX = np.linalg.inv(Kbi + np.kron(np.identity(T), CRiC))
    logdetM = -np.product(np.linalg.slogdet(VarX))
    Cb = np.kron(np.identity(T), C)
    KbCb = Kb.dot(Cb.T)
    Rbi = np.kron(np.identity(T), np.diagflat(1. / R))
    RbiCb = np.kron(np.identity(T), RiC)
    CKCRi = Rbi - RbiCb.dot(VarX).dot(RbiCb.T)
    EX = KbCb.dot(CKCRi).dot(np.reshape(Y, [q * T, N]))
    # calculate log-likelihood
    val = -T * np.log(R).sum() - logdetKb - logdetM - q * T * np.log(2 * np.pi)
    norm_y = (np.divide(1., np.sqrt(R))[:, np.newaxis] * Y).flatten()
    CRiY = np.reshape(RiC.T.dot(Y), [p * T, N])
    loglike = 0.5 * (N * val - norm_y.T.dot(norm_y) + np.sum(CRiY * VarX.dot(CRiY)))
    print 'Done estimation with ll = {}'.format(loglike)

    return EX.reshape([p, T, N]), VarX, loglike


def construct_rois(bin_shape, path, verbose=False, color=[1, 0, 0]):
    """
    Creates the rois along an instructive path, since the left and right paths
    share the central portion of the maze

    :return: centers of the rois and rois
    """
    import math as mt
    from matplotlib import pyplot as plt
    origin = path[:, 0]
    dist = np.cumsum(np.sqrt(np.sum(np.diff(path - origin[:, np.newaxis]) ** 2, axis=0)))
    segments = int(np.floor(dist[-1] / bin_shape[0]))
    border_old = 0
    connect = True
    rois = list()
    centers = list()
    old_roi = ([0, 0])
    for i in range(segments):
        dist_bin = bin_shape[0] * (0.5 + i)
        border = np.where(np.diff(dist <= dist_bin))[0][0]
        delta = path[:, border_old] - path[:, border]
        center = (path[:, border_old] + path[:, border]) / 2
        angle = mt.atan2(-delta[1], delta[0])
        roi = rect_roi(center, bin_shape, angle)
        if connect and i > 0:
            roi[0] = old_roi[0]
            roi[4] = old_roi[0]
            roi[1] = old_roi[1]

        border_old = border
        old_roi = [roi[3], roi[2]]
        rois.append(roi.T)
        if verbose:
            plt.plot(roi.T[0], roi.T[1], color=color)
        centers.append(center)
    return centers, rois
