import scipy.io as sio
import scipy.ndimage.filters as sfil
import numpy as np
import os


def find_files(folder_base):
    """
    Finds the matlab files in the hc-5 database containing the experimental data.
    see hc-5 description for further details

    :param folder_base: the main folder containing the hC-5 database files
    :return: names: (name, paths) to the matlab files with the processed data
    """
    pattern_folder = 'i01_maze'
    patters_file = 'BehavElectrData.mat'
    names = list()
    for case in os.listdir(folder_base):
        if case.__contains__(pattern_folder):
            for f in os.listdir(os.path.join(folder_base, case)):
                if f.__contains__(patters_file):
                    print os.path.join(folder_base, case, f)
                    names.append((case, os.path.join(folder_base, case, f)))
    return names


def get_cells(path, section=None, only_pyr=None, verbose=False, process=False):
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

    maze_regions = {'mid_arm': [0], 'pre_turn': [1], 'turn': [2, 3],
                    'lat_arm': [4, 5], 'reward': [6, 7], 'delay': [8, 11], 'wheel': [12]}
    alternations = ['left', 'right', 'errorLeft', 'errorRight']
    print 'Available sections of the maze {}'.format(maze_regions.keys())

    # data from the database
    data = sio.loadmat(path)
    clusters = np.squeeze(data['Spike']['totclu'][0, 0])
    spikes = np.squeeze(data['Spike']['res'][0, 0])
    isIntern = np.squeeze(data['Clu']['isIntern'][0, 0]) == 1
    sections = np.squeeze(data['Par']['MazeSectEnterLeft'][0, 0])
    X = np.squeeze(data['Track']['X'][0, 0])
    Y = np.squeeze(data['Track']['Y'][0, 0])
    x_spk = np.squeeze(data['Spike']['X'][0, 0])
    y_spk = np.squeeze(data['Spike']['Y'][0, 0])
    wheel = np.squeeze(data['Laps']['WhlSpeedCW'][0, 0])
    direction = np.squeeze(data['Laps']['TrialType'][0, 0])

    # Separate spikes by neuron number

    num_laps = len(sections)

    if section == 'Run':
        # for backwards compatibility
        maze_in = [0]
        maze_out = [5, 6]
    else:
        assert len(section) == 2, 'Sections should include start and end area: [turn, delay]'
        maze_in = maze_regions[section[0]]
        maze_out = maze_regions[section[1]]

    experiment = []
    for n_lap in range(num_laps):
        spk_data = []
        pos_data = []
        start_sect, end_sect = sum(sections[n_lap][maze_in, 0]), sum(sections[n_lap][maze_out, 1])
        t_type = alternations[direction[start_sect] - 1]
        animal_pos = (X[start_sect:end_sect], Y[start_sect:end_sect])
        wh_speed = wheel[start_sect:end_sect]

        for n_cell in range(1, max(clusters) + 1):
            if only_pyr and isIntern[n_cell - 1]:
                continue
            spk = spikes[clusters == n_cell]
            pos_xy = [x_spk[clusters == n_cell], y_spk[clusters == n_cell]]
            if verbose:
                print 'neuron {}-th with {} spks'.format(n_cell, len(spk))
            idx = np.where(np.logical_and(spk >= start_sect, spk <= end_sect))
            # spike times aligned to the start of the section
            spk_data.append(spk[idx] - start_sect)
            pos_data.append((pos_xy[0][idx], pos_xy[1][idx]))

        trial = Trials(t_id=n_lap, t_type=t_type, spikes=spk_data, spk_pos=pos_data,
                       animal_pos=animal_pos, section=section, wheel=wh_speed, process=process)
        experiment.append(trial)

    print '{} cells extracted'.format(experiment[0].n_cells)
    print '{} Loading completed'.format(path)
    return experiment


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


class Trials:
    n_trials = 0
    max_length = 0
    colors = {'left': [1, 0.8, 0.8], 'right': [0.8, 0.8, 1.0],
              'errorLeft': [0.2, 0.2, 0.2], 'errorRight': [0.2, 0.2, 0.2]}

    def __init__(self, t_id, t_type, spikes, spk_pos, animal_pos, section, wheel, process=False):
        self.id = t_id
        self.type = t_type
        self.spikes = spikes
        self.spikes_pos = spk_pos
        self.animal_pos = animal_pos
        self.section = section
        self.n_cells = len(spikes)
        self.color = self.colors[t_type]
        self.under_use = True
        self.wheel = wheel
        if process:
            self.spikes_lin = self.lin_spikes_pos()
        Trials.n_trials += 1

    def cumulative_dist(self):
        d = []
        X, Y = self.animal_pos
        x_org = X[0]
        y_org = Y[0]
        for x, y in zip(X, Y):
            d.append(np.sqrt((x - x_org) ** 2 + (y - y_org) ** 2))
            x_org = x
            y_org = y
        return np.cumsum(d)

    def disable(self):
        print 'Trial {}, type= {} disabled'.format(self.id, self.type)
        self.under_use = False

    def lin_spikes_pos(self):
        dn = np.linalg.norm(self.animal_pos, axis=0)
        [xa, ya] = self.animal_pos
        d = self.cumulative_dist()

        spikes_lin_pos = []
        for spk in self.spikes_pos:
            spkT = np.transpose(spk)
            index = []
            for xs, ys in spkT:
                idx = np.where(xa == xs) and np.where(ya == ys)
                index.extend(idx[0])
            spikes_lin_pos.append((index, d[index]))

        return spikes_lin_pos
