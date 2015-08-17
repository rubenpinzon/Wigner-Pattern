__author__ = 'ruben'
__doc__ = 'Creates synthetic place cells activity following a probabilistic approach with poisson distributions ' \
          'for the firing rates and Gaussian tuning curves. From https://github.com/nwilming/hpdecode.git' \
          'and Kechen Zhang et al. Interpreting Neuronal Population Activity by Reconstruction J of Neurophysiolog' \
          '1998'

from numpy import *
from scipy.stats import multivariate_normal, poisson
import matplotlib.pyplot as plt


def place_field(pos, covariance, firing_rate=10, baseline=1):
    """
    Creates a 2D Gaussian place field with center pos and
    covariance matrix. The max is scaled to desired firing_rate.
    Baseline gives the baseline firing rate.

    :return pdf: Probability density function
    """
    mv = multivariate_normal(pos, covariance)
    scale_constant = mv.pdf(pos)

    def pdf(arena):
        fr = firing_rate * mv.pdf(arena) / scale_constant + baseline
        try:
            fr[fr > firing_rate] = firing_rate
        except TypeError:
            if fr > firing_rate:
                fr = firing_rate
        return fr

    return pdf


def setup(num_objs=100):
    """
    Setup an arena with the same dimensions as Eva's experiment (100x120) and standard place fields
    modeled as Gaussians with random covariances. The arena is a T maze.

    :param samples  Number of samples from the X-Y length.
                    Determines the number of cells as X*Y/(samples**2)

    """
    # TODO: there is one section in the maze that is not included, the head of the "T". Should be included
    xmax = 100
    ymax = 100
    pfields = list()
    centers = list()
    for f in range(num_objs):
        origin = [random.normal(xmax / 2, xmax / 2), random.normal(xmax / 2, xmax / 2)]
        covariance = random.normal(loc=150., scale=10., size=2)
        pfields.append(place_field(origin, covariance))
        centers.append(origin)
    delta = 0.1
    x = arange(0, xmax, delta)
    y = arange(0, ymax, delta)
    X, Y = meshgrid(x, y)

    return pfields, X, Y, centers


def simulate_spikes(tuning_curve, rx, ry):
    """
    Compute firing rate for each neuron given place field center
    and sample number of observed spikes in one time unit.
    """
    rates = []
    obs_spikes = []
    for n, pfield in enumerate(tuning_curve):
        rate = pfield((rx, ry))
        spikes = poisson.rvs(rate)
        rates.append(rate)
        obs_spikes.append(spikes)
    return rates, obs_spikes


def extract_cells(path, radius=10.):
    """ Define a radius from the path inside which place fields are accepted
    # so that any place_field' center within a distance <= r will be included."""
    num_bins = len(path)

    # compute distance to each center
    cell_idx = list()
    for idx, center in enumerate(centers):
        for p in range(num_bins):
            distL = sqrt(sum((center - path_left[p, :]) ** 2))
            distR = sqrt(sum((center - path_right[p, :]) ** 2))
            if distL <= radius or distR <= radius:
                cell_idx.append(idx)
                break
    cells = [p_fields[i] for i in cell_idx]
    num_cells = shape(cells)[0]
    return num_cells, cells, num_bins


def sample(path, cells, verbose=False):
    spikes = list()
    for b, pos in enumerate(path):
        rates, obs_spikes = simulate_spikes(cells, pos[0], pos[1])
        spikes.append(obs_spikes)
        if verbose:
            print 'Number of spikes observed in bin {}-th, {}'.format(b, sum(obs_spikes))

    return array(spikes)


if __name__ == '__main__':
    show = True

    p_fields, ax, ay, centers = setup(500)
    # Show place fields distribution

    pos = dstack((ax, ay))

    path_left = loadtxt('stereotypicalLEft.txt')[::100, :] / 10
    path_right = loadtxt('stereotypicalRight.txt')[::100, :] / 10

    # Sample spikes from a given path
    num_cells, cells, num_bins = extract_cells(path_left)
    print '{} cells created'.format(num_cells)

    if show:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
        ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        for p in cells:
            ax2.contour(ax, ay, p(pos))
        plt.plot(path_left[:, 0], path_left[:, 1], linewidth=3.0, color='k')
        plt.plot(path_right[:, 0], path_right[:, 1], linewidth=3.0, color='k')

    # Stereotypical path taken from a real animal
    # TODO: below spaghetti must be coded
    spikesleft = sample(path_left, cells)
    spikesright = sample(path_right, cells)

    plt.show()
    savetxt('synthetic_spikesleft.txt', spikesleft)
    savetxt('synthetic_spikesright.txt', spikesright)

    spikesleft = sample(path_left, cells)
    spikesright = sample(path_right, cells)

    plt.show()
    savetxt('synthetic_spikesleft2.txt', spikesleft)
    savetxt('synthetic_spikesright2.txt', spikesright)
    # TODO: use the trajectories of the real animal to sample from the synthetic cells
