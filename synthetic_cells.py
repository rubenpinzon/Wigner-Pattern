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


def setup(num_cells=100):
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
    for f in range(num_cells):
        center = [random.normal(xmax / 2, xmax / 2), random.normal(xmax / 2, xmax / 2)]
        covariance = random.normal(loc=50., scale=10., size=2)
        pfields.append(place_field(center, covariance))
        centers.append(center)
    delta = 0.1
    x = arange(0, xmax, delta)
    y = arange(0, ymax, delta)
    X, Y = meshgrid(x, y)

    return pfields, X, Y, centers


def simulate_spikes(p_fields, rx, ry):
    """
    Compute firing rate for each neuron given place field center
    and sample number of observed spikes in one time unit.
    """
    rates = []
    obs_spikes = []
    for n, pfield in enumerate(p_fields):
        rate = pfield((rx, ry))
        spikes = poisson.rvs(rate)
        rates.append(rate)
        obs_spikes.append(spikes)
    return rates, obs_spikes


def generate_position(speed=5.):
    """
    Simulates the position of the animal in the maze given a constant speed
    Alternates between left and right arms of the "T" maze

    :param: speed: Speed at which the animal explore the arena.
    :return: x, y position of the rate at given speed
    """


# def close_fields(centers, path):


if __name__ == '__main__':

    p_fields, ax, ay, centers = setup(500)
    # Show place fields distribution
    fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
    ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    pos = dstack((ax, ay))

    path_left = loadtxt('stereotypicalLEft.txt')[::40, :] / 10
    path_right = loadtxt('stereotypicalRight.txt')[::40, :] / 10

    plt.plot(path_left[:, 0], path_left[:, 1], linewidth=3.0, color='k')
    plt.plot(path_right[:, 0], path_right[:, 1], linewidth=3.0, color='k')
    # Define a radius from the path inside which place fields are accepted
    # so that any place_field' center within a distance <= r will be included.
    radius = 10.
    # compute distance to each center
    cell_idx = list()
    for idx, center in enumerate(centers):
        for p in range(len(path_right)):
            distL = sqrt(sum((center - path_left[p, :]) ** 2))
            distR = sqrt(sum((center - path_right[p, :]) ** 2))
            if distL <= radius or distR <= radius:
                cell_idx.append(idx)
                break

    cell_selected = [p_fields[i] for i in cell_idx]
    print '{} cells created'.format(shape(cell_selected)[0])
    for p in cell_selected:
        ax2.contour(ax, ay, p(pos))
    # Stereotypical path taken from a real animal

    plt.show()
