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
    for f in range(num_cells):
        center = [random.normal(xmax / 2, xmax / 2), random.normal(xmax / 2, xmax / 2)]
        covariance = random.normal(loc=50., scale=10., size=2)
        pfields.append(place_field(center, covariance))
    delta = 0.1
    x = arange(0, xmax, delta)
    y = arange(0, ymax, delta)
    X, Y = meshgrid(x, y)

    return pfields, X, Y


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


# if __name__ == '__main__':

p_fields, ax, ay = setup(100)
# Show place fields distribution
fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
pos = dstack((ax, ay))
for p in p_fields:
    ax2.contour(ax, ay, p(pos))
plt.show()
