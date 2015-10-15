__author__ = 'ruben'
__doc__ = 'Based on: Interpreting Neuronal Population Activity by Reconstruction: Unified Framework With' \
          ' Application to Hippocampal Place Cells Kechen Zhang , Iris Ginzburg , Bruce L. McNaughton' \
          ' , Terrence J. Sejnowski Journal of NeurophysiologyPublished 1 February 1998Vol. 79no. ' \
          '1017-. and code from https://github.com/nwilming/hpdecode.git'

from numpy import *
from scipy.stats import norm, poisson, gamma
import matplotlib.pyplot as plt
import datetime as dt


def place_field(xmax=100, firing_rate=0.1, baseline=0.0001):
    """
    Creates a 1D Gaussian place field with center pos and
    covariance matrix. The max is scaled to desired firing_rate.
    Baseline gives the baseline firing rate.

    :return pdf: Probability density function
    """
    n_modes = floor(gamma.rvs(3, 0, 1))
    if n_modes < 1.:
        n_modes = 1
    if n_modes > 4:
        n_modes = 4

    pos = random.uniform(1, xmax, n_modes)
    var = random.uniform(1.5, xmax / 10, n_modes)

    gauss_m = list()
    for p, v in zip(pos, var):
        mv = norm(p, v)
        scale = mv.pdf(p)
        gauss_m.append((mv, scale))

    def pdf(arena):
        prob = 0.
        for g in gauss_m:
            prob += g[0].pdf(arena) / g[1]
        prob /= n_modes
        fr = firing_rate * prob + baseline
        return fr

    return pdf


def setup(n_p_fields=100, xmax=100):
    """
    Setup a linear track with dimmensions 1 x 100 mm
    :param n_p_fields: number of place fields
    """

    p_fields = list()
    for f in range(n_p_fields):
        p_fields.append(place_field(xmax))

    return p_fields


def simulate_spikes(tuning_curve, rx):
    """
    Compute firing rate for each neuron given place field center
    and sample number of observed spikes in one time unit.7

    :param tuning_curve: pdf of the tuning curve for each cell
    :param rx: position of the rat in the linear track
    :return rates, obs_spikes : firing rate as the rate of the Poisson distribution and spikes sampled from the pdf
    """
    rates = []
    obs_spikes = []
    for n, pfield in enumerate(tuning_curve):
        rate = pfield(rx)
        spikes = poisson.rvs(rate)
        rates.append(rate)
        obs_spikes.append(spikes)
    return rates, obs_spikes


def update_figure(fig, **kwargs):
    """updates the subplots in a figure and shares the x axis"""

    n = len(fig.axes)
    for i in range(n):
        fig.axes[i].change_geometry(n + 1, 1, i + 1)
    if 'link' in kwargs:
        link_to = kwargs['link']
        ax = fig.add_subplot(n + 1, 1, n + 1, sharex=fig.axes[link_to])
    else:
        ax = fig.add_subplot(n + 1, 1, n + 1)

    return ax


def save_spks(spikes, name):
    savetxt(name, spikes, '%3.5f\t%3d')


def sample(pos, time, p_fields, fig=None, **kwargs):
    ax2 = None
    if fig:
        ax2 = update_figure(fig, **kwargs)

    all_spikes = list()
    for p, t in zip(pos, time):
        rates, obs_spikes = simulate_spikes(p_fields, p)
        for c, spks in enumerate(obs_spikes):
            if spks != 0:
                times = linspace(t, t + 0.001, spks)
                cell = ones(spks) * c
                if fig:
                    ax2.plot(times, ones(spks) * c, '|', color='r')
                all_spikes.extend(column_stack((times, cell)))

    return ax2, all_spikes


if __name__ == '__main__':

    basepath = '/media/bigdata/synthetic/'
    show = True
    n_cells = 50
    p_fields = setup(n_cells)
    n_laps = 50

    speed = 5.  # cm/s
    max_pos = 80
    t_taken = max_pos / (speed * 10)
    Fs = 1000
    time = linspace(0, t_taken, int(t_taken * Fs))
    pos_rat = max_pos / (1 + exp(-2*time+1.))

    # Show place fields distribution
    if show:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        plt.subplots_adjust(hspace=0.5)

    for p in p_fields:
        ax.plot(p(range(0, max_pos)))
    plt.ylabel('Place fields')
    plt.xlabel('Linear track (mm)')

    ax1 = update_figure(fig)
    ax1.plot(time, pos_rat)
    plt.ylabel('Position animal (mm)')
    plt.xlabel('Time (s)')
    ax2, spikes = sample(pos_rat, time, p_fields, fig, link=1)
    plt.ylabel('Cell Num.')
    plt.xlabel('Time')
    key = 'spikes_{}_{}_lap_{}.'.format(n_cells, dt.date.today(), 0)
    fig.savefig(basepath + key + 'png')

    for n in range(n_laps):
        _, spikes = sample(pos_rat, time, p_fields, fig=None)
        key = 'spikes_{}_{}_lap_{}.'.format(n_cells, dt.date.today(), n)
        save_spks(spikes, basepath + key + 'txt')

    plt.show()
