__author__ = 'ruben'
__doc__ = 'Based on: Interpreting Neuronal Population Activity by Reconstruction: Unified Framework With' \
          ' Application to Hippocampal Place Cells Kechen Zhang , Iris Ginzburg , Bruce L. McNaughton' \
          ' , Terrence J. Sejnowski Journal of NeurophysiologyPublished 1 February 1998Vol. 79no. ' \
          '1017-. and code from https://github.com/nwilming/hpdecode.git'

from numpy import *
from scipy.stats import multivariate_normal, poisson
import matplotlib.pyplot as plt
import datetime as dt


def place_field(pos, covariance, firing_rate=0.2, baseline=0.001):
    """
    Creates a 1D Gaussian place field with center pos and
    covariance matrix. The max is scaled to desired firing_rate.
    Baseline gives the baseline firing rate.

    :return pdf: Probability density function
    """
    mv = multivariate_normal(pos, covariance)
    scale_constant = mv.pdf(pos)

    def pdf(arena):
        fr = firing_rate * mv.pdf(arena) / scale_constant + baseline
        # try:
        #     fr[fr > firing_rate] = firing_rate
        # except TypeError:
        #     if fr > firing_rate:
        #         fr = firing_rate
        return fr

    return pdf


def setup(n_objects=100):
    """
    Setup a linear track with dimmensions 1 x 100 mm
    :param n_objects: number of place fields
    """
    xmax = 100
    pfields = list()
    centers = list()
    for f in range(n_objects):
        origin = random.normal(xmax / 2, xmax / 2)
        covariance = random.normal(loc=20., scale=8., size=1)
        pfields.append(place_field(origin, covariance))
        centers.append(origin)

    return pfields, centers


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
        p_fields, centers = setup(n_cells)

        speed = 5.  # cm/s
        max_pos = 80
        t_taken = max_pos / (speed * 10)
        Fs = 1000
        time = linspace(0, t_taken, int(t_taken * Fs))
        pos_rat = linspace(0, max_pos, int(t_taken * Fs))

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

        key = 'spikes_{}_{}.'.format(n_cells, dt.date.today())
        fig.savefig(basepath + key + 'png')
        save_spks(spikes, basepath + key + 'txt')
        plt.show()
