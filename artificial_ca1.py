__author__ = 'ruben'
__doc__ = 'Based on: Interpreting Neuronal Population Activity by Reconstruction: Unified Framework With' \
          ' Application to Hippocampal Place Cells Kechen Zhang , Iris Ginzburg , Bruce L. McNaughton' \
          ' , Terrence J. Sejnowski Journal of NeurophysiologyPublished 1 February 1998Vol. 79no. ' \
          '1017-. and code from https://github.com/nwilming/hpdecode.git'

from numpy import *
from scipy.stats import norm, poisson, gamma
import matplotlib.pyplot as plt
import datetime as dt
import scipy.integrate as integral
import csv
import os

def place_field(xmax=100, firing_rate=0.1, baseline=0.0001, **kwargs):
    """
    Creates a 1D Gaussian place field with center pos and
    covariance matrix. The max is scaled to desired firing_rate.
    Baseline gives the baseline firing rate.

    :return pdf: Probability density function
    """
    if 'preloaded' in kwargs:
        pos = kwargs['preloaded'][0]
        var = kwargs['preloaded'][1]
        n_modes = len(pos)
    else:
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

    def info():
        parameters = (pos, var)
        return parameters

    return pdf, info


def setup(basepath, n_p_fields=100, xmax=100):
    """
    Setup a linear track with dimmensions 1 x 100 mm
    :param n_p_fields: number of place fields
    """
    if not os.path.isfile(basepath + 'place_fields.csv'):
        p_fields = list()
        for f in range(n_p_fields):
            p_fields.append(place_field(xmax))
    else:
        print 'Place fields CVS file found {}'.format(basepath + 'place_fields.csv')
        p_fields = load_fields(basepath)

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


def save_rates(rates, name):
    savetxt(name, rates)


def sample(pos, time, p_fields, fig=None, **kwargs):
    ax2 = None
    if fig:
        ax2 = update_figure(fig, **kwargs)

    all_spikes = list()
    all_rates = list()
    for p, t in zip(pos, time):
        rates, obs_spikes = simulate_spikes(p_fields, p)
        for c, spks in enumerate(obs_spikes):
            if spks != 0:
                times = linspace(t, t + 0.001, spks)
                cell = ones(spks) * c
                if fig:
                    ax2.plot(times, ones(spks) * c, '|', color='r')
                all_spikes.extend(column_stack((times, cell)))
        all_rates.append(rates)
    # TODO: firing rate should not be over fmax, and cell latency should be added
    return ax2, all_spikes, all_rates


def kern(x, xprime, variance=1.0, lengthscale=1.0):
    return exp(-variance * abs(x - xprime))


def save_fields(p_fields, basepath):
    if not os.path.isfile(basepath + 'place_fields.csv'):
        pars = [p[1]() for p in p_fields]
        with open(basepath + 'place_fields.csv', 'w') as f:
            writer = csv.writer(f)
            f.write("mean - var\n")
            for line in pars:
                writer.writerows(line)
    else:
        print "File not overwritten"


def load_fields(basepath):
    with open(basepath + 'place_fields.csv', 'r') as f:
        reader = csv.reader(f)
        # get rid of header
        reader.next()
        p_fields = list()
        pos_flag = 0
        for row in reader:
            pos_flag += 1
            if pos_flag == 1:
                pos = double(row)
            elif pos_flag == 2:
                var = double(row)
                p_fields.append(place_field(preloaded=(pos, var)))
                pos_flag = 0

    return p_fields


if __name__ == '__main__':

    basepath = '/media/bigdata/synthetic/db11/'
    extra_hd = 'spw'
    show = True
    n_cells = 100  # number of neurons
    n_laps = 50  # number of laps
    speed = 10.  # cm/s
    min_pos = 75
    max_pos = 90  # linear track is 100 mm long
    t_taken = (max_pos - min_pos) / (speed * 10)  # time taken to complete the task
    Fs = 1000  # sampling frequency
    time = linspace(0, t_taken, int(t_taken * Fs))

    # create the place fields of the neurons or load one
    cells = setup(basepath, n_cells)

    p_fields = [c[0] for c in cells]
    save_fields(cells, basepath)

    # velocity of the rat in the linear track as a Ornstein-Uhlenbeck process
    mean_vel = speed * ones(time.size)
    alpha = 5.0  # OU process is generated as a Gaussian Process with exponential kernel with parameter alpha
    K = zeros((time.size, time.size))
    for i in xrange(time.size):
        for j in xrange(time.size):
            K[i, j] = kern(time[i], time[j], variance=alpha)

    vel_rat = random.multivariate_normal(mean_vel, K)
    # position is computed as a time integral of the velocity
    pos_rat = 10 * integral.cumtrapz(vel_rat, time, initial=0.) + min_pos

    # save the position and velocity of the animal in a txt file
    key = '_{}.'.format(dt.date.today())
    name = 'time_position_velocity_'
    savetxt(basepath + name + key + 'txt', vstack((time, pos_rat, vel_rat)).T, fmt='%3.5f')

    for n in range(n_laps):
        print 'Lap {} out of {} completed'.format(n, n_laps)
        _, spikes, rates = sample(pos_rat, time, p_fields, fig=None)
        key = 'spikes_{}_{}_lap_{}.'.format(n_cells, dt.date.today(), n)
        save_spks(spikes, basepath + extra_hd + key + 'txt')
        # save_rates(rates, basepath + key.replace('spikes', 'rates') + 'txt')

    # Show place fields distribution, position and velocity, and spikes
    if show:
        fig = plt.figure(frameon=False, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        plt.subplots_adjust(hspace=0.5)

        for p in p_fields:
            ax.plot(p(range(0, 100)))
        plt.ylabel('Place fields')
        plt.xlabel('Linear track (mm)')

        ax1 = update_figure(fig)
        ax1.plot(time, pos_rat)
        ax1.plot(time, 10 * vel_rat)
        plt.ylabel('Pos./vel. rat (mm, mm/s)')
        plt.xlabel('Time (s)')
        ax2, spikes, rates = sample(pos_rat, time, p_fields, fig, link=1)
        plt.ylabel('Cell Num.')
        plt.xlabel('Time')
        fig.savefig(basepath + name + extra_hd + key + 'png')
        plt.show()
