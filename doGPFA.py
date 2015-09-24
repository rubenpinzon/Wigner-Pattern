__author__ = 'ruben'
__doc__ = 'GPFA implementation based on Aecker class and Byron Yu'

import scipy.io as sio
import numpy as np
import neurolib as nl
from matplotlib import pyplot as plt

folder = '/media/bigdata/'
names = nl.find_files(folder)

# Animal to process, length to produce uniform spike trains, and minimum firing rate;
animal = 3
window = 3.
fs = 1250
bin_size = 0.05
# get_cell extract the spike trains for the section specified (only Run currently implemented),
# and filters to only pyramidal cells
cells, trajectory, duration = nl.get_cells(names[animal][1], only_pyr=True, section='Run')
# convert the spike times to spike trains with duration = window and filter neurons with
# firing rate below the threshold
y = nl.spike_train(cells['spikes'], length=window * fs, threshold=0.2)
# bin the spike trains into segments of bin_size seconds, output is of dimensions N x Bins x Laps
y_binned = nl.binned(y, bin_size=bin_size)

# Following the notation of GPFA in Byron Yu et al. y:,t|x:,t = normal(C*x:,t + d, R)
# where y : partially observed spikes, x: hidden latent variables, d: mean vect, R, indepent. noise of cells
# y_pyr is the spike trains of only pyramidal cells. Those with average firing rate below thr will be filtered
# the spike trains will be binned according to bin_size
# just for development use only onle lap
y_obs = y_binned
gpfa = {}
q, t, laps = np.shape(y_obs)
gpfa['T'] = t
gpfa['N'] = laps
gpfa['q'] = q
gpfa['d'] = np.reshape(y_obs, [q, -1]).mean(axis=1)
# y is of dimension N x (T * laps)
gpfa['y'] = np.reshape(y_obs, [q, -1]) - gpfa['d'][:, np.newaxis]


# =================Data to test the GPFA implementation========================


import numpy as np
import neurolib as nl

from scipy.linalg import toeplitz

N = 1  # laps
T = 50
p = 2
q = 16
gamma = np.log(1 / np.array([4., 1.]) ** 2)
sigmaN = 1.e-3

K = toeplitz(nl.covFun(np.arange(0, T), gamma[0], sigmaN))
x1 = np.linalg.cholesky(K).T.dot(np.random.rand(T, N))
K = toeplitz(nl.covFun(np.arange(0, T), gamma[1], sigmaN))
x2 = np.linalg.cholesky(K).T.dot(np.random.rand(T, N))
X = np.array([x1.flatten(), x2.flatten()])

phi = np.arange(0, q) / float(q) * 2. * np.pi
C = np.array([np.cos(phi), np.sin(phi)]).T / np.sqrt(q / 2)
D = np.random.rand(q, T)
S = np.identity(T)
Sn = np.kron(np.ones(N), S)
R = 0.02 * np.identity(q)
Y = np.linalg.cholesky(R).T.dot(np.random.normal(size=[q, T * N])) + C.dot(X) + D.dot(Sn)
Y_true = np.reshape(Y, [q, T, N])
# GPFA based on Matlab class: GPFA by Alexander S. Ecker
# Initialize observation model parameters using PCA. In the
# orinal implementation by Byron Yu the initialization is made
# via FastFa.
# p is the dimension of the latent unobservable space, C is the
# factor loadings and Lambda are the PCA eigen-valor
# SigmaN is the Innovation noise variance for the latent GP default 10-3

gpfa = {'T': T, 'N': N, 'q': q, 'p': p, 'd': np.reshape(Y_true, [q, -1]).mean(axis=1)}
gpfa['y'] = np.reshape(Y_true, [q, -1]) - gpfa['d'][:, np.newaxis]
gpfa['sigmaN'] = 1.e-3
gpfa['tol'] = 1.e-3
gpfa['loglike'] = list()
# initialize the FA factors with PCA with the p principal components
Q = np.cov(gpfa['y'])
Lambda, C = np.linalg.eig(Q)
ord = np.argsort(Lambda)[::-1][:gpfa['p']]
gpfa['C'] = C[:, ord]
Lambda = Lambda[ord]
# initialize private noise as residual variance not accounted
# for by PCA
gpfa['R'] = np.diag(Q - np.dot(gpfa['C'].dot(np.diagflat(Lambda)), gpfa['C'].T))
# initialize gammas of the covariance matrix of the GP
gpfa['gamma'] = np.log(0.01) * np.ones([gpfa['p'], 1])

iteration = 0
loglike_Base = np.nan
Y = gpfa['y'].reshape([q, T, N])
# EM algorithm
while iteration <= 2 or (gpfa['loglike'][-1] - gpfa['loglike'][-2]) / (gpfa['loglike'][-1] - loglike_Base) > gpfa[
    'tol']:

    iteration += 1
    # E step
    EX, VarX, loglike = nl.estimate_hidden(gpfa)
    gpfa['loglike'].append(loglike)

    # M step
    T1 = np.zeros([q, p + T])
    T2 = np.zeros((p + T) ** 2).reshape([p + T, -1])
    for t in range(T):
        x = EX[:, t, :]
        y = Y[:, t, :]
        T1[:, range(p)] += y.dot(x.T)
        T2[0:p, 0:p] += N * VarX[p * t:p * t + p, p * t:p * t + p] + x.dot(x.T)
        T1[:, p + t] = np.sum(y, axis=1)
        sx = np.sum(x, axis=1)
        T2[0:p, p + t] = sx
        T2[p + t, 0:p] = sx.T
        T2[p + t, p + t] = N

    CD = T1.dot(np.linalg.inv(T2))
    gpfa['C'] = CD[:, 0:p]
    gpfa['D'] = CD[:, p:p + T]
