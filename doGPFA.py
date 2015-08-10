__author__ = 'ruben'

import scipy.io as sio
import numpy as np

files = {'rat005': '/media/bigdata/i01_maze05.005/i01_maze05_MS.005_BehavElectrData.mat',
         'rat006': '/media/bigdata/i01_maze06.002/i01_maze06_MS.002_BehavElectrData.mat'}

# Get the variables that we need from the MAT file, Spike, Cluster, isInhibitory,
data = {}
spikes = {}
spk_p_cell = {}
clusters = {}
isIntern = {}
sections = {}
for k, source in files.iteritems():
    data[k] = sio.loadmat(source)
    clusters[k] = np.squeeze(data[k]['Spike']['totclu'][0, 0])
    spikes[k] = np.squeeze(data[k]['Spike']['res'][0, 0])
    isIntern[k] = np.squeeze(data[k]['Clu']['isIntern'][0, 0]) == 1
    sections[k] = np.squeeze(data[k]['Par']['MazeSectEnterLeft'][0, 0])
    # Separate spikes by neuron number
    neuron = list()
    for n in range(0, max(clusters[k])):
        neuron.append(spikes[k][clusters[k] == n])
    spk_p_cell[k] = neuron
    spk_p_cell[k + '_numN'] = len(neuron)
print 'Completed'

# Get intervals of interest: Section 2 to 6 that correspond to the running sections.
# Section 1 seems to be between the running wheel and the central arm of the Maze.
#  It is not clear the boundary thus it is avoided
sect_run = {}
for k, val in sections.iteritems():
    sect_run[k + '_numLaps'] = len(val)
    win = list()
    for j in range(sect_run[k + '_numLaps']):
        win.append([val[j][1, 0], sum(val[j][6:8, 0])])
    sect_run[k + '_win'] = win
    sect_run[k + '_maxLen'] = np.diff(win).max()
    # TODO: add a max time for the running section.
    # TODO: Consider filtering the laps in which the animal failed the task
    # Extract spikes and create an array of spike trains with dims N x length x Laps
    spk_train = np.zeros([spk_p_cell[k + '_numN'], sect_run[k + '_maxLen'], sect_run[k + '_numLaps']])
    for i in range(spk_p_cell[k + '_numN']):
        for j in range(sect_run[k + '_numLaps']):
            idx = np.where(np.logical_and(spk_p_cell[k][i] >= win[j][0], spk_p_cell[k][i] < win[j][1]))
            spk_train[i, spk_p_cell[k][i][idx] - win[j][0], j] = 1
    sect_run[k + '_spk_train'] = spk_train

# Following the notation of GPFA in Byron Yu et al. y:,t|x:,t = normal(C*x:,t + d, R)
# where y : partially observed spikes, x: hidden latent variables, d: mean vect, R, indepent. noise of cells
# y_pyr is the spike trains of only pyramidal cells. Those with average firing rate below thr will be filtered
# the spike trains will be binned according to bin_size
M = {}
animal = 'rat006'
minRate = 0.2
bin_size = int(np.ceil(0.03 * 1250))
y = sect_run[animal + '_spk_train']
_, t, laps = np.shape(y)
T = t / bin_size
M['T'] = T
y_pyr = y[~isIntern[animal], :bin_size * T, :]
m = np.mean(y_pyr.reshape([len(y_pyr), -1]), axis=1) * 1250.
keep_cells = m > minRate
# number of cells q, Y is R(q x T)
M['q'] = sum(keep_cells)
y_pyr_fil = np.reshape(y_pyr[keep_cells], [M['q'], -1])

# q x (Bins * Laps)
y_obs = np.zeros([q, T * laps])
for n in range(M['q']):
    for i in range(T * laps):
        bin_start, bin_end = i * T, (i + 1) * T - 1
        y_obs[n, i] = sum(y_pyr_fil[n, bin_start:bin_end])
M['means'] = y_obs.mean(axis=1)
M['y'] = y_obs - M['means'][:, np.newaxis]
# q x Bins x Laps
# y_obs = y_obs.reshape([q, T, laps])

# GPFA based on Matlab class: GPFA by Alexander S. Ecker
# Initialize observation model parameters using PCA. In the
# orinal implementation by Byron Yu the initialization is made
# via FastFa.
# p is the dimension of the latent unobservable space, C is the
# factor loadings and Lambda are the PCA eigen-valor
# SigmaN is the Innovation noise variance for the latent GP default 10-3
M['sigmaN'] = 1.e-3
M['tol'] = 1.e-3
M['p'] = 10
Q = np.cov(M['y'])
Lambda, C = np.linalg.eig(Q)
ord = np.argsort(Lambda)[::-1][:M['p']]
M['C'] = C[:, ord]
Lambda = Lambda[ord]
# initialize private noise as residual variance not accounted
# for by PCA
M['R'] = np.diag(Q - np.dot(M['C'].dot(np.diagflat(Lambda)), M['C'].T))
# initialize gammas of the covariance matrix of the GP
M['gamma'] = np.log(0.01) * np.ones([M['p'], 1])


def covFun(time, tau, sigma):
    """ Gaussian process covariance function square exp"""
    sn = sigma
    sf = 1 - sn
    kernel = sf * np.exp(-0.5 * np.exp(tau) * time ** 2) + sn * (time == 0)
    return kernel


def extX(M):
    """  Estimate latent factors (and log-likelihood).
       EX, VarX, logLike = estX(Y) returns the expected
       value (EX) of the latent state X, its variance (VarX), and
       the log-likelihood (logLike) of the data Y.
    """
    R = M['R']
    C = M['C']
    T = M['T']
    p = M['p']
    gamma = M['gamma']
    sigmaN = M['sigmaN']
    Y = M['y']

    from scipy.linalg import toeplitz

    # compute GP covariance and its inverse T*p x T*p
    Kb = Kbi = np.zeros([T * p, T * p])
    logdetKb = 0.
    idx = np.arange(0, (T - 1) * p * T, p)

    for i in range(M['p']):
        K = toeplitz(covFun(np.arange(0, T - 1), gamma[i], sigmaN))
        # TODO: iplement a fast inverse method for symmetric Toeplitz matrices
        Kbi.flat[idx + 1] = np.linalg.inv(K)
        sig, det = np.linalg.slogdet(K)
        logdetKb += sig * det
    # Perform E step (1) C'inv(R)*C
    RiC = np.divide(1., R)[:, np.newaxis] * C
    CRiC = np.dot(C.T, RiC)
    VarX = np.linalg.inv(Kbi + np.kron(np.identity(T), CRiC))
    Cb = np.kron(np.identity(T), C)
    KbCb = Kb.dot(Cb.T)
    Rbi = np.kron(np.identity(T), np.diagflat(1./ R))
    RbiCb = np.kron(np.identity(T), RiC)
    CKCRi = Rbi - RbiCb.dot(VarX).dot(RbiCb.T)
    EX = KbCb.dot(CKCRi).dot(np.reshape(Y,[]))

def em(M):
    iter = 0
    ll_base = np.nan

    while iter <= 2 or (M['ll'][-1] - M['ll'][-2]) / (M['ll'][-1] - ll_base) > M['tol']:
        iter += 1
        # perform E step
        EX, VarX, ll = estX(M)
