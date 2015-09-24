__author__ = 'ruben'
__doc__ = 'Python implementation of Gaussian Process FA/PCA (latent variable) based on gpfa by Byron Yu and ' \
          'generalized version by Alexander S. Ecker'
import numpy as np
from sklearn.utils import check_array


class GpLv(object):
    """
    Gaussian Process Latent Variable (PCA/FA)

    Based on Yu et al. 2009, J. Neurophys.
    The implementation is based on and heavily influenced by Byron Yu and John Cunningham's code
    and the generalized version by Alexander S. Ecker.
    """

    def __init__(self, n_latents=None, tol=1e-4, copy=True, max_iter=1000,
                 gp_noise=1e-3, gp_kernel='rbf', verbose=False):
        """
        Constructs a GpLv object with the following default parameters:

        :param gp_noise: smoothness of the functions priors (Default: 0.001)
        :param tolerance: Stopping criterion used for EM fitting (Default 0.0001)
        :param verbose: verbosed output (Default: false)
        :param cov_type: name of the covariance function (Default: rbf)
        :param kwargs: additional key:value options

        """
        self.n_latents = n_latents
        self.gp_noise = gp_noise
        self.tolerance = tol
        self.verbose = verbose
        self.copy = copy
        self.max_iter = max_iter

        if gp_kernel not in ['rbf', 'exp']:
            raise ValueError('Kernel {} for the Gaussian process not implemented'.format(gp_kernel))
        self.kernel = gp_kernel

    def fit(self, Y, x=None):
        """
        Fit the gplv model to observations Y using EM

        :param: Y: array-like, shape (n_cells, n_samples, n_trials)

        :return self
        """
        from scipy.optimize import minimize

        n_cells, n_samples, n_laps = Y.shape
        #  laps with same number of samples
        self.mean_ = np.mean(Y, axis=3)
        self.n_cells_ = n_cells
        self.n_laps_ = n_laps
        self.n_samples_ = n_samples
        Y0 = Y - self.mean_[:, np.newaxis]

        Y0 = np.reshape(Y0, [n_cells, -1])
        n_latents = self.n_latents

        if n_latents is None:
            n_latents = n_cells



        # initialize the latent variables with PCA with the n_latent components
        Q = np.cov(Y0)
        eigen_values, C = np.linalg.eig(Q)
        order = np.argsort(eigen_values)[::-1][:n_latents]
        self.mapping_ = C[:, order]
        eigen_values = eigen_values[order]
        # initialize private noise as residual variance not accounted for by PCA
        self.observation_noise_ = np.diag(Q - np.dot(C.dot(np.diagflat(eigen_values)), C.T))
        # initialize gammas of the covariance matrix of the GP
        self.gp_gamma_ = np.log(0.01) * np.ones([n_latents, 1])

        loglike = []
        loglike_Base = -np.inf
        tol = self.tolerance
        iteration = 0

        for i in xrange(self.max_iter):

            # E step
            EX, VarX, ll = self.inference(Y)
            loglike.append(ll)

            # M step
            T1 = np.zeros([n_cells, n_latents + n_samples])
            T2 = np.zeros((n_cells + n_samples) ** 2).reshape([n_cells + n_samples, -1])
            for t in range(n_samples):
                x = EX[:, t, :]
                y = Y[:, t, :]
                T1[:, range(n_latents)] += y.dot(x.T)
                idx = (n_latents * t) + np.arange(n_latents)
                T2[0:n_latents - 1, 0:n_latents - 1] += n_laps * VarX[idx, idx] + x.dot(x.T)
                T1[:, n_latents + t] = np.sum(y, axis=1)
                sx = np.sum(x, axis=1)
                T2[0:n_latents - 1, n_latents + t] = sx
                T2[n_latents + t, 0:n_latents - 1] = sx.T
                T2[n_latents + t, n_latents + t] = n_laps

            CD = T1.dot(np.linalg.inv(T2))
            self.mapping_ = CD[:, 0:n_latents - 1]
            self.mean_ = CD[:, n_latents:n_latents + n_samples]

            Y0 = Y - self.mean_[:, np.newaxis]
            Y0 = Y0.reshape([n_cells, -1])

            self.observation_noise_ = np.diag(np.mean(Y0 ** 2, axis=1) - np.sum(
                Y0.dot(EX.reshape([n_latents, -1]).T) * self.mapping_)/(n_samples*n_laps))

            #  optimize the GPs time scale
            self.gp_gamma_ = np.zeros(n_latents)

            for i in range(n_latents):
                ndx = range(i, n_samples * n_latents, n_latents)
                EXi = EX[i, :, :]
                EXX = VarX[ndx, ndx] + EXi.dot(EXi.T)/n_laps

    def Egamma(self, gamma, EXX):
        """
        Function to optimize the time scale of the GPs
        :param gamma:
        :param EXX:
        :return:
        """
        sigmaf = 1 - self.gp_noise
        t = np.arange(self.n_samples_)
        logdetKb, Kb, Kbi = rbf(t, gamma, self.gp_noise)


    def inference(self, Y):
        """
        Estimate latent factors (and log-likelihood).
           EX, VarX, logLike = estX(Y) returns the expected
           value (EX) of the latent state X, its variance (VarX), and
           the log-likelihood (logLike) of the data Y.

        :param model:
        :return: Expentation hidden, variance and lok-likelihood
        """
        R = self.observation_noise_, C = self.mapping_, T = self.n_samples_
        p = self.n_latents, q = self.n_cells_, N = self.n_laps_

        # compute GP covariance and its inverse T*p x T*p
        logdetKb, Kb, Kbi = self.covFuc()

        # Perform E step (1) C'inv(R)*C
        RiC = np.divide(1., R)[:, np.newaxis] * C
        CRiC = C.T.dot(RiC)
        VarX = np.linalg.inv(Kbi + np.kron(np.identity(T), CRiC))
        logdetM = -np.product(np.linalg.slogdet(VarX))
        Cb = np.kron(np.identity(T), C)
        KbCb = Kb.dot(Cb.T)
        Rbi = np.kron(np.identity(T), np.diagflat(1. / R))
        RbiCb = np.kron(np.identity(T), RiC)
        CKCRi = Rbi - RbiCb.dot(VarX).dot(RbiCb.T)

        Y0 = Y - self.mean_[:, np.newaxis]
        EX = KbCb.dot(CKCRi).dot(np.reshape(Y0, [q * T, N]))
        # calculate log-likelihood
        val = -T * np.log(R).sum() - logdetKb - logdetM - q * T * np.log(2 * np.pi)
        norm_y = (np.divide(1., np.sqrt(R))[:, np.newaxis] * Y0).flatten()
        CRiY0 = np.reshape(RiC.T.dot(Y0), [p * T, N])
        loglike = 0.5 * (N * val - norm_y.T.dot(norm_y) + np.sum(CRiY0 * VarX.dot(CRiY0)))

        if self.verbose:
            print 'Done inference with ll = {}'.format(loglike)

        return EX.reshape([p, T, N]), VarX, loglike

    def covFuc(self):
        """
        Calculate the Covarince of GP and its inverse
        :return: Kb : p covariance functions in a Block Matrix K
        :return: Kbi : inverse of Kb
        :return: logdetKb : log determinant of Kb

        """
        from scipy.linalg import toeplitz

        p = self.n_latents, cov_type = self.kernel, gp_gamma = self.gp_gamma_, gp_sigma = self.gp_noise
        T = self.n_samples_

        Kb = np.zeros(2 * T * p).reshape([T * p, -1])
        Kbi = np.zeros(2 * T * p).reshape([T * p, -1])
        logdetKb = 0.
        for i in range(p):
            # TODO : include a dynamic call for covariance functions
            variance = rbf(np.arange(0, T), gp_gamma[i], gp_sigma)
            K = toeplitz(variance)
            ndx = np.arange(i, T * p, p) + (T * p) * np.arange(i, T * p, p)[:, np.newaxis]
            Kbi.flat[ndx] = np.linalg.inv(K)
            Kb.flat[ndx] = K
            sig, det = np.linalg.slogdet(K)
            logdetKb += sig * det

        return logdetKb, Kb, Kbi

    @staticmethod
    def toy_example(num_cells=100, num_laps=10):
        """
        Create a toy example to test the gplv implementation including the number of cells indicated in the parameters
           two hidden variables and 100 samples per Lap
        :param num_cells:
        :param num_laps:

        """
        from scipy.linalg import toeplitz

        N = num_laps  # laps
        T = 100  # samples per lap
        p = 2  # hidden variables
        q = num_cells
        gamma = np.log(1 / np.array([4., 1.]) ** 2)
        gp_noise = 1.e-3

        K = toeplitz(rbf(np.arange(0, T), gamma[0], gp_noise))
        x1 = np.linalg.cholesky(K).T.dot(np.random.rand(T, N))
        K = toeplitz(rbf(np.arange(0, T), gamma[1], gp_noise))
        x2 = np.linalg.cholesky(K).T.dot(np.random.rand(T, N))
        X = np.array([x1.flatten(), x2.flatten()])

        phi = np.arange(0, q) / float(q) * 2. * np.pi
        C = np.array([np.cos(phi), np.sin(phi)]).T / np.sqrt(q / 2)
        d = np.random.rand(q, T)
        S = np.identity(T)
        Sn = np.kron(np.ones(N), S)
        R = 0.02 * np.identity(q)  # isotopic noise
        Y = np.linalg.cholesky(R).T.dot(np.random.normal(size=[q, T * N])) + C.dot(X) + d.dot(Sn)
        Y = np.reshape(Y, [q, T, N])

        return Y, X, S


# Covariance functions

def rbf(time, gamma, gp_noise):
    """ Exponential quadratic covariance function square exp"""

    sn = gp_noise
    sf = 1 - sn
    kernel = sf * np.exp(-0.5 * np.exp(gamma) * time ** 2) + sn * (time == 0)
    return kernel
