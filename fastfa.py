# coding=utf-8
__author__ = 'ruben'

import numpy as np
from scipy.stats.mstats import gmean
import warnings


class FastFA(object):
    """ Factor analysis and probabilistic PCA

    Parameters
    ----------
        data :  data matrix (xDim x N)

                Each row of `data` represents a variable, and each column a single
                observation of all those variables
        based on ffa.m by Zoubin Ghahramani and fastfa.m by Byron Yu -- byronyu@stanford.edu
    """

    def __init__(self, data=0., z=None, typ=None, path=None):
        np.random.seed(0)
        if path is None:
            self.Y = data
        else:
            self.Y = np.loadtxt(path)
            print np.shape(self.Y)
        self.xDim = self.Y.shape[0]
        self.zDim = self.xDim if z is None else z
        self.N = self.Y.shape[1]
        if typ is None:
            self.typ = 'fa'
        else:
            self.typ = 'ppca'
        self.tol = 1.e-8
        self.itm = 1
        self.minVar = 0.01
        self.verbose = True
        self.covX = np.cov(m=self.Y)
        self.scale = self.cholesky_norm()
        if path is None:
            self.Lambda = np.random.normal(loc=0., scale=1., size=(self.xDim, self.zDim)) * np.sqrt(
                self.scale / self.zDim)
        else:
            self.Lambda = np.loadtxt(path.replace('data.txt', 'randdata.txt')) * np.sqrt(self.scale / self.zDim)

        self.Phi = np.diag(self.covX)
        self.d = np.mean(self.Y, axis=1)
        self.varFloor = self.minVar * np.diag(self.covX)
        self.trace_loglike = []

    def cholesky_norm(self, val=0., det=None):
        """
            In Bayesian data analysis, the log determinant of symmetric positive definite matrices often pops up as a
            normalizing constant in MAP. Oftentimes, the determinant of A will evaluate as infinite in Matlab although
            the log det is finite, so one can’t use log(det(A)). However using the properties of the logs and cholesky

        :return: Norm: scale based on determinant of the data
        """
        if det is None:
            if np.linalg.matrix_rank(self.Y) == self.xDim:
                factor = np.linalg.cholesky(self.covX)
                _, logdet = np.linalg.slogdet(factor)
                norm = np.exp(2. * logdet / self.xDim)
            else:
                warnings.warn('Data matrix is not full rank. Scaling with geometrical mean')
                rnk = np.linalg.matrix_rank(self.covX)
                eigvals = np.linalg.eigvals(a=rnk).sort()
                norm = gmean(eigvals[1:rnk])
        else:
            factor = np.linalg.cholesky(val)
            _, norm = np.linalg.slogdet(factor)

        return norm

    def emfa(self):
        """ expectation maximization steps for fa
        :return: log likelihood
        """
        i = np.identity(self.zDim)
        const = -self.xDim / 2. * np.log(2 * np.pi)
        loglike_i = 0.
        for k in range(self.itm):
            """ E step """
            i_phi = np.diag(1./self.Phi)
            i_phi_lambda = np.dot(i_phi, self.Lambda)
            mm = i_phi - np.dot(np.dot(i_phi_lambda, np.linalg.inv(i + np.dot(self.Lambda.T, i_phi_lambda))),
                                i_phi_lambda.T)
            beta = np.dot(self.Lambda.T, mm)
            covx_beta = np.dot(self.covX, beta.T)
            ezz = i - np.dot(beta, self.Lambda) + np.dot(beta, covx_beta)

            """ compute log likelihood"""
            loglike_old = loglike_i
            logdetm = self.cholesky_norm(val=mm, det=True)
            loglike_i = self.N * const + self.N * logdetm - .5 * self.N * np.sum(np.sum(mm * self.covX))

            if self.verbose:
                print 'EM iteration {} lik {} \n'.format(k, loglike_i)

            self.trace_loglike.append(loglike_i)

            """ M Step """
            self.Lambda = np.dot(covx_beta, np.linalg.inv(ezz))
            self.Phi = np.diag(self.covX) - np.sum(covx_beta * self.Lambda, axis=1)
            if self.typ is 'fa':
                self.Phi = np.maximum(self.varFloor, self.Phi)
            else:
                self.Phi = np.dot(np.mean(self.Phi), np.identity(self.xDim))

            if k <= 1:
                loglikebase = loglike_i
            elif loglike_i < loglike_old:
                print "The log likelihood did not increase ({} and {})".format(loglike_i, loglike_old)
            elif (loglike_i - loglikebase) < (1 + self.tol) * (loglike_old - loglikebase):
                print 'No further change in Log Likelihood value. Terminating'
                return
        return
