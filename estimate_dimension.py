__author__ = 'ruben pinzon'
__doc__ = 'Dimension of the intrinsic latent space using PCA and FA for the HC-5 database.' \
          'based on Alexandre Gramfort and Denis A. Engemann plot_pca_vs_fa_model_selection.py'

import numpy as np
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.covariance import ShrunkCovariance, LedoitWolf
from sklearn.cross_validation import cross_val_score
from sklearn.grid_search import GridSearchCV
import matplotlib.pyplot as plt
import scipy.ndimage.filters as fil

y = np.load('i01_maze05.005_firings.npy')

# ###############################Evaluate PCA and FA#################
max_dims = 40
n_components = np.arange(0, max_dims, 2)  # options for n_components


def firing_rate(X, bin_size=100):
    """Convolution of spike trains with Gaussian function"""
    cells, bins = X.shape
    Y = list()
    for i in range(cells):
        Y.append(fil.gaussian_filter1d(X[i, :], sigma=bin_size))
    return np.array(Y)


def compute_scores(X):
    pca = PCA()
    fa = FactorAnalysis()

    pca_scores, fa_scores = [], []
    for n in n_components:
        print 'Processing dimension {}'.format(n)
        pca.n_components = n
        fa.n_components = n
        pca_scores.append(np.mean(cross_val_score(pca, X)))
        fa_scores.append(np.mean(cross_val_score(fa, X)))

    return pca_scores, fa_scores


def shrunk_cov_score(X):
    shrinkages = np.logspace(-2, 0, 30)
    cv = GridSearchCV(ShrunkCovariance(), {'shrinkage': shrinkages})
    return np.mean(cross_val_score(cv.fit(X).best_estimator_, X))


def lw_score(X):
    return np.mean(cross_val_score(LedoitWolf(), X))


q, T = y.shape
X = y - np.mean(y, axis=1)[:, np.newaxis]

pca_scores, fa_scores = compute_scores(X.T)
n_components_pca = n_components[np.argmax(pca_scores)]
n_components_fa = n_components[np.argmax(fa_scores)]
plt.plot(X.T)
pca = PCA(n_components='mle')
pca.fit(X.T)
n_components_pca_mle = pca.n_components_

print("best n_components by PCA CV = %d" % n_components_pca)
print("best n_components by FactorAnalysis CV = %d" % n_components_fa)
print("best n_components by PCA MLE = %d" % n_components_pca_mle)

plt.figure()
plt.plot(n_components, pca_scores, 'b', label='PCA scores')
plt.plot(n_components, fa_scores, 'r', label='FA scores')
plt.axvline(n_components_pca, color='b',
            label='PCA CV: %d' % n_components_pca, linestyle='--')
plt.axvline(n_components_fa, color='r',
            label='FactorAnalysis CV: %d' % n_components_fa, linestyle='--')
plt.axvline(n_components_pca_mle, color='k',
            label='PCA MLE: %d' % n_components_pca_mle, linestyle='--')

# compare with other covariance estimators
plt.axhline(shrunk_cov_score(X), color='violet',
            label='Shrunk Covariance MLE', linestyle='-.')
plt.axhline(lw_score(X), color='orange',
            label='LedoitWolf MLE' % n_components_pca_mle, linestyle='-.')

plt.xlabel('nb of components')
plt.ylabel('CV scores')
plt.legend(loc='lower right')

plt.show()
