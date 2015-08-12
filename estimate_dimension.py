__author__ = 'ruben pinzon'
__doc__ = 'Dimension of the intrinsic latent space using PCA and FA for the HC-5 database.' \
          'based on Alexandre Gramfort and Denis A. Engemann plot_pca_vs_fa_model_selection.py'

import numpy as np
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.cross_validation import cross_val_score
import matplotlib.pyplot as plt

y_mat = np.loadtxt('binned_spks.txt')
# ###############################Evaluate PCA and FA#################
max_dims = 40
n_components = np.arange(0, max_dims, 2)  # options for n_components


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


q, T = y_mat.shape
X = y_mat

pca_scores, fa_scores = compute_scores(X.T)
n_components_pca = n_components[np.argmax(pca_scores)]
n_components_fa = n_components[np.argmax(fa_scores)]
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
plt.legend(loc='lower right')
plt.show()
