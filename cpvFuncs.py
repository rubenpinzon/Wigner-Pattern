__author__ = 'ruben'
__doc__ = 'File to generate observations from a Gaussian process with arbitrary covariance kernels'

import numpy as np
from matplotlib import pyplot as plt


def kern(x, xprime, variance=1.0, lengthscale=1.0):
    return np.exp(-variance * (xprime - x) / lengthscale)


num_cells = 3
num_hidden = 2
num_samples = 200
x_pred = np.linspace(0, 4, num_samples)[:, None]

alpha = 1.0
lengthscale = 1
K = np.zeros((x_pred.size, x_pred.size))
for i in xrange(x_pred.size):
    for j in xrange(x_pred.size):
        K[i, j] = kern(x_pred[i], x_pred[j], variance=alpha, lengthscale=lengthscale)

plt.imshow(K, interpolation='none')
plt.colorbar()

plt.figure(1)
X = np.zeros((2, x_pred.size))
for i in xrange(2):
    x = np.random.multivariate_normal(mean=np.zeros(x_pred.size), cov=K)
    plt.plot(x.flatten())
    X[i, :] = x

plt.figure(2)
plt.plot(X[0, :], X[1, :])

plt.figure(3)
phi = np.arange(num_cells) / (num_cells * 2. * np.pi)
C = np.array([-0.9, 0.8, -0.1, 0.7, -0.8, 0.1]).reshape([num_cells, 2])
D = np.random.rand(num_cells, num_samples)
R = np.linalg.cholesky(0.1 * np.identity(num_cells)).T  # isotopic noise
Y = R.dot(np.random.rand(num_cells, num_samples)) + C.dot(X)

plt.plot(Y.T)
