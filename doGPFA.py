__author__ = 'ruben'
__doc__ = 'MCMC GPFA step-by-step based on PyMC'

import neurolib as nl
import numpy as np
from sklearn.decomposition import FactorAnalysis as FA
import matplotlib.pyplot as plt

folder = '/media/bigdata/'
names = nl.find_files(folder)

# Animal to process, length to produce uniform spike trains, and minimum firing rate;
# Neurons below that firing rate are removed. Bin_size is in ms. Latent_dim must be
# smaller than the num of cells.
animal = 1
window = 3.
fs = 1250
min_firing = 0.2
bin_size = 0.05
latent_dim = 10
cells, space = nl.get_cells(names[animal][1], only_pyr=True, section='Run')
# sanity check: raster one lap

# sanity check: show the position of the animal per lap
nl.plot_position(space, title='{} Lap 0'.format(names[animal][0]))
nl.raster([x[0] for x in cells], title='{} Lap 0'.format(names[animal][0]))

y = nl.spike_train(cells, length=int(window * fs), threshold=min_firing)
num_cells, _, num_laps = np.shape(y)
# bin the spike trains and do square root transform, remove mean
# X is the data to perform the FA per lap
X = np.sqrt(nl.binned(y, bin_size))
d = np.mean(X, axis=1)
assert latent_dim < num_cells, 'Latent dimension must be smaller than the number of cells'

# perform FA on centered data, per lap, to get the latent variables that will be regressed with GP
_, samples, _ = np.shape(X)
X_centered = X - d[:, np.newaxis]
C = np.zeros([latent_dim, samples, num_laps])
for lap in range(num_laps):
    fa = FA()
    fa.n_components = latent_dim
    C[:, :, lap] = fa.fit_transform(X_centered[:, :, lap].T).T

fig = plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
plt.plot(C[4, :, 0], 'x')

# regression over the latent variables using GP
