__author__ = 'ruben'
__doc__ = 'MCMC GPFA step-by-step based on PyMC'

import neurolib as np

folder = '/media/bigdata/'
names = np.find_files(folder)

animal = 0
data = np.get_cells(names[animal][1], only_pyr=True, section='Run')
# sanity check: raster one lap
np.raster([x[0] for x in data], title='{} Lap 0'.format(names[animal][0]))