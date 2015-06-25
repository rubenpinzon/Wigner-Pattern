__author__ = 'ruben'

import neurolib as nl


path = '/media/LENOVO/HAS/DB/hc-5/i01_maze06_MS.002_BehavElectrData.mat'
data = nl.LoadBuszaki(path)
data.path = path
lap = 3
data.raster(lap=lap)
data.smooth_spikes(lap, plot=False)
data.sort_firing()
data.plot('eeg', lap=lap, color='m')
# data.plot('whspeed', lap=lap)
# data.plot('speed', lap=lap)
# data.autocorrelogram(width=0.05, bin=1./data.fs)
# data.plotxy('ave_firing', 'autocorr_mean')
#
# #
# data.get_psth(bin_size=0.05, lap=lap)
# for i in range(data.numN):
#     print 'Neuron {}, fr {}Hz, mean autocorrelogram {}ms'.format(i, data.ave_firing[i], data.autocorr_mean[i]*1000)
# #   Analyze by sections
data.show()
