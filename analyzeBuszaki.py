__author__ = 'ruben'

import neurolib as nl


path = '/media/LENOVO/HAS/DB/hc-5/i01_maze06_MS.002_BehavElectrData.mat'
data = nl.LoadBuszaki(path)
data.path = path
lap = 3
data.raster(lap=lap)
data.plot('eeg', lap=lap, color='m')
#data.plot('whspeed', lap=lap)
#data.plot('speed', lap=lap)
data.autocorr = data.autocorrelogram(width=0.05, bin=1./data.fs)
#
# data.plot('eeg', lap=lap, color='r')
#data.get_psth(bin_size=0.05, lap=lap)
print data.firing_rate()
#   Analyze by sections
data.show()


