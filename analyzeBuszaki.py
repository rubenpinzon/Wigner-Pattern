__author__ = 'ruben'
__doc__ = 'This script shows how to construct progressively a screen of analysis for the HC-5 database by calling the ' \
          'library neurolib.' \
          'USAGE: (Version 1.0)' \
          '1- Set the path to the math file animal_BehavElectrData.mat of interest in the HC-5 Database' \
          '2- Select the lap of interest' \
          '3- Add the analysis of interes as shown in the example. Each anaysis will have a plot in the screen'


import neurolib as nl

path = '/media/bigdata/i01_maze05.005/i01_maze05_MS.005_BehavElectrData.mat'
# Create a new object to read the database
data = nl.LoadBuszaki(path)
data.path = path
# Lap of interest
lap = 5
# If inhibitory neurons are to be filter out
data.no_inh()

# ======== Example of analysis of interest ===========
# windows will be stacked up one after the other
# rater plot for lap
data.raster(lap=lap)
# Compute smoother firing rates
data.smooth_spikes(lap, plot=True)
data.sort_firing()
# show eeg
data.plot('eeg', lap=lap, color='m')
# show wheel speed
data.plot('whspeed', lap=lap)
# show theta phase
data.plot('theta', lap=lap)
# show the marker for sections in the maze
data.plot('section', lap=lap)
# shows the animal's speed
data.plot('speed', lap=lap)
# data.autocorrelogram(width=0.05, bin=1./data.fs)
# show window
data.show()
