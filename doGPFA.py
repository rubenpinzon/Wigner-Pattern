__author__ = 'ruben'

import scipy.io as sio
import numpy as np

files = {'rat005': '/media/bigdata/i01_maze05.005/i01_maze05_MS.005_BehavElectrData.mat',
         'rat006': '/media/bigdata/i01_maze06.002/i01_maze06_MS.002_BehavElectrData.mat'}

# Get the variables that we need from the MAT file, Spike, Cluster, isInhibitory,
data = {}
spikes = {}
spk_p_cell = {}
clusters = {}
isIntern = {}
sections = {}
for k, source in files.iteritems():
    data[k] = sio.loadmat(source)
    clusters[k] = np.squeeze(data[k]['Spike']['totclu'][0, 0])
    spikes[k] = np.squeeze(data[k]['Spike']['res'][0, 0])
    isIntern[k] = np.squeeze(data[k]['Clu']['isIntern'][0, 0]) == 1
    sections[k] = np.squeeze(data[k]['Par']['MazeSectEnterLeft'][0, 0])
    # Separate spikes by neuron number
    neuron = list()
    for n in range(0, max(clusters[k])):
        neuron.append(spikes[k][clusters[k] == n])
    spk_p_cell[k] = neuron
    spk_p_cell[k + '_numN'] = len(neuron)
print 'Completed'

# Get intervals of interest: Section 2 to 6 that correspond to the running sections.
# Section 1 seems to be between the running wheel and the central arm of the Maze.
#  It is not clear the boundary thus it is avoided
sect_run = {}
for k, val in sections.iteritems():
    sect_run[k + '_numLaps'] = len(val)
    win = list()
    for j in range(sect_run[k + '_numLaps']):
        win.append([val[j][1, 0], sum(val[j][6:8, 0])])
    sect_run[k + '_win'] = win
    sect_run[k + '_maxLen'] = np.diff(win).max()
    # TODO: add a max time for the running section.
    # TODO: Consider filtering the laps with wrong performance
    # Extract spikes and create an array of spike trains with dims N x length x Laps
    spk_train = np.zeros([spk_p_cell[k + '_numN'], sect_run[k + '_maxLen'], sect_run[k + '_numLaps']])
    for i in range(spk_p_cell[k + '_numN']):
        for j in range(sect_run[k + '_numLaps']):
            idx = np.where(np.logical_and(spk_p_cell[k][i] >= win[j][0], spk_p_cell[k][i] < win[j][1]))
            spk_train[i, spk_p_cell[k][i][idx] - win[j][0], j] = 1
    sect_run[k + '_spk_train'] = spk_train

# Following the notation of GPFA in Byron Yu et al. y:,t|x:,t = normal(C*x:,t + d, R)
# where y : partially observed spikes, x: hidden latent variables, d: mean vect, R, indepent. noise of cells
# y_pyr is the spike trains of only pyramidal cells. Those with average firing rate below thr will be filtered
animal = 'rat006'
minRate = 0.2
y = sect_run[animal + '_spk_train']
y_pyr = y[~isIntern[animal]]
m = np.mean(y_pyr.reshape([len(y_pyr), -1]), axis=1) * 1250.
keep_cells = m > minRate
y_pyr = y_pyr[keep_cells]