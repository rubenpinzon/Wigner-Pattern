# from __future__ import print_function
#
# __author__ = 'ruben'
#
# """
# Show how to connect to keypress events
# """
# import sys
# import numpy as np
# import matplotlib.pyplot as plt
#
#
# def press(event):
#     print('press', event.key)
#     sys.stdout.flush()
#     if event.key=='x':
#         visible = xl.get_visible()
#         xl.set_visible(not visible)
#         fig.canvas.draw()
#
# fig, ax = plt.subplots()
#
# fig.canvas.mpl_connect('key_press_event', press)
#
# ax.plot(np.random.rand(12), np.random.rand(12), 'go')
# xl = ax.set_xlabel('easy come, easy go')
#
# plt.show()
# import numpy as np
# import matplotlib.pyplot as plt
#
# x = np.arange(-10,10)
# y = x**2
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(x,y)
#
# coords = []
#
# def onclick(event):
#     global ix, iy
#     ix, iy = event.xdata, event.ydata
#     print 'x = %d, y = %d'%(
#         ix, iy)
#
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# plt.show()
# fig.canvas.mpl_disconnect(cid)
# from matplotlib import pyplot as plt
#
# class LineBuilder:
#     def __init__(self, line):
#         self.line = line
#         self.xs = list(line.get_xdata())
#         self.ys = list(line.get_ydata())
#         self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
#
#     def __call__(self, event):
#         print 'click', event
#         if event.inaxes!=self.line.axes: return
#         self.xs.append(event.xdata)
#         self.ys.append(event.ydata)
#         self.line.set_data(self.xs, self.ys)
#         self.line.figure.canvas.draw()
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('click to build line segments')
# line, = ax.plot([0], [0])  # empty line
# linebuilder = LineBuilder(line)
#
# plt.show()

import numpy
from neuronpy.graphics import spikeplot
from neuronpy.math import kernel
import matplotlib.pyplot as plt


def get_psth(data, bin_size, **kwargs):
    if 'keys' in kwargs:
        keys = kwargs['keys']
    else:
        keys = data.keys()

    if len(keys) != 0:
        bin_min, bin_max = (0, numpy.ceil(max(data[max(data)])))
        psth = numpy.zeros(numpy.ceil(bin_max/bin_size), dtype=int)
        for k in keys:
            spike_bin = data[k]
            hist, _ = numpy.histogram(a=spike_bin, bins=numpy.ceil(bin_max/bin_size), range=(bin_min, bin_max))
            psth += hist
        return psth/(len(keys)*bin_size)
    else:
        return None


spikes = []
num_cells = 10
num_spikes_per_cell = 5
frequency = 20
dt = 5.0
# Make the spike data. Use a simple Poisson-like spike generator
# (just for illustrative purposes here. Better spike generators should
# be used in simulations).
spikes_dict = {}
for i in range(num_cells):
    isi = numpy.random.poisson(frequency, i + num_spikes_per_cell)
    spikes.append(numpy.cumsum(isi))
    spikes_dict['neuron {}'.format(i)] = spikes[i]

psth = get_psth(data=spikes_dict, bin_size=dt)
# spikes is now a list of lists where each cell has a list of spike
# times. Now, let's plot these spikes.

sp = spikeplot.SpikePlot(sth_ratio=0.3, savefig=False)
sth = sp.get_sth()
sth.set_dt(dt)
# k = kernel.gauss_1d(sigma=2., dt=dt)
# sth.set_kernel(k)
# sth.set_style('lineto')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(psth, marker='o')

sp.plot_spikes(spikes)
