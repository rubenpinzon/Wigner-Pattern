__author__ = 'ruben'

# from pylab import *
#
# t = arange(0.01, 5.0, 0.01)
# s1 = sin(2*pi*t)
# s2 = exp(-t)
# s3 = sin(4*pi*t)
# ax1 = subplot(311)
# plot(t,s1)
# setp( ax1.get_xticklabels(), fontsize=6)
#
# ## share x only
# ax2 = subplot(312, sharex=ax1)
# plot(t, s2)
# # make these tick labels invisible
# setp( ax2.get_xticklabels(), visible=False)
#
# # share x and y
# ax3 = subplot(313,  sharex=ax1, sharey=ax1)
# plot(t, s3)
# xlim(0.01,5.0)
# show()


import matplotlib.cm as cm
import pylab

import pywt
#
# x = pylab.arange(0, 1, 1. / 512)
# data = pylab.sin((5 * 50 * pylab.pi * x ** 2))
#
# wavelet = 'db4'
# level = 4
# order = "freq"  # "normal"
# interpolation = 'nearest'
# cmap = cm.cool
#
# wp = pywt.WaveletPacket(data, wavelet, 'sym', maxlevel=level)
# nodes = wp.get_level(level, order=order)
# labels = [n.path for n in nodes]
# values = pylab.array([n.data for n in nodes], 'd')
# values = abs(values)
#
# f = pylab.figure()
# f.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
# pylab.subplot(2, 1, 1)
# pylab.title("linchirp signal")
# pylab.plot(x, data, 'b')
# pylab.xlim(0, x[-1])
#
# ax = pylab.subplot(2, 1, 2)
# pylab.title("Wavelet packet coefficients at level %d" % level)
# pylab.imshow(values, interpolation=interpolation, cmap=cmap, aspect="auto",
#     origin="lower", extent=[0, 1, 0, len(values)])
# pylab.yticks(pylab.arange(0.5, len(labels) + 0.5), labels)
# #pylab.setp(ax.get_xticklabels(), visible=False)
#
# #pylab.figure(2)
# #pylab.specgram(data, NFFT=64, noverlap=32, cmap=cmap)
# #pylab.imshow(values, origin='upper', extent=[-1,1,-1,1],
# # interpolation='nearest')
#
# pylab.show()

import numpy as np
import pylab as P

#
# The hist() function now has a lot more options
#

#
# first create a single histogram
#
mu, sigma = 200, 25
x = mu + sigma*P.randn(10000)

# the histogram of the data with histtype='step'
n, bins, patches = P.hist(x, 50, normed=1, histtype='stepfilled')
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# add a line showing the expected distribution
y = P.normpdf( bins, mu, sigma)
l = P.plot(bins, y, 'k--', linewidth=1.5)
P.show()