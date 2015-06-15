#
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

import argparse

parser = argparse.ArgumentParser(description='Function to plot files from hc-5 database')
parser.add_argument('PATH', type=str, nargs='+',
                    help='Path to the file to be plotted.')

args = parser.parse_args()
data = sio.loadmat(args.PATH[0])

print 'Number of signals inside file:', len(data)
print "Number of samples {l}".format(l=len(data['sh2']), t=type(data))

fig = plt.figure(frameon=False, figsize=(17, 6), dpi=80, facecolor='w', edgecolor='k')
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')

for key, value in data.iteritems():
    try:
        plt.plot(value, linewidth=1, label=key)
        print "Field {} plotted".format(key)
    except ValueError:
        print "Field {} omitted".format(key)

# Place a legend above this legend, expanding itself to
# fully use the given bounding box.
plt.legend(loc=5, borderaxespad=0.)
plt.show()
