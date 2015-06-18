__author__ = 'ruben'

import argparse
import pywt
import scipy.io as sio
import numpy as np
import matplotlib.cm as cm
import pylab

parser = argparse.ArgumentParser(description='Function to calculate the wavelet packet transform')
parser.add_argument("-p", '--PATH', type=str, nargs='*',
                    help='Path to the mat file ##.BehavElectrData.mat')
parser.add_argument("-s", '--signal', type=str, nargs='*',
                    help='Name of the field to analyze (e.g., eeg)')
parser.add_argument("-w", "--wavelet", type=str, nargs=1,
                    help="name of the wavelet function. Default db4. Check families with get_wpt.py -f",
                    default='db4')
parser.add_argument("-l", "--level", type=int, nargs=1,
                    help="decomposition level (integer positive)",
                    default=4)
parser.add_argument("-f", "--wavenames", action='store_true',
                    help="names of wavelet function available")

color_ascii = {'green': '\033[1;32m{t}\033[1;m', 'red': '\033[1;31m{t}\033[1;m'}

if __name__ == '__main__':
    args = parser.parse_args()
    if args.wavenames:
        for family in pywt.families():
            print "%s family:" % family, ', '.join(pywt.wavelist(family))
    else:

        data = sio.loadmat(args.PATH[0])
        fields = args.signal
        signals = list()
        for key, value in data.iteritems():
            if type(value).__module__ == np.__name__:
                variables = value[0].dtype.names
                for f in fields:
                    if f in variables:
                        signals.append(value[0][f][0])
                        print 'Signal {s} found in field {f}'.format(s=f, f=color_ascii['green'].format(t=key))
        if signals.__sizeof__() == 0:
            print 'Signal not found. Terminating'
        else:
            wavelet = args.wavelet
            level = args.level[0]
            order = "freq"  # "normal"
            interpolation = 'nearest'
            cmap = cm.cool

            for s in signals:
                wp = pywt.WaveletPacket(data=np.squeeze(s), wavelet=wavelet, mode='zpd')
                nodes = wp.get_level(level, order=order)
                labels = [n.path for n in nodes]
                values = pylab.array([n.data for n in nodes], 'd')
                values = abs(values)
                f = pylab.figure()
                f.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
                pylab.subplot(2, 1, 1)
                pylab.title("Original signal")
                pylab.plot(s, 'b')

                ax = pylab.subplot(2, 1, 2)
                pylab.title("Wavelet packet coefficients at level %d" % level)
                pylab.imshow(values, interpolation=interpolation, cmap=cmap, aspect="auto",
                    origin="lower", extent=[0, 1, 0, len(values)])
                pylab.yticks(pylab.arange(0.5, len(labels) + 0.5), labels)

                pylab.show()