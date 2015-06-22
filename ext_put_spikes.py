__author__ = 'ruben'

import scipy.io as sio
import numpy as np
import argparse
import os

import numpy as np

parser = argparse.ArgumentParser(description='Shows the putative spike waveforms')
parser.add_argument('-p', '--PATH', type=str, nargs='+',
                    help='Path to the folder containing the .spk files')
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action='store_true')
color_ascii = {'green': '\033[1;32m{t}\033[1;m', 'red': '\033[1;31m{t}\033[1;m'}


if __name__ == '__main__':
    args = parser.parse_args()
    folder = os.path.dirname(args.PATH[0])
    for files in os.listdir(folder):
        if 'spk' in files:
            f = np.fromfile(os.path.join(folder, files), dtype=np.int16, )
            print files, len(f)

#f = h5py.File("mytestfile.hdf5", "w")

