__author__ = 'ruben'

import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Returns the electrode number for each cluster in the db h-c5')
parser.add_argument('PATH', type=str, nargs='+',
                    help='Path to any .clu.N file, other will be located in the same folder')
parser.add_argument("--verbosity", help="increase output verbosity")
color_ascii = {'green': '\033[1;32m{t}\033[1;m', 'red': '\033[1;31m{t}\033[1;m'}

if __name__ == '__main__':
    args = parser.parse_args()
    folder = os.path.dirname(args.PATH[0])
    for files in os.listdir(folder):
        if 'clu' in files:
            elec_values = np.loadtxt(os.path.join(folder, files), delimiter='\n', dtype=int)
            np.delete(elec_values, 0)
            hist, bins = np.histogram(elec_values, bins=range(min(elec_values), max(elec_values)))
            print hist, bins



