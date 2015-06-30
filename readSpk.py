__author__ = 'ruben'

import os
import tables
import time
import xml.etree.ElementTree as ET
import re
import numpy as np
import matplotlib.pyplot as plt


def read_xml(filename_xml, fileindex=1):
    """Read the XML file associated to the current dataset,
    and return a metadata dictionary."""

    tree = ET.parse(filename_xml)
    root = tree.getroot()

    d = {}

    ac = root.find('acquisitionSystem')
    if ac is not None:
        nc = ac.find('nChannels')
        if nc is not None:
            d['total_channels'] = int(nc.text)
        sr = ac.find('samplingRate')
        if sr is not None:
            d['rate'] = float(sr.text)

    sd = root.find('spikeDetection')
    if sd is not None:
        cg = sd.find('channelGroups')
        if cg is not None:
            # find the group corresponding to the fileindex
            g = cg.findall('group')[fileindex - 1]
            if g is not None:
                ns = g.find('nSamples')
                if ns is not None:
                    d['nsamples'] = int(ns.text)
                nf = g.find('nFeatures')
                if nf is not None:
                    d['fetdim'] = int(nf.text)
                c = g.find('channels')
                if c is not None:
                    d['nchannels'] = len(c.findall('channel'))

    if 'nchannels' not in d:
        d['nchannels'] = d['total_channels']

    if 'nsamples' not in d:
        ne = root.find('neuroscope')
        if ne is not None:
            sp = ne.find('spikes')
            if sp is not None:
                ns = sp.find('nSamples')
                if ns is not None:
                    d['nsamples'] = int(ns.text)

    # If no nFeatures, default to 3 (really old XML from Neuroscope).
    if 'fetdim' not in d:
        d['fetdim'] = 3

    # klusters tests
    metadata = dict(
        nchannels=d['nchannels'],
        nsamples=d['nsamples'],
        fetdim=d['fetdim'],
        freq=d['rate'])

    return metadata

class MemMappedBinary(object):
    def __init__(self, filename, dtype, rowsize=None):
        self.filename = filename
        self.dtype = dtype

        # Number of bytes of each item.
        self.itemsize = np.dtype(self.dtype).itemsize
        # Number of items in each row.
        self.rowsize = rowsize
        # Number of bytes in each row.
        self.rowsize_bytes = self.rowsize * self.itemsize
        # Current row.
        self.row = 0

        # Open the file in binary mode, even for text files.
        self.f = open(filename, 'rb')

    def next(self):
        """Return the values in the next row."""
        self.f.seek(self.rowsize_bytes * self.row, os.SEEK_SET)
        values = np.fromfile(self.f, dtype=self.dtype, count=self.rowsize)
        self.row += 1
        return values

    def close(self):
        self.f.close()

    def __del__(self):
        self.close()

t0 = time.clock()
# Get the filenames.
folder = r"/media/bigdata/i01_maze05.005"
basename = "i01_maze05_MS.005"

dir = os.path.join(folder, basename)
files = os.listdir(dir)

filenames = {}
for ext in ['xml', 'spk']:
    fileindex_set = set()
    for file in files:
        r = re.search(ext, file)
        if r:
            r = re.search(r"([^\n]+)\.([^\.]+)\.([0-9]+)$", file)
            if r:
                fileindex_set.add(int(r.group(3)))
    filenames[ext] =sorted(fileindex_set)
print filenames

metadata = read_xml(dir + '.xml', 1)

data = MemMappedBinary(dir + '.spk.2', np.int16,
            rowsize=metadata['nchannels'] * metadata['nsamples'])

plt.figure(frameon=False, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
print data.rowsize
for i in range(1000):
    spk = data.next()
    plt.plot(spk)
    print i
plt.show()