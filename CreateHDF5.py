#!/usr/bin/env python

# Load external modules
import h5py
import numpy as np
import argparse
import os.path
from datetime import datetime

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--dataset", help="Name of the dataset in the HDF5 file")
parser.add_argument("--chrLen", help="Chromosome length (bp)")
parser.add_argument("--binSize", help="Size of the bins (bp)")

# Get the matching parameters from the command line
args = parser.parse_args()
hdf5 = args.hdf5
dataset = args.dataset
chrLen = args.chrLen
binSize = args.binSize

t = datetime.now().time()
print "Loaded params: " + str(t)

# Open the files. The HDF5 file is created is it does not already exist
f = h5py.File(hdf5, "a")

if dataset in f:
  raise NameError("Dataset with name " + dataset + " already exists.")

xy = int(np.ceil(float(chrLen)/float(binSize)))
if xy > 10000:
  bxy = 10000
else:
  bxy = xy

tmpCount = 0

d1 = np.zeros([bxy,bxy], dtype='int32')
d2 = np.zeros([bxy,bxy], dtype='int32')
dset = f.create_dataset(dataset, (xy, xy), dtype='int32', chunks=True, compression="gzip")

t = datetime.now().time()
t1 = datetime.now()
print "Created dataset: " + str(t)

lCount = 0
t = datetime.now()

cells = int(np.ceil(float(xy)/float(bxy)))
for x_cell in xrange(cells):
  x_start = x_cell*bxy
  x_end   = x_cell*bxy + bxy
  for y_cell in xrange(x_cell,cells):
    y_start = y_cell*bxy
    y_end   = y_cell*bxy + bxy
    
    if x_end > xy:
      x_end = xy
    if y_end > xy:
      y_end = xy
    
    x_limit = 0
    y_limit = 0
    if x_start+bxy > xy:
      x_limit = x_start + bxy - xy
    if y_start+bxy > xy:
      y_limit = y_start + bxy - xy
    
    print ""
    print "  x-coords: x_cell:" + str(x_cell) + ", x_start:" + str(x_start) + ", x_end:" + str(x_end) + ", x_limit:" + str(x_limit)
    print "  y-coords: y_cell:" + str(y_cell) + ", y_start:" + str(y_start) + ", y_end:" + str(y_end) + ", y_limit:" + str(y_limit)
    
    d1 = np.zeros([bxy-x_limit,bxy-y_limit], dtype='int32')
    d2 = np.zeros([bxy-y_limit,bxy-x_limit], dtype='int32')
    if os.path.isfile(dataset + "/" + str(x_cell) + "_" + str(y_cell) + ".tmp") == False:
      continue
    fi = open(dataset + "/" + str(x_cell) + "_" + str(y_cell) + ".tmp", "r")
    t0 = datetime.now()
    
    print "  Iterator start"
    lCount = 0
    for line in fi:
      line = line.rstrip()
      line = line.split("\t")
      
      x = (int(line[0])/int(binSize)) - 1
      y = (int(line[1])/int(binSize)) - 1
      v = int(float(line[2]))
      
      if x>=x_start and x<x_end and y>=y_start and y<y_end:
        d1[x-x_start, y-y_start] = v
        d2[y-y_start, x-x_start] = v
      else:
        print "Miss-placed value: " + str(line)
        exit(0)
      
      lCount+=1
      if lCount % 1000000 == 0:
        dset[x_start:x_end, y_start:y_end] += d1
        if x_start != y_start:
          dset[y_start:y_end, x_start:x_end] += d2
        d1 = np.zeros([bxy-x_limit,bxy-y_limit], dtype='int32')
        d2 = np.zeros([bxy-y_limit,bxy-x_limit], dtype='int32')
        print "    ... " + str(lCount)
    fi.close()
    
    dset[x_start:x_end, y_start:y_end] += d1
    if x_start != y_start:
      dset[y_start:y_end, x_start:x_end] += d2
    
    t1 = datetime.now()
    print "  Parsed " + str(x_cell) + " : " + str(y_cell) + " - Run time: " + str(t1-t0)

f.close()

