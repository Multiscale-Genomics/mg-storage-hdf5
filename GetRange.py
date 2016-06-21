#!/usr/bin/env python

# Load external modules
import h5py
import numpy as np
import argparse

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--dataset", help="Name of the dataset in the HDF5 file")
parser.add_argument("--chrLen", help="Chromosome length (bp)")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--start", help="Chromosome Start (bp)")
parser.add_argument("--end", help="Chromosome End (bp)")

# Get the matching parameters from the command line
args = parser.parse_args()
hdf5 = args.hdf5
dataset = args.dataset
chrLen = args.chrLen
binSize = args.binSize
start = float(args.start)
end = args.end

# Defines columns to get extracted from the array
x = int(np.floor(float(start)/float(binSize)))
y = int(np.ceil(float(end)/float(binSize)))

# Open the hdf5 file
f = h5py.File(hdf5, "r")
dset = f[dataset] 
print dset[x:y,:]

f.close()


