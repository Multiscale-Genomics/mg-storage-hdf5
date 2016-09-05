#!/usr/bin/env python

"""
Copyright 2016 EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# Load external modules
import h5py
import numpy as np
import argparse

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--dataset", help="Name of the dataset in the HDF5 file")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--chr", help="Chromosome")
parser.add_argument("--start", help="Chromosome Start (bp)")
parser.add_argument("--end", help="Chromosome End (bp)")

# Get the matching parameters from the command line
args = parser.parse_args()
hdf5 = args.hdf5
dataset = args.dataset
binSize = args.binSize
chromo = args.chr
start = float(args.start)
end = args.end

# Load the chromosome sizes for the whole genome
chrSizes = open("chrom.size", "r")
chr_param_list = []
chr_param = {}
genomeLen = 0
for line in chrSizes:
  line = line.rstrip()
  line = line.split("\t")
  
  c = line[0]
  l = int(line[1])
  
  genomeLen += l
  
  chr_param_list.append([c, l, genomeLen])
  chr_param[c] = [c, l, genomeLen]
chrSizes.close()

# Defines columns to get extracted from the array
x = int(np.floor(float(start)/float(binSize)))
y = int(np.ceil(float(end)/float(binSize)))

xy_offset = chr_param[c][2]-chr_param[c][1]

# Open the hdf5 file
f = h5py.File(hdf5, "r")
dset = f[dataset]
print dset[(x+xy_offset):(y+xy_offset),:]

f.close()


