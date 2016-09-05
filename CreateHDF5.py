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
import numpy as np
import argparse, os.path, sys
from datetime import datetime

from hdf5_loader import hdf5
from datasets import datasets


# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--taxon", help="Taxon ID")
parser.add_argument("--genome", help="Genome name")
parser.add_argument("--dataset", help="Name of the dataset")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--chrA", help="Chromosome A")
parser.add_argument("--chrB", help="Chromosome B")


# Get the matching parameters from the command line
args = parser.parse_args()
hdf5_file = args.hdf5
taxon_id = args.taxon
genome = args.genome
dataset = args.dataset
binSize = args.binSize
chrA = args.chrA
chrB = args.chrB

# Load the chromosome sizes for the whole genome
h5l = hdf5()
ds = datasets()
chr_param = ds.getChr_param()

t = datetime.now().time()
print "Loaded params:                " + str(t)

# Defines the 2 genomes where A <= B
chrA_meta = ds.getChromosome(genome, int(binSize), chrA)
chrB_meta = ds.getChromosome(genome, int(binSize), chrB)

# Max number of potential bins. Requires pre-splitting the adjacency files
if chrA_meta["bins"] > 10000:
    axy = 10000
else:
    axy = chrA_meta["bins"]

if chrB_meta["bins"] > 10000:
    bxy = 10000
else:
    bxy = chrB_meta["bins"]

t = datetime.now().time()
print "Start loading adjacency data: " + str(t)

# Number of cells the adjacency list was split into for chrA
x_cells = int(np.ceil(float(chrA_meta["bins"])/float(axy)))
for x_cell in xrange(x_cells):
  
    # Number of cells the adjacency list was split into for chrB
    y_cells = int(np.ceil(float(chrB_meta["bins"])/float(bxy)))
    for y_cell in xrange(y_cells):
        file_name = dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp"
        if os.path.isfile(file_name) == False:
            continue
    
    h5l.load_area(ds, genome, dataset, file_name, binSize, chrA, chrB, x_cell, y_cell)


