#!/usr/bin/env python

# Load external modules
import h5py
import numpy as np
import argparse
import os.path
from datetime import datetime

from hdf5_loader import hdf5

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--genome", help="Genome name")
parser.add_argument("--dataset", help="Name of the dataset")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--chrA", help="Chromosome A")
parser.add_argument("--chrB", help="Chromosome B")


# Get the matching parameters from the command line
args = parser.parse_args()
hdf5_file = args.hdf5
genome = args.genome
dataset = args.dataset
binSize = args.binSize
chrA = args.chrA
chrB = args.chrB

# Load the chromosome sizes for the whole genome
h5 = hdf5()
chr_param = h5.get_chromosome_sizes(genome)
"""
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
"""



t = datetime.now().time()
print "Loaded params: " + str(t)

# Defines the 2 genomes where A <= B
x_chrA = chr_param[genome][chrA]["bins"][int(binSize)][0]
y_chrB = chr_param[genome][chrB]["bins"][int(binSize)][0]

# Max size of the array (size ofthe whole genome)
xy = chr_param[genome]["meta"]["totalBinCount"][int(binSize)][1]

#print chr_param[genome][chrA]
#print chr_param[genome][chrB]

if x_chrA > 10000:
  axy = 10000
else:
  axy = x_chrA

if y_chrB > 10000:
  bxy = 10000
else:
  bxy = y_chrB

tmpCount = 0

d1 = np.zeros([axy,bxy], dtype='int32')
d2 = np.zeros([bxy,axy], dtype='int32')

# Open the files. The HDF5 file is created is it does not already exist
f = h5py.File(hdf5_file, "a")

if binSize in f:
  # raise NameError("Dataset with name " + dataset + " already exists.")
  dset = f[binSize]
else:
  dset = f.create_dataset(binSize, (xy, xy), dtype='int32', chunks=True, compression="gzip")

#print axy, bxy, xy

t = datetime.now().time()
t1 = datetime.now()
print "Created dataset: " + str(t)

lCount = 0
t = datetime.now()

x_cells = int(np.ceil(float(x_chrA)/float(axy)))
for x_cell in xrange(x_cells):
  x_start = x_cell*axy
  x_end   = x_cell*axy + axy
  
  y_cells = int(np.ceil(float(y_chrB)/float(axy)))
  for y_cell in xrange(y_cells):
    if os.path.isfile(dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp") == False:
      continue
    
    y_start = y_cell*bxy
    y_end   = y_cell*bxy + bxy
    
    if x_end > xy:
      x_end = xy
    if y_end > xy:
      y_end = xy
    
    x_limit = 0
    y_limit = 0
    
    if x_end > chr_param[genome][chrA]["bins"][int(binSize)][0]:
      x_limit = x_end - chr_param[genome][chrA]["bins"][int(binSize)][0]
    if y_end > chr_param[genome][chrB]["bins"][int(binSize)][0]:
      y_limit = y_end - chr_param[genome][chrB]["bins"][int(binSize)][0]
    
    d1 = np.zeros([axy-x_limit,bxy-y_limit], dtype='int32')
    d2 = np.zeros([bxy-y_limit,axy-x_limit], dtype='int32')
    
    print ""
    print "  File Ready: " + dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp"
    print "    x-coords: x_cell:" + str(x_cell) + ", x_start:" + str(x_start) + ", x_end:" + str(x_end) + ", x_limit:" + str(x_limit)
    print "    y-coords: y_cell:" + str(y_cell) + ", y_start:" + str(y_start) + ", y_end:" + str(y_end) + ", y_limit:" + str(y_limit)
    
    #continue
    
    fi = open(dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp", "r")
    t0 = datetime.now()
    
    print "  Iterator start"
    lCount = 0
    for line in fi:
      line = line.rstrip()
      line = line.split("\t")
      
      x = (int(line[0])/int(binSize))
      y = (int(line[1])/int(binSize))
      v = int(float(line[2]))
      
      if x>=x_start and x<x_end and y>=y_start and y<y_end:
        d1[x-x_start, y-y_start] = v
        d2[y-y_start, x-x_start] = v
      else:
        print "Mis-placed value: " + str(line)
        print x, x_start, x_end, y, y_start, y_end
        exit(0)
      
      lCount+=1
      if lCount % 1000000 == 0:
        #x_offset = int(np.ceil(float(chr_param[genome][str(chrA)][1] - chr_param[genome][str(chrA)][0])/float(binSize)))
        x_offset = chr_param[genome][chrA]["bins"][int(binSize)][1]
        x_matrix_start = x_start + x_offset
        x_matrix_end   = x_end  + x_offset - x_limit
        
        #y_offset = int(np.ceil(float(chr_param[genome][str(chrB)][1] - chr_param[genome][str(chrB)][0])/float(binSize)))
        y_offset = chr_param[genome][chrB]["bins"][int(binSize)][1]
        y_matrix_start = y_start + y_offset
        y_matrix_end   = y_end  + y_offset - y_limit
        
        dset[x_matrix_start:x_matrix_end, y_matrix_start:y_matrix_end] += d1
        if x_start != y_start:
          dset[y_matrix_start:y_matrix_end, x_matrix_start:x_matrix_end] += d2
        d1 = np.zeros([axy-x_limit,bxy-y_limit], dtype='int32')
        d2 = np.zeros([bxy-y_limit,axy-x_limit], dtype='int32')
        print "    ... " + str(lCount)
    fi.close()
    
    #x_offset = int(np.ceil(float(chr_param[str(chrA)][2] - chr_param[str(chrA)][1])/float(binSize)))
    x_offset = chr_param[genome][chrA]["bins"][int(binSize)][1]
    x_matrix_start = x_start + x_offset
    x_matrix_end   = x_end  + x_offset - x_limit
    
    #y_offset = int(np.ceil(float(chr_param[str(chrB)][2] - chr_param[str(chrB)][1])/float(binSize)))
    y_offset = chr_param[genome][chrB]["bins"][int(binSize)][1]
    y_matrix_start = y_start + y_offset
    y_matrix_end   = y_end  + y_offset - y_limit
    
    dset[x_matrix_start:x_matrix_end, y_matrix_start:y_matrix_end] += d1
    if x_start != y_start:
      dset[y_matrix_start:y_matrix_end, x_matrix_start:x_matrix_end] += d2
    
    t1 = datetime.now()
    print "  Parsed: " + dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp - Run time: " + str(t1-t0)

f.close()

