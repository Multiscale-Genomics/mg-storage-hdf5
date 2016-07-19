#!/usr/bin/env python

# Load external modules
import h5py
import numpy as np
import argparse
import os.path
import sys
from datetime import datetime

#from hdf5_loader import hdf5
from datasets import datasets

def countArray(d, diag=True):
  dCount = 0
  if diag == True:
    for i in xrange(len(d)):
      dCount += np.count_nonzero(d[i,0:i+1])
    return dCount
  return np.count_nonzero(d)
  

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
#h5 = hdf5()
ds = datasets()
#chr_param = h5.get_chromosome_sizes(genome)
chr_param = ds.getChr_param()

t = datetime.now().time()
print "Loaded params: " + str(t)

# Defines the 2 genomes where A <= B
chrA_meta = ds.getChromosome(genome, int(binSize), chrA)
chrB_meta = ds.getChromosome(genome, int(binSize), chrB)
x_chrA = chrA_meta["bins"]
y_chrB = chrB_meta["bins"]

# Max size of the array (size ofthe whole genome)
xy = chr_param[genome]["meta"]["totalBinCount"][int(binSize)][1]

# Max number of potential bins. Requires pre-splitting the adjacency files
if chrA_meta["bins"] > 10000:
  axy = 10000
else:
  axy = chrA_meta["bins"]

if chrB_meta["bins"] > 10000:
  bxy = 10000
else:
  bxy = chrB_meta["bins"]

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

loading_break_point = 1000000
lCount = 0
total_lCount = 0
t = datetime.now()

# Number of cells the adjacency list was split into for chrA
x_cells = int(np.ceil(float(chrA_meta["bins"])/float(axy)))
for x_cell in xrange(x_cells):
  # Define the start and end points of the array for ChrA
  x_start = x_cell*axy
  x_end   = x_cell*axy + axy
  
  # Number of cells the adjacency list was split into for chrB
  y_cells = int(np.ceil(float(chrB_meta["bins"])/float(bxy)))
  for y_cell in xrange(y_cells):
    if os.path.isfile(dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp") == False:
      continue
    
    # Define the start and end points of the array for ChrB
    y_start = y_cell*bxy
    y_end   = y_cell*bxy + bxy
    
    # Ensure the end x or y value is not greater than the total size of the genome
    if x_end > xy:
      x_end = xy
    if y_end > xy:
      y_end = xy
    
    x_limit = 0
    y_limit = 0
    
    if x_end > chrA_meta["bins"]:
      x_limit = x_end - chrA_meta["bins"]
    if y_end > chrB_meta["bins"]:
      y_limit = y_end - chrB_meta["bins"]
    
    x_offset = ds.getOffset(genome, binSize, chrA)
    y_offset = ds.getOffset(genome, binSize, chrB)
    
    # Ensure that this sub section of the array is set to zero before loading
    x_matrix_start = x_start + x_offset
    x_matrix_end   = x_end   + x_offset - x_limit
    
    y_matrix_start = y_start + y_offset
    y_matrix_end   = y_end   + y_offset - y_limit
    
    dset[x_matrix_start:x_matrix_end, y_matrix_start:y_matrix_end] = np.zeros([axy-x_limit,bxy-y_limit], dtype='int32')
    if x_matrix_start != y_matrix_start:
      dset[y_matrix_start:y_matrix_end, x_matrix_start:x_matrix_end] = np.zeros([bxy-y_limit,axy-x_limit], dtype='int32')
    ###
    
    d1 = np.zeros([axy-x_limit,bxy-y_limit], dtype='int32')
    d2 = np.zeros([bxy-y_limit,axy-x_limit], dtype='int32')
    
    file_name = dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp"
    print ""
    print "  File Ready: " + file_name
    print "    x-coords: x_cell:" + str(x_cell) + ", x_start:" + str(x_start) + ", x_end:" + str(x_end) + ", x_limit:" + str(x_limit)
    print "    y-coords: y_cell:" + str(y_cell) + ", y_start:" + str(y_start) + ", y_end:" + str(y_end) + ", y_limit:" + str(y_limit)
    
    #continue
    
    fi = open(file_name, "r")
    
    t0 = datetime.now()
    
    lCount = 0
    for line in fi:
      line = line.rstrip()
      line = line.split("\t")
      
      x = (int(line[0])/int(binSize))
      y = (int(line[1])/int(binSize))
      v = int(float(line[2]))
      
      if v == 0:
        print "WARNING: VALUE IS ZERO"
        print x, x_start, x_end, y, y_start, y_end, v
        sys.exit(1)
      
      if x>=x_start and x<x_end and y>=y_start and y<y_end:
        d1[x-x_start, y-y_start] = v
        if chrA == chrB and x_cell == y_cell:
          d1[y-y_start, x-x_start] = v
        else:
          d2[y-y_start, x-x_start] = v
      else:
        print "ERROR: MIS-PLACED VALUE: " + str(line)
        print x, x_start, x_end, y, y_start, y_end
        exit(1)
      
      if d1[x-x_start, y-y_start] != v:
        print "ERROR: INSERTED VALUE DOES NOT MATCH RETRIEVED VALUE"
        print x, x_start, x_end, y, y_start, y_end, v, d1[x-x_start, y-y_start]
        sys.exit(1)
      
      lCount+=1
      
      if lCount % loading_break_point == 0:
        
        x_matrix_start = x_start + x_offset
        x_matrix_end   = x_end   + x_offset - x_limit
        
        y_matrix_start = y_start + y_offset
        y_matrix_end   = y_end   + y_offset - y_limit
        
        if x_matrix_start != y_matrix_start:
          d1Count = countArray(d1, False)
        else:
          d1Count = countArray(d1)
        
        if d1Count != loading_break_point:
          print "ERROR: VALUES NOT LOADED"
          print np.nonzero(d1), np.count_nonzero(d1)
          print "x:", x, "y:", y, "v:", v
          print "x_start:", x_start, "x_end:", x_end
          print "y_start:", y_start, "y_end:", y_end
          print "d1Count:", d1Count, "lCount:", lCount
          print d1[x-x_start, y-y_start]
          sys.exit(1)
        
        
        
        dset[x_matrix_start:x_matrix_end, y_matrix_start:y_matrix_end] += d1
        if x_matrix_start != y_matrix_start:
          dset[y_matrix_start:y_matrix_end, x_matrix_start:x_matrix_end] += d2
        d1 = np.zeros([axy-x_limit,bxy-y_limit], dtype='int32')
        d2 = np.zeros([bxy-y_limit,axy-x_limit], dtype='int32')
        print "    ... " + str(lCount)
        
        
      
    fi.close()
    
    x_matrix_start = x_start + x_offset
    x_matrix_end   = x_end  + x_offset - x_limit
    
    y_matrix_start = y_start + y_offset
    y_matrix_end   = y_end  + y_offset - y_limit
    
    dset[x_matrix_start:x_matrix_end, y_matrix_start:y_matrix_end] += d1
    diag=True
    if x_matrix_start != y_matrix_start:
      diag=False
      dset[y_matrix_start:y_matrix_end, x_matrix_start:x_matrix_end] += d2
    
    print "    ... " + str(lCount)
    #print "      d1Count:    " + str(d1Count)
    
    t1 = datetime.now()
    
    #stopCount=countArray(dset)
    #d1Count = countArray(d1)
    
    #total_lCount = lCount + startCount
    
    #if stopCount != total_lCount:
    #    print ""
    #    print "ERROR: MISSING VALUES"
    #    print "  startCount: " + str(startCount)
    #    print "  stopCount:  " + str(stopCount)
    #    print "  d1Count:    " + str(d1Count)
    #    print "  lCount:     " + str(lCount)
    #    print "  totalCount: " + str(total_lCount)
    #    print "  Coordinates:"
    #    print "    x_matrix_start: " + str(x_matrix_start)
    #    print "    x_matrix_end:   " + str(x_matrix_end)
    #    print "    y_matrix_start: " + str(y_matrix_start)
    #    print "    y_matrix_end:   " + str(y_matrix_end)
    #    sys.exit(1)
     
    print "  Parsed: " + dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp - Run time: " + str(t1-t0)
    print "    Loaded " + str(lCount) + " lines"

f.close()

