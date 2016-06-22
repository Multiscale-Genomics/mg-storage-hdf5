#!/usr/bin/env python

# Load external modules
import numpy as np
import argparse
from datetime import datetime

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--fin", help="Input file")
parser.add_argument("--dataset", help="Name of the dataset in the HDF5 file")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--chrA", help="")
parser.add_argument("--chrB", help="")

# Get the matching parameters from the command line
args = parser.parse_args()
fin = args.fin
dataset = args.dataset
binSize = args.binSize
chrA = args.chrA
chrB = args.chrB

t = datetime.now().time()
print "Loaded params: " + str(t)

maxBinSize = 10000

fi = open(fin, "r")
lCount = 0
lineDict = {}
for line in fi:
  originalLine = line
  line = line.rstrip()
  line = line.split("\t")
  
  x = (int(line[0])/int(binSize))
  y = (int(line[1])/int(binSize))
  v = int(float(line[2]))
  
  xBin = int(np.floor(x/maxBinSize))
  yBin = int(np.floor(y/maxBinSize))
  
  if str(xBin) + "_" + str(yBin) in lineDict:
    lineDict[str(xBin) + "_" + str(yBin)].append(originalLine)
  else:
    lineDict[str(xBin) + "_" + str(yBin)] = [originalLine]
  
  lCount += 1
  
  if lCount % 1000000 == 0:
    for k in lineDict.keys():
      fo = open(dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + k + ".tmp", "a")
      for l in lineDict[k]:
        fo.write(l)
      fo.close()
    lineDict.clear()

fi.close()

for k in lineDict.keys():
  fo = open(dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + k + ".tmp", "a")
  for l in lineDict[k]:
    fo.write(l)
  fo.close()

