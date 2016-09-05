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

import h5py
import numpy as np
import os.path
import sys
from datetime import datetime

class hdf5:

    def __init__(self):
        """
        Initialise the chr_param as an empty index when the service is started
        """
        self.chr_param = {}
    
    def countArray(d, diag=True):
        """
        Function to analysis the number of values in an array
        """
        
        dCount = 0
        if diag == True:
            for i in xrange(len(d)):
                dCount += np.count_nonzero(d[i,0:i+1])
            return dCount
        return np.count_nonzero(d)
    
    def load_area(self, ds, genome, dataset, file_name, binSize, chrA, chrB, x_cell = 0, y_cell = 0):
        """
        Function to load a sub-dataset into an HDF5 files
        """
        
        binSizeInv = 1.0/float(binSize)
        loading_break_point = 1000000
        
        # Defines the 2 genomes where A <= B
        chrA_meta = ds.getChromosome(genome, int(binSize), chrA)
        chrB_meta = ds.getChromosome(genome, int(binSize), chrB)
        
        # Max size of the array (size ofthe whole genome)
        xy = ds.getTotalBinCount(genome, int(binSize))
        
        # Max number of potential bins. Requires pre-splitting the adjacency files
        if chrA_meta["bins"] > 10000:
            axy = 10000
        else:
            axy = chrA_meta["bins"]

        if chrB_meta["bins"] > 10000:
            bxy = 10000
        else:
            bxy = chrB_meta["bins"]
        
        resource_package = __name__
        resource_path = os.path.join(os.path.dirname(__file__), dataset + '.hdf5')
        f = h5py.File(resource_path, "a")
        if binSize in f:
            dset = f[binSize]
        else:
            dset = f.create_dataset(binSize, (xy, xy), dtype='int32', chunks=True, compression="gzip")
        
        # Define the start and end points of the array for ChrA
        x_start = x_cell*axy
        x_end   = x_cell*axy + axy
        
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
        
        fi = open(file_name, "r")
        
        t0 = datetime.now()
        lCount = 0
        for line in fi:
            line = line.rstrip()
            line = line.split("\t")
            
            x = int(int(line[0]) * binSizeInv)
            y = int(int(line[1]) * binSizeInv)
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
         
        print "  Parsed: " + dataset + "/" + "chr" + str(chrA) + "_chr" + str(chrB) + "_" + str(x_cell) + "_" + str(y_cell) + ".tmp"
        print "    Run time:                 " + str(t1-t0)
        print "    Loaded " + str(lCount) + " lines"
        
        f.close()
        return {"lines": lCount}
