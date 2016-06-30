import h5py
import numpy as np

class hdf5:
    
    def __init__(self):
        # Initialise the chr_param as an empty index when the service is started
        self.chr_param = {}
    
    """
    Load the chromosomes for a given genomes along with the potential bins and 
    offsets for each bin.
    """
    def load_chromosome_sizes(self, genome):
        # Set the base bin sizes and initial offsets
        binSizes = [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
        binCount = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        # Load the chromosome sizes for the whole genome
        chrSizes = open(genome + ".size", "r")
        genomeLen = 0
        for line in chrSizes:
            line = line.rstrip()
            line = line.split("\t")
            
            c = line[0]
            l = int(line[1])
            
            genomeLen += l
            
            # Calculate the number of bins for a chromosome and then join with
            # the offset values for teh start in the array
            binS = [int(np.ceil(l/float(y))) for y in binSizes]
            binC = dict(zip(binSizes, [[binS[i], binCount[i]]for i in range(len(binCount))]))
            if self.chr_param.has_key(genome):
                self.chr_param[genome][c] = {'size': [l, genomeLen], 'bins': binC}
            else:
                self.chr_param[genome] = {c: {'size':[l, genomeLen], 'bins': binC}}
            
            # Calculate the new offset values.
            binCount = [binCount[i]+binS[i] for i in range(len(binCount))]
        chrSizes.close()
        
        totalBinCount = dict(zip(binSizes, [[binS[i], binCount[i]]for i in range(len(binCount))]))
        self.chr_param[genome]["meta"] = {"genomeSize": genomeLen, "totalBinCount": totalBinCount}
    
    def get_chromosome_sizes(self, genome):
        # Check that the genome indexes are loaded
        if self.chr_param.has_key(genome) == False:
             self.load_chromosome_sizes(genome)
        return self.chr_param
