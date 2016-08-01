# mg-storage-hdf5

Is for the storage and retrieval of interactions from an adjacency matrix stored in HDF5.

# Requirements
- Python 2.7+
- Python Modules:
  - h5py
  - NumPy

# Code Example
## Split the list
Reads through the input adjacency list and splits it into managable chunks. This script is also able to handle interactions between different chromosomes.

```
python SplitAdjacencyFile.py --fin=chr1_1kb.RAWobserved --dataset=I1k --binSize=1000 --chrA=1 --chrB=1
```

## Load HDF5
Loads each of the chunks into the HDF5 file

```
python CreateHDF5.py --hdf5=test.hdf5 --genome=hg37 --dataset=I10k --binSize=10000 --chrA=1 --chrB=20
```
This relies on the presence of a chrom.size file that contains a tab separated list of the chromosome number and the length of the chromosome. This way is is able to calculate the corect size matrix in teh HDF5 file.

## Get data from HDF5
Extract interactions that occur between a start and stop on a chromosome

```
python GetRange.py --hdf5=test.hdf5 --dataset=I1k --binSize=1000 --chr=1 --start=1 --end=2000
```
