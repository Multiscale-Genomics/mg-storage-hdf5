# mg-storage-hdf5

ls for the storage and retrieval of interactions from an adjacency matrix stored in HDF5.

# Requirements
- Python 2.7+
- Python Modules:
  - h5py
  - NumPy

# Code Example
## Split the list
Reads through the input adjacency list and splits it into managable chunks

```
python SplitAdjacencyFile.py  --fin=chr1_1kb.RAWobserved --dataset=I1k --binSize=1000
```

## Load HDF5
Loads each of the chunks into the HDF5 file

```
python CreateHDF5.py  --hdf5=test.hdf5 --dataset=I10k --chrLen=249250621 --binSize 10000
```

