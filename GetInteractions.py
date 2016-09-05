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
import argparse
import operator

from datasets import datasets
from hdf5_reader import hdf5

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
#parser.add_argument("--hdf5", help="HDF5 file name")
parser.add_argument("--taxon", help="Taxon ID")
parser.add_argument("--accession", help="Genome accession")
parser.add_argument("--dataset", help="Name of the dataset in the HDF5 file")
parser.add_argument("--binSize", help="Size of the bins (bp)")
parser.add_argument("--chr", help="Chromosome")
parser.add_argument("--start", help="Chromosome Start (bp)")
parser.add_argument("--end", help="Chromosome End (bp)")
parser.add_argument("--limit_chr", help="Limit the second chromosome", default=None)
parser.add_argument("--limit_region", help="Limit the to inter or intra chromosomal interactions", default=None)

# Get the matching parameters from the command line
args = parser.parse_args()
#hdf5 = args.hdf5
taxon_id = args.taxon
accession_id = args.accession
dataset = args.dataset
resolution = args.binSize
chr_id = args.chr
start = float(args.start)
end = args.end
limit_chr = args.limit_chr
limit_region = args.limit_region

ds = datasets()
h5 = hdf5()
x = h5.get_range(ds, dataset, resolution, accession_id, chr_id, start, end, limit_region, limit_chr)

#print x["log"]
for v in x["results"]:
    print str(v["chrA"]) + "\t" + str(v["startA"]) + "\t" + str(v["chrB"]) + "\t" + str(v["startB"]) + "\t" + str(v["value"])

