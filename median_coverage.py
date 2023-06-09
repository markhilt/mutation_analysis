#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
median_coverage.py
Author: Markus Hiltunen
E-mail: markus.hiltunen@su.se

This script is used to calculate the median coverage per sample from
a file containing the sequencing coverage for each base in a reference
genome.

Copyright (c) 2023, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

import argparse
import numpy as np
import gzip

__version__ = "0.1"

parser = argparse.ArgumentParser(description='	Calculates the median coverage both per sample \
												and per base pair, from a file containing the \
												sequencing coverage for each base in a reference \
												genome. The coverage file is gzipped and can be \
												the output from e.g. samtools depth. WARNING: \
												script may use a lot of memory, depending on \
												genome size and number of samples.')
parser.add_argument("input", \
                    help="Input gzipped tab delimited file of depth of coverage.\n \
                    	  Format is: 'scaffold\tcoordinate\tdepth_sample_1\t...\tdepth_sample_n' .\n \
                    	  The first line is a header in the format: '# CHROM\tPOS\tname_sample_1\t...\tname_sample_n'", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Output prefix.", \
                    default="median_coverage_out", \
                    type = str)
parser.add_argument("-v","--version", \
                    help="Print version and quit.", \
                    action = "version", \
                    version = "median_coverage v.{}".format(__version__))
args = parser.parse_args()

def main():
	with gzip.open(args.input, "rt") as infile:
		with open(args.output+".txt", "w") as outfile:
			for line in infile:
				line = line.strip()
				fields = line.split("\t")

				# First line is the header. Here we define the samples
				# Order is important, since the columns are in the same order
				if line.startswith("#"):
					# Keep the samples in a list, so we know the order
					samples = [ f for f in fields[2:] ]
					# Also define a dict, which is unordered, but we will use the list to fill it.
					cov_per_sample = { f:[] for f in fields[2:] }

				else:
					# For each data line we want to do two things
					# 1. Calculate the median for the line = genomic position
					dat = np.asarray(fields[2:], dtype = int)
					median = np.median(dat)
					sd = np.std(dat)
					outfile.write("{}\t{}\t{}\t{}\n".format(fields[0], fields[1], str(median), str(sd)))
					# 2. Save each value to its respective dict entry. Use enumerate to get the index for the sample
					for idx, val in enumerate(fields[2:]):
						cov_per_sample[samples[idx]].append(int(val))

	# After done looping through the file, convert lists in the dict to arrays and calculate median, sd
	print("# Median coverage per sample")
	for k,v in cov_per_sample.items():
		dat = np.asarray(v)
		median, sd = np.median(dat), np.std(dat)
		print("\t".join([k, str(median), str(sd)]))

if __name__ == "__main__":
    main()
