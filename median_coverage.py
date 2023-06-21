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

parser = argparse.ArgumentParser(description='')
parser.add_argument("input", \
                    help="Input gzipped tab delimited file of coverage.", \
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


def write_out(outlines):
	''' Write the output.
	'''
	with open(args.output+".txt", "w") as out:
		out.write("\n".join(outlines))
		out.write("\n")

def main():
	outlines = [] # Save each output line as a string in a list
	with gzip.open(args.input, "rt") as infile:
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
				# 2. Save each value to its respective dict entry
				#import ipdb; ipdb.set_trace()
				# 1.
				dat = np.asarray(fields[2:], dtype = int)
				median = np.median(dat)
				sd = np.std(dat)
				outlines.append("\t".join([fields[0], fields[1], str(median), str(sd)]))
				# 2.
				for idx, val in enumerate(fields[2:]):
					cov_per_sample[samples[idx]].append(int(val))

		# After done looping through the file, convert lists in the dict to arrays and calculate median, sd
		print("# Median coverage per sample")
		for k,v in cov_per_sample.items():
			dat = np.asarray(v)
			median, sd = np.median(dat), np.std(dat)
			print("\t".join([k, str(median), str(sd)]))
		# And write the output file.
		write_out(outlines)

if __name__ == "__main__":
    main()
