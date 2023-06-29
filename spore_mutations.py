#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spore_mutations.py
Author: Markus Hiltunen
E-mail: markus.hiltunen@su.se

Copyright (c) 2023, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

import argparse

__version__ = "0.1"

parser = argparse.ArgumentParser(description='')
parser.add_argument("input", \
                    help="Input vcf.", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Output prefix.", \
                    type = str)
parser.add_argument("-v","--version", \
                    help="Print version and quit.", \
                    action = "version", \
                    version = "spore_mutations v.{}".format(__version__))
args = parser.parse_args()

def main():
    with open(args.input, "r") as vcffile:
        for line in vcffile:
            line = line.strip()
            fields = line.split("\t")
            if line.startswith("##"):
                print(line)
            elif line.startswith("#CHROM"):
                # Divide the input based on prefix of sample names.
                # Samples start at column 10
                # Build an index of in which fields different fruiting bodies occur
                sample_fields = {}
                for idx,sample in enumerate(fields[10:]):
                    sample_fields[idx] = sample[0]
                print(line)
            else:
                # Now we look for entries in the vcf file where we have variant calls from exclusively one
                # fruiting body.
                variant_fruitbodies = set([sample_fields[idx] for idx, gt in enumerate(fields[10:]) if "/1" in gt or "|1" in gt])
                if not len(variant_fruitbodies) > 1:
                    print(line)

if __name__ == "__main__":
    main()
