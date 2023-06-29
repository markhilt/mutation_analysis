#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_variant_fasta.py
Author: Markus Hiltunen
E-mail: markus.hiltunen@su.se

Extract the sequence flanking a variant.

Description: Extract the sequence in fasta format that
is directly flanking a variant of interest. Output
is a fasta file with variant bases in IUPAC coding
and indels masked as multi-N.

Copyright (c) 2023, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import vcfpy
import pandas as pd

__version__ = "0.1"

parser = argparse.ArgumentParser(description='  Extract the sequence in fasta format that \
                                                is directly flanking a variant of interest. \
                                                Output is a fasta file with variant bases in \
                                                IUPAC coding and indels masked as multi-N.')
parser.add_argument("reference", \
                    help="Input fasta file used as reference during variant calling.", \
                    type = str)
parser.add_argument("-c","--vcf", \
                    help="VCF file with variants.", \
                    type = str)
parser.add_argument("-t","--table", \
                    help="  Tab delimited file of variants of interest. \
                            First two columns are read, and should be \
                            formatted a VCF file (also accepted directly).", \
                    type = str)
parser.add_argument("-f","--flanks", \
                    help="Number of flanking bp to output [1000].", \
                    type = int,
                    default = 1000)
parser.add_argument("-o","--output", \
                    help="Output prefix.", \
                    type = str)
parser.add_argument("-v","--version", \
                    help="Print version and quit.", \
                    action = "version", \
                    version = "extract_variant_fasta v.{}".format(__version__))
args = parser.parse_args()

def getOut():
    '''
    To create a prefix for output files.
    '''
    if args.output:
        outprefix = args.output
    else:
        outprefix = "output"
    return outprefix

def readVariants(tablefile):
    ''' Read the table of variants and return it as a list of tuples.
    '''
    variants_list = [] # List of tuples
    with open(tablefile, "r") as table:
        for line in table:
            line = line.strip()
            if not line.startswith("#"):
                ctg,crd = line.split("\t")[0],int(line.split("\t")[1])
                variants_list.append((ctg,crd))
    return variants_list

def extractSeq(fasta, variant, flank):
    ''' Extract the sequence surrounding variant.
    '''
    ctg,crd = variant[0], variant[1] # Position of the variant
    seq = fasta[ctg][crd-flank:crd+flank] # Extract the sequence
    return seq.seq

def iupacify():
    ''' Create a pandas dataframe for the IUPAC coding.
    R   Purine (A or G)
    Y   Pyrimidine (C or T)
    M   C or A
    K   T or G
    W   T or A
    S   C or G
    '''
    iupac_df = pd.DataFrame(None, index=["A", "T", "C", "G"], columns=["A", "T", "C", "G"])
    iupac_df["A"]["G"] = "R"
    iupac_df["A"]["C"] = "M"
    iupac_df["A"]["T"] = "W"
    iupac_df["C"]["G"] = "S"
    iupac_df["C"]["A"] = "M"
    iupac_df["C"]["T"] = "Y"
    iupac_df["G"]["A"] = "R"
    iupac_df["G"]["T"] = "K"
    iupac_df["G"]["C"] = "S"
    iupac_df["T"]["A"] = "W"
    iupac_df["T"]["C"] = "Y"
    iupac_df["T"]["G"] = "K"
    # Fill in remaining NaNs
    for base in ["A", "T", "C", "G"]:
        iupac_df[base][base] = base
    return iupac_df


def ambiguify(pos, seq, variants, flank):
    ''' Swap out bases from seq at the positions in variants
        for IUPAC ambiguities.
    '''
    ctg,crd = pos[0], pos[1] # Position of the variant of interest
    rel_start = pos[1]-flank # Relative start position in the sequence
    ambseq = list(seq.upper()) # Make a list to mutate
    iupac = iupacify()
    for var in variants:
        variant_pos = var.POS - rel_start - 1 # relative position of variant
        if len(var.ALT) > 1: # Keep only biallelic sites
            continue
        zygosity = "hom" if var.INFO['AF'][0] == 1.0 else "het" # Track homo/heterozygosity
        ref_allele = ambseq[variant_pos]
        if var.is_snv():
            if zygosity == "hom":
                # If vcf is homozygous for alternate allele, replace with that base
                ambseq[variant_pos] = var.ALT[0].value
            else:
                # Else we replace with the IUPAC coded base.
                #import ipdb; ipdb.set_trace()
                ambiguous = iupac[ambseq[variant_pos]][var.ALT[0].value]
                ambseq[variant_pos] = ambiguous
        else:
        # If variant is indel, but sample is fixed for alternate allele,
        # we can insert the alternate allele
            if zygosity == "hom":
                ambseq[variant_pos] = var.ALT[0].value
            else:
                ambseq[variant_pos] = "N"
        # Finally, check if the variant is the one of interest,
        # and if so, add [] around it
        if var.POS == crd:
            ambseq[variant_pos] = "["+ambseq[variant_pos]+"]"
    ambseq = "".join(ambseq)
    return ambseq 

def main():
    outfile = getOut() + ".fasta"
    fasta_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    outfasta = {}
    variants = readVariants(args.table)
    reader = vcfpy.Reader.from_path(args.vcf)
    for variant in variants:
        ctg,start,end = variant[0], variant[1]-args.flanks, variant[1]+args.flanks
        seq = extractSeq(fasta_dict, variant, args.flanks)
        vcf_subset = reader.fetch(ctg, begin=start, end=end)
        variantseq = ambiguify(variant, seq, vcf_subset, args.flanks)
        # The output fasta header:
        header = "{}:{}-{}".format(ctg, str(start), str(end))
        outfasta[header] = variantseq
    reader.close()

    with open(outfile, "w") as out:
        for k,v in outfasta.items():
            out.write(">{}\n{}\n".format(k,v))

if __name__ == "__main__":
    main()
