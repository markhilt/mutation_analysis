#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
call_mutations.py
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

This script is used to filter out background variants and retain potentially
new mutations in samples originating from different tissues of the same individual.
It takes as input a variant table produced by GATK variantsToTable.
The principle is to look for variants that are present
in a subset of tissues. Variants that conform to this principle can be optionally
further filtered, as there might be many false positives left, depending on
how the user filtered the inital vcf file.

Copyright (c) 2022, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

from operator import truediv # To be able to divide lists
from scipy import stats as ss
import argparse

__version__ = "0.3"

parser = argparse.ArgumentParser(description='Filter variants to find potential \
                                            mutations. Input: variant table from \
                                            GATK VariantsToTable. The first 5 \
                                            columns contain data about the variant, \
                                            data for each sample start at the 6th. \
                                            Each element in the table is comma separated \
                                            and contains number of reads with reference \
                                            allele and number of reads with alternate allele.')
parser.add_argument("input", \
                    help="Input table.", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Output prefix.", \
                    type = str)
parser.add_argument("-n","--minimum_coverage", \
                    help="Minimum coverage to keep a variant [20].", \
                    type = int, \
                    default = 20)
parser.add_argument("-x","--maximum_coverage", \
                    help="Maximum coverage to keep a variant [60].", \
                    type = int, \
                    default = 60)
parser.add_argument("-r","--minimum_reads", \
                    help="Minimum number of reads to keep a variant. At least \
                    one sample needs to have this many reads to keep the variant [2].", \
                    type = int, \
                    default = 2)
parser.add_argument("-f","--frequency_cutoff", \
                    help="Minimum variant read frequency to keep a variant [0.2].", \
                    type = float, \
                    default = 0.2)
parser.add_argument("-s","--skip_indels", \
                    help="Add to skip indels (useful for pileup data).", \
                    action='store_true')
parser.add_argument("-b","--background", \
                    help="Add to quit after filtering out background variants.", \
                    action='store_true')
parser.add_argument("-t","--test_binomial", \
                    help="Add to test the reads for a binomial distribution. \
                    Only reads from the sample(s) containing the variant are considered.", \
                    action='store_true')
parser.add_argument("-v","--version", \
                    help="Print version and quit.", \
                    action = "version", \
                    version = "call_mutations v.{}".format(__version__))
args = parser.parse_args()

def getOut():
    '''
    To create a prefix for output files.
    '''
    if args.output:
        outprefix = args.output
    else:
        outprefix = "potential_mutations"
    return outprefix

def test_signif(var,cov):
    '''
    Tests for significant skew of reads towards one variant
    Allows one error per 100 000 variants
    True means site passed check (= is not significantly different
    from a 0.5 frequency)
    '''
    alpha = 0.00001
    p = ss.binom_test(var,cov,1.0/2,alternative="two-sided")
    if p < alpha:
        return False
    else:
        return True

def main():
    outfile = getOut() + ".table"
    outlines = []

    with open(args.input, "r", encoding="utf-8") as input:
        for line in input:
            line = line.strip()
            fields = line.split("\t")
            fields = list(filter(None, fields))

            # If the first field starts with CHROM we're at the beginning of the file
            if fields[0] == "CHROM": # Define fields where we have samples.
                samples = fields[6:]
                genotype_positions = [i for i in range(6, 6+len(samples))]
                outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(samples))
                continue

            # If skipping indels
            if args.skip_indels == True and len(fields[3]) != len(fields[4]):
                continue

            # Skip lines with missing data.
            if "." in "".join(fields[6:]) or "NA" in "".join(fields[6:]):
                continue

            qual = float(fields[5])
            pos = int(fields[1])
            tig = fields[0]
            minpos = min(genotype_positions)
            maxpos = max(genotype_positions) + 1
            genotypes = fields[minpos:maxpos]

            # Convert to list of tuples
            genotypes = [ (int(i.split(",")[0]), int(i.split(",")[1])) \
                        if i != "NA" else (0,0) for i in genotypes ]
            variant_alleles, ref_alleles, cov = [], [], []

            for g in genotypes:
                variant_alleles.append(g[1])
                ref_alleles.append(g[0])
                cov.append(sum(g))

            # Divides two lists element if cov > 0
            variant_frequencies = [ai/bi if bi > 0 else ai for ai,bi in zip(variant_alleles,cov)]

            # If all frequencies are the same, continue. This is the case if
            # all samples are homozygous for either allele.
            if all(x == variant_frequencies[0] for x in variant_frequencies):
                continue

            # Look for homozygous sites. Require also more than 1 read
            # of the minor allele in at least one sample
            if (1.0 in variant_frequencies or 0.0 in variant_frequencies) and \
            max(ref_alleles) >= args.minimum_reads and \
            max(variant_alleles) >= args.minimum_reads:

                # If only filtering out background variants, skip remaining filters
                # and write out this variant.
                if args.background == True:
                    # Convert genotypes back to string
                    genotypes = [ "{},{}".format(i[0],i[1]) for i in genotypes ]
                    outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(genotypes))
                    continue

                # Otherwise, continue filtering
                counter = 0
                homoz_checker, heteroz_checker = False, False
                minor_allele_checker, freq_checker = False, False
                totcov = 0
                tot_ref, tot_alt = 0, 0

                for i in variant_frequencies:

                    # Check total ref/alt frequency at this site,
                    # for heterozygous samples with more than 1 read
                    # containing variant allele
                    if i != 1.0 and i != 0.0:
                        totcov = totcov + cov[counter]

                        # Add number of variant reads to total variant
                        # allele counter
                        tot_alt += variant_alleles[counter]
                        tot_ref += ref_alleles[counter]

                        # Check that at least one heterozygous site has
                        # allele frequency within the boundaries
                        if i > args.frequency_cutoff \
                        and i < 1 - args.frequency_cutoff:
                            freq_checker = True

                        # Check that at least one heterozygous site has
                        # sufficient coverage to consider
                        elif int(cov[counter]) >= args.minimum_coverage \
                        and int(cov[counter]) <= args.maximum_coverage:
                            heteroz_checker = True

                    # Check that at least one homozygous site has
                    # sufficient coverage to consider
                    elif i == 0 or i == 1:
                        if int(cov[counter]) >= args.minimum_coverage \
                        and int(cov[counter]) <= args.maximum_coverage:
                            #print(int(n_o_reads[counter]))
                            homoz_checker = True

                    counter += 1

                # If checkers are true, continue to check minor_allele
                # read frequency
                if homoz_checker == True \
                and heteroz_checker == True \
                and freq_checker == True:
                    # False positives often have several samples with a
                    # low frequency of the variant in the reads.
                    # Real variants, not caused by sporadic read mapping,
                    # should follow a binomial distribution of frequency
                    # in the reads, assuming an equal nuclear ratio in the sample
                    # This is solved by calculating a p-value for the
                    # read frequency across all samples where the
                    # variant is found. Samples with a single read
                    # containing the variant are ignored as this could be
                    # the result of sequencing error or index hopping
                    minor_allele = tot_alt if tot_alt < tot_ref else tot_ref
                    if minor_allele > 1:
                        if args.test_binomial:
                            if test_signif(minor_allele, totcov) == True:
                                # Convert genotypes back to string
                                genotypes = [ "{},{}".format(i[0],i[1]) for i in genotypes ]
                                outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(genotypes))
                        else:
                            genotypes = [ "{},{}".format(i[0],i[1]) for i in genotypes ]
                            outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(genotypes))

        with open(outfile, "w") as out:
            out.write("\n".join(outlines))

if __name__ == "__main__":
    main()
