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

__version__ = "0.41"

parser = argparse.ArgumentParser(description='Filter variants to find potential \
                                            mutations. Input: variant table from \
                                            GATK VariantsToTable. Should be of the \
                                            following format: \
                                            "gatk VariantsToTable -V $VCF -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -GF AD -O $OUTPUT" \
                                            data for each sample start at the 6th. \
                                            The AD field in the table should be comma separated \
                                            and contain the number of reads with reference \
                                            allele and number of reads with alternate allele, \
                                            e.g. "23,21".')
parser.add_argument("input", \
                    help="Input table.", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Output prefix.", \
                    type = str)
parser.add_argument("-c","--coverage_file", \
                    help="Tab separated file containing each sample in the \
                        input table in the first column, and the median coverage \
                        for that sample in the second.", \
                    type = str)
parser.add_argument("-n","--minimum_coverage", \
                    help="Minimum coverage to keep a variant [20].", \
                    type = int, \
                    default = 20)
parser.add_argument("-x","--maximum_coverage", \
                    help="Maximum coverage to keep a variant [60].", \
                    type = int, \
                    default = 60)
parser.add_argument("-a","--coverage_factor", \
                    help="Allow up to this much deviation from the median \
                        coverage defined in the file given by -c. Higher => less \
                        conservative [0.2].", \
                    type = float, \
                    default = 0.2)
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

class Genotype:
    ''' Describe the genotype of each sample on the reference and variant read coverage
    '''
    def __init__(self, sample, ref_reads, var_reads, cov, var_freq):
        self.sample = sample
        self.ref_reads = ref_reads
        self.var_reads = var_reads
        self.cov = cov
        self.var_freq = var_freq

def getOut():
    '''
    To create a prefix for output files.
    '''
    if args.output:
        outprefix = args.output
    else:
        outprefix = "potential_mutations"
    return outprefix

def readCov():
    ''' Read the median coverage for each sample in the dataset from a separate
        file.
    '''
    median_covs = {}
    with open(args.coverage_file, "r", encoding="utf-8") as cov:
        for line in cov:
            line = line.strip()
            sample = line.split("\t")[0]
            coverage = int(line.split("\t")[1])
            median_covs[sample] = coverage
    return median_covs

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

            # Collect data into a list of Genotypes
            reads_per_sample = fields[min(genotype_positions):max(genotype_positions) + 1] # It's still a list of strings
            reads_per_sample = [ (int(i.split(",")[0]), int(i.split(",")[1])) \
                                for i in reads_per_sample ] # Convert to list of tuples
            genotypes = []
            variant_frequencies = []
            for sample,gt in zip(samples,reads_per_sample):
                ref_reads,variant_reads = gt[0],gt[1]
                cov = variant_reads+ref_reads
                freq = variant_reads/cov if cov > 0 else 0
                genotypes.append(Genotype(sample, ref_reads,variant_reads, cov, freq))
                variant_frequencies.append(freq)

            # If all frequencies are the same, continue. This is the case if
            # all samples are homozygous for either allele.
            if all(x == variant_frequencies[0] for x in variant_frequencies):
                continue

            # Look for homozygous sites. Require also more than 1 read
            # of the minor allele in at least one sample
            if (1.0 in variant_frequencies or 0.0 in variant_frequencies) and \
            max([x.ref_reads for x in genotypes]) >= args.minimum_reads and \
            max([x.var_reads for x in genotypes]) >= args.minimum_reads:

                # If only filtering out background variants, skip remaining filters
                # and write out this variant.
                if args.background == True:
                    outlines.append(line)
                    continue

                # Otherwise, continue filtering
                homoz_checker, heteroz_checker = False, False
                indv_cov_checker, freq_checker = True, False
                tot_ref, tot_alt, totcov = 0,0,0

                for gt in genotypes:

                    # Check total ref/alt frequency at this site,
                    # for heterozygous samples with more than 1 read
                    # containing variant allele
                    if gt.var_freq != 1.0 and gt.var_freq != 0.0 \
                    and gt.var_reads > 1:
                        totcov = totcov + gt.cov

                        # Add number of variant reads to total variant
                        # allele counter
                        tot_alt += gt.var_reads
                        tot_ref += gt.ref_reads

                        # Check that at least one heterozygous site has
                        # allele frequency within the boundaries
                        if gt.var_freq > args.frequency_cutoff \
                        and gt.var_freq < 1 - args.frequency_cutoff:
                            freq_checker = True

                        # Check that at least one heterozygous site has
                        # sufficient coverage to consider
                        elif int(gt.cov) >= args.minimum_coverage \
                        and int(gt.cov) <= args.maximum_coverage:
                            heteroz_checker = True

                        # If user gives a file with median coverage values (-c),
                        # we check that the coverage for all heterozygous
                        # samples is within the defined boundaries (-a)
                        if args.coverage_file:
                            median_covs = readCov()
                            median = median_covs[gt.sample]
                            # Check if coverage is outside the boundaries
                            # if so we skip this site by turning the
                            # checker to false
                            if gt.cov > median + median * args.coverage_factor \
                            or gt.cov < median - median * args.coverage_factor:
                                indv_cov_checker = False
                                break

                    # Check that at least one homozygous site has
                    # sufficient coverage to consider
                    elif gt.var_freq == 0 or gt.var_freq == 1:
                        if int(gt.cov) >= args.minimum_coverage \
                        and int(gt.cov) <= args.maximum_coverage:
                            #print(int(n_o_reads[counter]))
                            homoz_checker = True

                # If checkers are true, continue to check minor_allele
                # read frequency
                if heteroz_checker == True \
                and freq_checker == True \
                and indv_cov_checker == True \
                and homoz_checker == True:

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
                    tot_minor_allele = tot_alt if tot_alt < tot_ref else tot_ref

                    if tot_minor_allele > 1:
                        if args.test_binomial:
                            if test_signif(tot_minor_allele, totcov) == True:
                                outlines.append(line)
                        else:
                            outlines.append(line)

        with open(outfile, "w") as out:
            out.write("\n".join(outlines))

if __name__ == "__main__":
    main()
