#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
dnds.py
Author: Markus Hiltunen Thorén
E-mail: markus.hiltunen@su.se

Calculate the dN/dS ratio for a set of genes that carry genetic variation.
Assumes closely related taxa, such as individuals from a population, and
does no correction for sequence evolution.

In short, the script reads a vcf file and collects single nucleotide variants,
compares these to coding regions in a gff file and predicts amino acid
changes. Also requires the genome in a fasta file to get initial gene
sequences from. The genes in the gff file have to have entries for both gene
and CDS (e.g. maker output).

Originally called analyze_coding_SNPs.py when it was used for the first
Marasmius mutation paper (Hiltunen et al. 2019 Curr. Biol.)

Copyright (c) 2023, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

__version__ = "1.0"

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Analyzes SNPs in coding regions \
                                to predict their amino acid changes, then uses this \
                                information to calculate dN/dS statistics.")
parser.add_argument("reference", \
                    help="Reference genome fasta file. Required.", \
                    type = str)
parser.add_argument("annotation", \
                    help="Genome annotation gff file. Has to be coordinate sorted! Required.", \
                    type = str)
parser.add_argument("variants", \
                    help="Variant call file in tab delimited format. \
                    Can be a vcf file, but all variants are assumed \
                    to be present in the query sample, and vcf genotype \
                    fields are NOT taken into account. Has to have the \
                    following fields: \n\
                    <scaffold> <coordinate> <anything> <ref_allele> <alt_allele> \
                    Required", \
                    type = str)
parser.add_argument("-i", "--include_all", \
                    help="Add to include fragmented gene annotations, \
                    i.e. where the spliced CDS length is not divisible \
                    by three, or where start/stop codons deviate", \
                    action="store_true")
parser.add_argument("-s", "--write_seq", \
                    help="Add to write fasta files for the genes that were \
                    included in the analysis", \
                    action="store_true")
parser.add_argument("-v","--version", \
                    help="Print version and quit.", \
                    action = "version", \
                    version = "dnds.py v.{}".format(__version__))
args = parser.parse_args()

def readGeneticCode():
    '''
    Creates a dict with the genetic code.
    '''
    global genetic_code
    genetic_code = {"TTT":"Phe",
                    "TCT":"Ser",
                    "TAT":"Tyr",
                    "TGT":"Cys",
                    "TTC":"Phe",
                    "TCC":"Ser",
                    "TAC":"Tyr",
                    "TGC":"Cys",
                    "TTA":"Leu",
                    "TCA":"Ser",
                    "TAA":"Ter",
                    "TGA":"Ter",
                    "TTG":"Leu",
                    "TCG":"Ser",
                    "TAG":"Ter",
                    "TGG":"Trp",
                    "CTT":"Leu",
                    "CCT":"Pro",
                    "CAT":"His",
                    "CGT":"Arg",
                    "CTC":"Leu",
                    "CCC":"Pro",
                    "CAC":"His",
                    "CGC":"Arg",
                    "CTA":"Leu",
                    "CCA":"Pro",
                    "CAA":"Gln",
                    "CGA":"Arg",
                    "CTG":"Leu",
                    "CCG":"Pro",
                    "CAG":"Gln",
                    "CGG":"Arg",
                    "ATT":"Ile",
                    "ACT":"Thr",
                    "AAT":"Asn",
                    "AGT":"Ser",
                    "ATC":"Ile",
                    "ACC":"Thr",
                    "AAC":"Asn",
                    "AGC":"Ser",
                    "ATA":"Ile",
                    "ACA":"Thr",
                    "AAA":"Lys",
                    "AGA":"Arg",
                    "ATG":"Met",
                    "ACG":"Thr",
                    "AAG":"Lys",
                    "AGG":"Arg",
                    "GTT":"Val",
                    "GCT":"Ala",
                    "GAT":"Asp",
                    "GGT":"Gly",
                    "GTC":"Val",
                    "GCC":"Ala",
                    "GAC":"Asp",
                    "GGC":"Gly",
                    "GTA":"Val",
                    "GCA":"Ala",
                    "GAA":"Glu",
                    "GGA":"Gly",
                    "GTG":"Val",
                    "GCG":"Ala",
                    "GAG":"Glu",
                    "GGG":"Gly"}


def calcdnds(syn_sites,nonsyn_sites,syn_subst,nonsyn_subst):
    '''
    Calculates the dN/dS ratio from given files.
    '''
    pS, pN = syn_subst/syn_sites, nonsyn_subst/nonsyn_sites
    dS, dN = -3/4*np.log(1-(4*pS)/3), -3/4*np.log(1-(4*pN)/3)
    dnds = dN/dS
    return pS, pN, dS, dN, dnds

def collectVariants():
    '''
    Reads the vcf to collect the single nucleotide variants in a dict.
    '''
    with open(args.variants, "r") as vcf:
        variants = {}
        n_indels, n_snps = 0,0
        for line in vcf:
            line = line.strip()
            if not line.startswith("#"):
                fields = line.split("\t")
                tig, coord = fields[0], fields[1]
                ref_allele, alt_allele = fields[3], fields[4]
                if len(ref_allele) == 1 and len(alt_allele) == 1: # Ignore indels
                    variants[(tig, coord)] = (ref_allele, alt_allele)
                    n_snps += 1
                else:
                    n_indels += 1

        return variants, n_indels, n_snps

def readFasta():
    '''
    Reads the input fasta file into a dict.
    '''
    fa = {}
    with open(args.reference, "r") as fasta:
        sequence = None
        for line in fasta:
            line = line.strip()

            if line.startswith(">") and sequence == None:
                header = line[1:]
                sequence = []

            elif line.startswith(">") and sequence != None:
                # If new fasta entry, add old one to dict,
                # pick new header and reset sequence
                fa[header] = "".join(sequence)
                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        # Last passthrough won't have any new entries,
        # just add the remaining sequence
        fa[header] = "".join(sequence)

    return fa

def writeFasta(fasta, fileprefix):
    with open(fileprefix+".fasta", "w") as out:
        for k,v in fasta.items():
            header = ">{} {}_{}_{}_{}".format(k[4],k[0],k[1],k[2],k[3])
            out.write(header+"\n")
            out.write(v+"\n")


def reverse_complement(nuclstring):
    '''
    Reverse complements input nucleotide sequence.
    '''
    rev_comped = ""
    for l in reversed(nuclstring):
        if l == "A" or l == "a":
            rev_comped += "T"
        elif l == "T" or l == "t":
            rev_comped += "A"
        elif l == "C" or l == "c":
            rev_comped += "G"
        elif l == "G" or l == "g":
            rev_comped += "C"
        elif l == "N" or l == "n":
            rev_comped += "N"
    return rev_comped

def collectTranscripts(genes, fasta, ref_bool, all_bool):
    '''
    Collects spliced gene sequences.
    '''
    spliced_genes = {}

    for k,exons in genes.items():
        exon_sequence_list = []
        tig,coord1,coord2,strand,ID = k[0],k[1],k[2],k[3],k[4]# Unpack k

        for ex in exons:
            seq = fasta[ k[0] ][ int(ex[0])-1:int(ex[1]) ] # Spliced coordinates assume sorted gff format.
            exon_sequence_list.append(seq.upper())

        # After collecting all the exons, splice the resulting CDS
        spliced_seq = "".join(exon_sequence_list)
        if strand == "-": # Reverse complement if we're on the minus strand
            spliced_seq = reverse_complement(spliced_seq)

        if not ref_bool:
            # ref_bool is to say if we are splicing reference or variant cds
            # if false, we are splicing var genes, and have already gone through this
            # function once with the reference genes, and filtered the gene set.
            # Thus we can simply collect all remaining genes.
            spliced_genes[k] = spliced_seq
            continue

        else:
            # If we are dealing with ref genes, we need to do some quality control
            # to make sure we are looking at good annotations.
            # Collect the start codon, we need this for quality control
            startcodon = spliced_seq[:3]

            # We need to check that each gene:
            # 1. Is divisible by three
            # 2. Has ATG as start codon
            # 3. Has a valid stop codon

            if not len(spliced_seq)%3 == 0:
                # If spliced gene is an uneven number of bases, skip it,
                # unless the user gave the --include_all parameter.
                # If so, we check if we can keep it.
                if not all_bool:
                    print("WARNING Gene {}: uneven number of bases - skipping (add --include_all to try to keep it).".format(str(ID)))
                    continue
                elif startcodon == "ATG":
                    # If correct start codon, we simply give a warning and remove the last incomplete codon,
                    # forcing the gene to be divisible, and then keep it.
                    remaining = len(spliced_seq)%3
                    spliced_seq = spliced_seq[:-remaining] # Remove the last remaining bases
                    spliced_genes[k] = spliced_seq
                    print("WARNING Gene {}: uneven number of bases - reading frame determined by start codon.".format(str(ID)))
                else:
                    print("WARNING Gene {}: uneven number of bases - unable to determine reading frame, skipping gene".format(str(ID)))
                    continue
            else:
                # Then we check for annotations with faulty start or stop codons.
                stopcodon = spliced_seq[-3:] # Now that we know that gene spliced sequence is divisible by three, 
                                             # we can collect the stop codon.
                if startcodon == "ATG" and stopcodon in ("TAG","TAA","TGA"):
                    spliced_genes[k] = spliced_seq
                elif all_bool: # If user wants to keep as many genes as possible, warn about problematic genes.
                    if startcodon != "ATG":
                        print("WARNING Gene {}: non-canonical start codon ({}). Is the GFF sorted?".format(str(ID), startcodon))
                    if stopcodon not in ("TAG","TAA","TGA"):
                        print("WARNING Gene {}: non-canonical stop codon ({}). Is the GFF sorted?".format(str(ID), stopcodon))
                    spliced_genes[k] = spliced_seq

    return spliced_genes

def swapbases(variants, genome):
    '''
    Swaps bases at variant positions in the genome fasta.
    '''
    variant_genome = {}
    for k,v in variants.items():
        tig = k[0]
        coord = int(k[1])
        ref_base = v[0]
        alt_base = v[1]

        if tig not in variant_genome.keys():
            variant_genome[tig] = "".join([ genome[tig][:coord-1], \
                                            alt_base, \
                                            genome[tig][coord:] ])
        else:
            s = variant_genome[tig]
            variant_genome[tig] = "".join([ s[:coord-1], \
                                            alt_base, \
                                            s[coord:] ])

    return variant_genome

def findVariantGenes(variants, genes):
    '''
    Filters the genes dict to only include genes with variants in them.
    Also returns the variants in coding regions.
    '''
    filtered_genes = {}
    coding_variants = []

    for key in variants.keys():
        tig = key[0]
        coord = int(key[1])

        for k,v in genes.items():
            # Check if tig is the same and coord is within the gene boundaries
            # Both VCF and GFF are 1-based
            if k[0] == tig:
                if int(k[1]) <= coord and int(k[2]) >= coord:

                    # If so, go through exons to check if variant
                    # is in an exon or intron
                    for exon in v:
                        c1, c2 = int(exon[0]), int(exon[1])

                        if c1 <= coord and c2 >= coord:
                            filtered_genes[k] = v
                            coding_variants.append(key)

    return filtered_genes, coding_variants

def codonsplit(dnaseq):
    '''
    Splits a DNA sequence string into a list of codons.
    '''
    return [dnaseq[i:i+3] for i in range(0, len(dnaseq), 3)]

def translate(geneseq):
    '''
    Translates a gene sequence into a protein sequence.
    '''
    protseq = []
    for codon in codonsplit(geneseq):
        protseq.append( genetic_code[codon] )
    return "".join(protseq)

def calcPossibleChanges(geneseq):
    '''
    Calculates the number of possible synonymous and nonsynonymous changes.
    in a gene sequence.
    '''
    totsyn = 0
    totnonsyn = 0
    for codon in codonsplit(geneseq):
        nonsyn = 0
        for idx,base in enumerate(codon):
            bases = ["A", "T", "C", "G"]
            bases.remove(base)
            for b in bases:
                mutated_codon = "".join( [codon[:idx], b, codon[idx+1:]] )
                if translate(mutated_codon) != translate(codon):
                    nonsyn += 1/3
        syn = 3 - nonsyn
        totsyn += syn
        totnonsyn += nonsyn

    return totsyn, totnonsyn


def collectCDS():
    '''
    Collects coding region coordinates.
    '''
    with open(args.annotation, "r") as gff:
        genes = {}
        coord1,coord2,strand = 0,0,"+"
        for line in gff:
            line = line.strip()
            # Some gff files have the genome sequence at the end of the 
            # file, in fasta format. Break if we get to that point.
            if line.startswith(">"): 
                break
            if not line.startswith("#"):
                fields = line.split("\t")
                tig = fields[0]
                # Here we collect the coding region coordinates into a dict
                if fields[2] == "gene":
                    coord1,coord2,strand = fields[3],fields[4],fields[6]
                    ID = [ t[3:] for t in fields[8].split(";") if t[:3] == "ID=" ][0] # Find the ID tag
                    genes[(tig,coord1,coord2,strand,ID)] = []
                elif fields[2] == "CDS":
                    genes[(tig,coord1,coord2,strand,ID)].append((fields[3], \
                                                              fields[4]))

        return genes

def calcObserved(cds,variant_cds):
    ''' Calculates the number of synonymous and nonsynonymous changes
        in two gene sequences.
    '''
    syn,nonsyn = 0,0
    per_gene = {}

    for key,ref_seq in cds.items():
        alt_seq = variant_cds[key]
        ref_codons = codonsplit(ref_seq)
        alt_codons = codonsplit(alt_seq)
        per_gene[key] = np.zeros(9)

        for idx, codon in enumerate(ref_codons):
            #import ipdb;ipdb.set_trace()
            if codon != alt_codons[idx]: # If codon is the same, continue
                # If codons are not the same, and the amino acid is also not the same,
                # we have a nonsyn change
                if genetic_code[codon] != genetic_code[ alt_codons[idx] ]:
                    nonsyn += 1
                    per_gene[key][3] += 1
                # Else we have a syn change
                else:
                    syn += 1
                    per_gene[key][2] += 1
    return syn, nonsyn, per_gene

def calcNonCDSGenomeSize(genome, n_coding_bases):
    '''
    Calculates the noncoding size of the genome.
    '''
    genomesize = 0
    for value in genome.values():
        genomesize += len(value)
    return genomesize - n_coding_bases

def drawDnDs(variant_genes):
    '''Draw dN vs dS per gene as a scatter plot.
    '''
    dNs, dSs = [], []

    for v in variant_genes.values():
        pS, pN, dS, dN, dnds = calcdnds(v[0],v[1],v[2],v[3])
        dNs.append(dN)
        dSs.append(dS)

    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.scatter(dNs, dSs)
    ax.set_xlabel("dN")
    ax.set_ylabel("dS")
    plt.show()


def main():
    variants, n_indels, n_snps = collectVariants()
    genome = readFasta()
    genes = collectCDS()
    readGeneticCode()
    filtered_genes, coding_variants = findVariantGenes(variants,genes) # Collect the genes with SNPs
    n_coding = len(coding_variants) # Number of SNPs from the vcf file that were located in CDS
    n_noncoding = len(variants) - n_coding # Number of SNPs not in CDS
    # Collect the spliced coding sequences for all genes with one or more SNPs
    # Here, genes with non-canonical start/stop codon annotations are removed
    cds = collectTranscripts(filtered_genes, genome, True, args.include_all)
    # Keep only those genes in the filtered_genes
    hq_genes = { key: filtered_genes[key] for key in cds.keys() }
    # Count the number of remaining variants. I use the same function, which is probably inefficient.
    x, filtered_variants = findVariantGenes(variants,hq_genes)
    n_filtered_variants = len(filtered_variants)

    # Create the sequence for the variant genome
    variant_genome = swapbases(variants, genome) 
    variant_cds = collectTranscripts(hq_genes,variant_genome, False, False)
    n_variant_genes = len(filtered_genes)
    n_variant_hq_genes = len(hq_genes)

    if args.write_seq:
        ref_fname = args.reference.split("/")[-1].rstrip("fasta") + "CDS"
        var_fname = args.variants.split("/")[-1].rstrip("vcf") + "CDS"
        writeFasta(cds, ref_fname)
        writeFasta(variant_cds, var_fname)
        print("\nSpliced reference CDS written to file: {}".format(ref_fname+".fasta"))
        print("Spliced variant CDS written to file: {}\n".format(var_fname+".fasta"))

    # variant_genes_table collects numbers of syn/nonsyn changes to write to a file
    obs_syn, obs_nonsyn, variant_genes_table = calcObserved(cds,variant_cds)

    tot_synonymous = 0
    tot_nonsynonymous = 0

    # Calculate number of possible changes for each gene
    for gene, seq in cds.items():
        syn, nonsyn = calcPossibleChanges(seq)
        tot_synonymous += syn
        tot_nonsynonymous += nonsyn
        variant_genes_table[gene][0] = syn
        variant_genes_table[gene][1] = nonsyn
        pS, pN, dS, dN, dnds = calcdnds(syn,nonsyn,variant_genes_table[gene][2],variant_genes_table[gene][3])
        variant_genes_table[gene][4] = pS
        variant_genes_table[gene][5] = pN
        variant_genes_table[gene][6] = dS
        variant_genes_table[gene][7] = dN
        variant_genes_table[gene][8] = dnds

    # Full dataset concatenated
    pS, pN, dS, dN, dnds = calcdnds(tot_synonymous,tot_nonsynonymous,obs_syn,obs_nonsyn)

    # Write statistics per gene to file
    outfilename = args.variants.split("/")[-1].rstrip("vcf") + "dnds.tsv"
    with open(outfilename, "w") as variant_genes_out:
        variant_genes_out.write("Gene\tScaffold\tStart_coord\tEnd_coord\tStrand\tPoss_syn\tPoss_nonsyn\tObs_syn\tObs_nonsyn\tpS\tpN\tdS\tdN\tdN/dS\n")

        for k,v in variant_genes_table.items():
            tig,coord1,coord2,strand,ID = k[0],k[1],k[2],k[3],k[4]
            outstr = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ID,tig,coord1,coord2,strand,v[0], v[1], int(v[2]), int(v[3]), v[4], v[5],v[6],v[7],v[8])
            variant_genes_out.write(outstr)

    print("\nStatistics per gene written to file: {}".format(outfilename))
    print("\n#### Full dataset statistics ####")
    print("Number of variants: {}".format(n_indels+n_snps))
    print(" of which are indels: {}".format(n_indels))
    print(" of which are SNPs: {}".format(n_snps))
    print("  of which are in CDS: {}".format(n_coding))
    print("   of which are in genes passing QC: {}".format(n_filtered_variants))
    print("  of which are outside CDS: {}".format(str(n_noncoding)))
    print("Number of genes with at least one SNP: {}".format(n_variant_genes))
    print(" of which passed QC: {}".format(n_variant_hq_genes))
    print("\n#### dN/dS statistics ####")
    print("Calculated based on genes that passed QC (n={}) and SNPs within the CDS of these genes (n={})\n".format(n_variant_hq_genes,n_filtered_variants))
    print("Synonymous sites: {}".format(tot_synonymous))
    print("Nonsynonymous sites: {}".format(tot_nonsynonymous))
    print("Synonymous substitutions: {}".format(obs_syn))
    print("Nonsynonymous substitutions: {}".format(obs_nonsyn))
    print("dN: " + str(dN))
    print("dS: " + str(dS))
    print("dN/dS: " + str(dnds))

    #import matplotlib.pyplot as plt
    #drawDnDs(variant_genes_table)

if __name__ == "__main__":
    main()
