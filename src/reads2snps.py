# Script for associating reads with variants reported by discoSNP++
# The output is written to stdout
# Usage:
# python3 reads2snps.py -V <variant fasta file> -F <read file>

import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--variants', metavar = 'FILE', type = str, required = True, 
                     help = 'Fasta file with variants produce by DiscoSNP++')
parser.add_argument('-F', '--fasta', metavar = 'FILE', type = str, required = True, 
                     help = 'File with reads')

args = parser.parse_args()

varFile = args.variants
fastaFile = args.fasta

variants = dict()    # Dictionary of mapping k-mers to variant names
varlen = -1

# Read the variant fasta file and put the variants into a dictionary
def read_variants(varFile):
    fasta_seqs = SeqIO.parse(open(varFile), 'fasta')
    for fasta in fasta_seqs:
        name,seq = fasta.id, str(fasta.seq)
        snp_name = name.split("|")[0].split("_")[3]
        allele = name.split("|")[0].split("_")[1]
        snpornot = name.split("|")[0].split("_")[0]
        if not snpornot == "SNP":
            continue
        if allele == 'higher':
            allele = 'H'
        elif allele == 'lower':
            allele = 'L'
        else:
            sys.stderr.write("Unknown allele: " + allele)
            exit(1)
        rseq = str(Seq.Seq(seq, generic_dna).reverse_complement())
        variants[seq] = (snp_name, 'f',  allele)
        variants[rseq] = (snp_name, 'r',  allele)
        varlen = len(seq)
    return varlen

# Read the read fasta file and find the snps/alleles occurring in the read
def find_snps_in_reads(fastaFile):
    fasta_seqs = SeqIO.parse(open(fastaFile), 'fasta')
    for fasta in fasta_seqs:
        name,seq = fasta.id, str(fasta.seq)
        snps = []
        for i in range(0, len(seq)-varlen):
            kmer = seq[i:i+varlen]
            if kmer in variants:
                snps.append(variants[kmer])
        print(name)
        print(snps)
        
varlen=read_variants(varFile)

find_snps_in_reads(fastaFile)
