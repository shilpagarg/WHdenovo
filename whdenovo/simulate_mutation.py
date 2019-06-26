import math
import sys
import random as rn
from collections import defaultdict
from Bio import SeqIO

bases = ["A", "C", "T", "G"]
indel_prob = 0.01
max_indel_length = 4
subs_prob = 0.1


def sample_sites(ref_id, ref_seq, num_sites):
   
    site_list = []
## Remove ambiguous sites
    is_amb = lambda genome, pos : genome[pos] == "N"

    mut_count = 0
## Remove duplicate sites
    x_d = defaultdict(int)

    rlen = len(ref_seq)

    while mut_count < num_sites:
        p = rn.randint(0, rlen - 1)
        if p not in x_d and not is_amb(ref_seq, p):
            x_d[p] = 1
            site_list.append(p)
            mut_count += 1

    return site_list

def get_random_base(orig):
    base_ind = rn.randint(0,3)
    while bases[base_ind] == orig:
        base_ind = rn.randint(0,3)

    return bases[base_ind]


def insert_mutations(ref, sites):
    
    new_seq_list = [i for i in ref]
    pos_nuc_tuples = []
    for s in sorted(sites):
        b = get_random_base(ref[s])
        new_seq_list[s] = b
        pos_nuc_tuples.append( (s, ref[s], b) )
        

    return "".join(new_seq_list), pos_nuc_tuples

def usage():
    print("simple_mutator.py <genomeFA> <divergence> [optional basename]")

if __name__ == "__main__":

    ref_d = {}
    mod_d = {}
    len_d = {}
    num_mutations_d = {}
    site_d = {}
    vcf_d = {}

    ## Read in HPV16 ref
    fa = SeqIO.parse(sys.argv[1], "fasta")
    for rec in fa:
        ref_d[rec.id] = str(rec.seq)
    ## read divergence
    div = float(sys.argv[2]) / float(100)

    ## optional basename
    namebase = ""
    if len(sys.argv) > 3:
        namebase = str(sys.argv[3])

    ## Calculate number of mutations to insert
    len_d = {i : len(ref_d[i]) for i in ref_d}
    num_mutations_d = {i : math.ceil( float(div) * len_d[i] ) for i in len_d}

    ## get sites to modify
    site_d = {i : sample_sites(i, ref_d[i], num_mutations_d[i]) for i in ref_d}
    
    ## modify sites
    for i in ref_d:
        nseq = insert_mutations(ref_d[i], site_d[i])
        mod_d[i] = nseq[0]
        vcf_d[i] = nseq[1]

    if namebase is "":
        namebase = "mod"
    ## write new genome
    if namebase is not "":
        with open(namebase + ".fasta", "w") as ofi:
            for i in mod_d:
                ofi.write(">" + namebase + "_" + i + "\n")
                ofi.write(mod_d[i] + "\n")
        ## write pseudo-VCF of mutations introduced
        with open(namebase + ".pvcf", "w") as vcfi:
            for i in mod_d:
                vcfi.write("##newname=" + namebase + "_" + i + "\n")
                vcfi.write("##oldname=" + i + "\n")
                vcfi.write("##percentDivergence=" + str(div) + "\n")
                vcfi.write("##numMutations=" + str(num_mutations_d[i]) + "\n")

                for s in vcf_d[i]:
                    vcfi.write( "\t".join([str(t) for t in s]) + "\n")
    


