"""
create association between reads and bubbles.
"""
import pyfaidx

from xopen import xopen
import stream
import logging
from . import vg_pb2
from collections import Counter
import sys
from collections import defaultdict
from .core import ReadSet, Read
from functools import reduce
import operator as op
from itertools import groupby
from .graph import ComponentFinder

import random
import collections
from collections import OrderedDict, namedtuple


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
from .core import ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .variants import ReadSetReader, ReadSetError
from heapq import heappush, heappop
from itertools import count


__author__ = "Shilpa Garg"

logger = logging.getLogger(__name__)


class CoverageMonitor:
        '''TODO: This is a most simple, naive implementation. Could do this smarter.'''
        def __init__(self, length):
                self.coverage = [0] * length

        def max_coverage_in_range(self, begin, end):
                return max(self.coverage[begin:end])

        def add_read(self, begin, end):
                for i in range(begin, end):
                        self.coverage[i] += 1

def setup_pedigree(ped_path, numeric_sample_ids, samples):
    """
    Read in PED file to set up list of relationships.

    Return a pair (trios, pedigree_samples), where trios is a list of Trio
    objects and pedigree_samples is the set of all samples that are mentioned
    in the PED file (as individual, mother or father).

    ped_path -- path to PED file
    samples -- samples to be phased
    """
    trios = []
    pedigree_samples = set()
    for trio in PedReader(ped_path, numeric_sample_ids):
        if (trio.child is None or trio.mother is None or
                    trio.father is None):
            logger.warning('Relationship %s/%s/%s ignored '
                           'because at least one of the individuals is unknown',
                           trio.child, trio.mother, trio.father)
        else:
            # if at least one individual is not is samples, skip trio
            if( (trio.mother in samples) and (trio.father in samples) and (trio.child in samples) ):
                trios.append(trio)
                pedigree_samples.add(trio.child)
                pedigree_samples.add(trio.father)
                pedigree_samples.add(trio.mother)
            else:
                # happens in case --ped and --samples are used
                logger.warning('Relationship %s/%s/%s ignored because at least one of the '
                        'individuals was not given by --samples.', 
                        trio.child, trio.mother, trio.father)

    return trios, pedigree_samples


def find_components(phased_positions, reads, master_block=None, heterozygous_positions=None):
    """
    Return a dict that maps each variant position to the component it is in.
    Variants are considered to be in the same component if a read exists that
    covers both. A component is identified by the position of its leftmost
    variant.
    master_block -- List of positions in a "master block", i.e. all blocks containing
                    any of these positions are merged into one block.
    heterozygous_positions -- A dictionary mapping numeric sample ids to sets of
                              positions. Component building is then restricted to variants
                              at these positions. If none, all variants are used.
    """
    logger.debug('Finding connected components ...')
    assert phased_positions == sorted(phased_positions)

    # Find connected components.
    # A component is identified by the position of its leftmost variant.
    component_finder = ComponentFinder(phased_positions)
    phased_positions = set(phased_positions)
    for read in reads:
        if heterozygous_positions is None:
            positions = [ variant.position for variant in read if variant.position in phased_positions ]
        else:
            positions = [ variant.position for variant in read \
                if (variant.position in phased_positions) and (variant.position in heterozygous_positions[read.sample_id])
            ]
        for position in positions[1:]:
            component_finder.merge(positions[0], position)
    if not master_block is None:
        for position in master_block[1:]:
            component_finder.merge(master_block[0], position)
    components = { position : component_finder.find(position) for position in phased_positions }
    return components


def haplotag(pred_superreads, read_set, components, ap, locus_file, iteration, f):
    phases = []
    
    for s1,s2 in zip(*pred_superreads):
        VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'position', 'phase1', 'phase2'])
        extract_phase = VariantCallPhase(components[s1.position], s1.position, s1.allele, s2.allele)
        phases.append(extract_phase) #TODO: check it
        
    variantpos_to_allele1 = {
        v.position:int(v.allele) for v, v2 in zip(*pred_superreads) if v.allele!=-2
    }
    variantpos_to_allele2 = {
        v2.position:int(v2.allele) for v, v2 in zip(*pred_superreads) if v2.allele!=-2
    }
    
    variantpos_to_phaseset1 = {
        v.position:components[v.position] for v, v2 in zip(*pred_superreads) if v.allele!=-2
    }
    variantpos_to_phaseset2 = {
        v2.position:components[v2.position] for v, v2 in zip(*pred_superreads) if v2.allele!=-2
    }
    print('i am here')
    print(components)
    read_to_haplotype = {}
    #read_set = read_reads(readset_reader, chromosome_name, variants, sample, fasta)
    for read in read_set:
        # mapping: phaseset --> phred scaled difference between costs of assigning reads to haplotype 0 or 1
        haplotype_costs = defaultdict(int)
        haplotype_costs[0] = 0
        haplotype_costs[1] = 0
        for v in read:
            if v.position not in variantpos_to_allele1 or v.position not in variantpos_to_allele2:
                continue
            phaseset1 = variantpos_to_allele1[v.position]
            phaseset2 = variantpos_to_allele2[v.position]
    
            if v.allele != phaseset1:
                haplotype_costs[0] += 10
            #else:
                #haplotype_costs[0] -= 10
            if v.allele != phaseset2:
                haplotype_costs[1] += 10
            #else:
                #haplotype_costs[1] -= 10
            
        l =[]
        for k,v in haplotype_costs.items():
            l.append(k)
            l.append(v)
        #l = list(haplotype_costs.items())
        if l[1] < l[3]:
            read_to_haplotype[read.name] = (0, l[1], l[0])
        else:
            read_to_haplotype[read.name] = (1, l[3], l[2])
        #l.sort(key=lambda t:-abs(t[1]))
        #phaseset, quality = l[0]
        #if quality != 0:
            #haplotype = 0 if quality > 0 else 1
            #read_to_haplotype[read.name] = (haplotype, abs(quality), phaseset)
    #print(read_to_haplotype)
    #f = open(str(locus_file)+".allreads", 'a')

    for read in read_set:
        #print(read.name)
        for variant in read:
            if variant.allele!=-2:
                if variant.position in ap:
                    phaseset = components[variant.position] + 1
                    break
        print(read.name, str(locus_file)+ str(phaseset), read_to_haplotype[read.name][2], len(read), file=f)

def find_largest_component(components):
    """
    Determine the largest component and return a sorted list of positions
    contained in it.
    components -- dictionary mapping positin to block_id as returned by find_components.
    """
    blocks = defaultdict(list)
    for position, block_id in components.items():
        blocks[block_id].append(position)
    largest = []
    for block in blocks.values():
        if len(block) > len(largest):
            largest = block
    largest.sort()
    return largest
      

"""
output the possible allele-pairs for a bubble.
"""
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom
  

"""
build node-sequence list for vg graph
"""
def vg_graph_reader(vg_file):
    node_seq_list= defaultdict()
    edge_connections = defaultdict(list)
    with stream.open(str(vg_file), "rb") as istream:
        for data in istream:
            l = vg_pb2.Graph()
            l.ParseFromString(data)
            for i in range(len(l.node)):
                index = l.node[i].id
                seq = l.node[i].sequence
                node_seq_list[index]=seq
            for j in range(len(l.edge)):
                from_edge = getattr(l.edge[j], "from")
                edge_connections[from_edge].append(l.edge[j].to)
    return node_seq_list, edge_connections


"""
Output phased SnarlTraversal using superreads_list. 
Then take input vg graph and phased SnarlTraversal to reconstruct underlying two sequences.
"""
"""
Input: Phase variants from Locus file and aligned reads from GAM file.

It creates an association between het variants and read alignments. 

Output: The optimal partitioning is written to standard output.
"""
def vg_reader(locus_file, gam_file, sample):
    """
    input: sorted locus and sorted GAM file output from vg.
    output: sorted readset for core DP.
    assumptions: 
    1. locus file consists of linear ordering of simple bubbles only and hence sorted. Each locus file does not contain start and end vertex.
    2. paths in the locus should be covered by atleast one pacbio read.
    2. GAM file is sorted and restricted to locus file.
    3. files consists of all DAG connected components.
    4. add variant only when it identifies the branch uniquely.
    """

    locus_count=0
    prev_startsnarl = 0
    prev_endsnarl = 0
    locus_branch_mapping=OrderedDict()
    locus_count=0
    prev_startsnarl = 0
    prev_startsnarl_orientation = -1
    prev_endsnarl = 0
    prev_endsnarl_orientation = -1

    with stream.open(str(locus_file), "rb") as istream:
        for data in istream:
            l = vg_pb2.SnarlTraversal()
            l.ParseFromString(data)
            #TODO: make ordered doctionary locus_branch_mapping
            # handle forward and backward case of nodes
            current_startsnarl = l.snarl.start.node_id
            current_startsnarl_orientation = l.snarl.start.backward
            current_endsnarl = l.snarl.end.node_id
            current_endsnarl_orientation = l.snarl.end.backward
            path_in_bubble =[]

            if len(l.visits) ==0:
                #TODO: for now, assumed, all nodes in path are either forward or backward
                if l.snarl.start.backward == True:
                    path_in_bubble.append(tuple ((l.snarl.end.node_id,l.snarl.start.node_id)))
                else:
                    path_in_bubble.append(tuple ((l.snarl.start.node_id,l.snarl.end.node_id)))
            else:
                #TODO: for now, assumed, all nodes in path are either forward or backward
                if l.snarl.start.backward == True:
                    path_in_bubble.append(tuple ((l.snarl.end.node_id, l.visits[-1].node_id)))
                    for i in range(0,len(l.visits)-1):
                        path_in_bubble.append(tuple((l.visits[i+1].node_id, l.visits[i].node_id)))
                    path_in_bubble.append(tuple ((l.visits[0].node_id,l.snarl.start.node_id)))
                else:
                    path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
                    for i in range(0,len(l.visits)-1):
                        path_in_bubble.append(tuple((l.visits[i].node_id, l.visits[i+1].node_id)))
                    path_in_bubble.append(tuple ((l.visits[-1].node_id, l.snarl.end.node_id))) 

            if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
                per_locus.append(path_in_bubble)
            else:
                locus_count=locus_count+1
                per_locus = []
                per_locus.append(path_in_bubble)
            prev_startsnarl = current_startsnarl
            prev_startsnarl_orientation = current_startsnarl_orientation
            prev_endsnarl = current_endsnarl
            prev_endsnarl_orientation = current_endsnarl_orientation
            if len(per_locus) < 3:
                locus_branch_mapping[locus_count]=per_locus
    print('The number of hets:')
    het_count= 0
    for k,v in locus_branch_mapping.items():
        if len(v) >1:
            het_count = het_count +1
    print(het_count)
    # keep branch of paths in each bubble.
    alleles_per_pos= defaultdict()
    for k,v in locus_branch_mapping.items():
        alleles_per_pos[k]=len(v)

    # both simple and complex bubbles: key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles)
    reverse_mapping= defaultdict(list)
    for k,v in locus_branch_mapping.items():
        if len(v) > 1: # more than one branch
            for i,b in enumerate(v):
                if len(b) > 0:
                    for p,j in enumerate(b):
                        reverse_mapping[j].append([k,i, len(v)]) # in complex bubbles, a node can map to multiple branches.


    # both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
    # in complex bubble, set of nodes uniquely determine the path. 
    readset=ReadSet()
    count =0
    duplicated = 0
    #TODO: consider reads with only positive score.
    with stream.open(str(gam_file), "rb") as istream:
        for data in istream:
            g = vg_pb2.Alignment()
            g.ParseFromString(data) 
            # hard-coded source id, mapping quality and other values.
            val1 = True
            val2 = False

            count1 =0
            count2=0
            #score = g.score/len(g.sequence)

            #if score > 0.2:
            #   continue
            read=Read(g.name, 0, 0, sample) # create read for each read alignment

            prev_tmp=[]
            prev_locus= -1
            locus = -1

            for i in range(0,len(g.path.mapping)-1):
            #for i in g.path.mapping: # go over the mapping in a read
                # TODO: check for forward or reverse strand, we may not need it for DAG.
                edge1 = tuple((int(g.path.mapping[i].position.name), int(g.path.mapping[i+1].position.name))) # go over nodes in a mapping
                edge2 = tuple((int(g.path.mapping[i+1].position.name), int(g.path.mapping[i].position.name))) # go over nodes in a mapping
                if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
                    if edge1 in reverse_mapping:
                        #qualities = [10]* reverse_mapping[edge1][0][2]
                        qualities= 1 
                        node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge1]] # consider (locus, branch)
                    else:
                        # qualities = [10]* reverse_mapping[edge2][0][2]
                        qualities= 1 
                        node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge2]]
                    tmp = [x for x in node_inf]
                    if prev_locus != tmp[0][0]:
                        prev_tmp = tmp
                        prev_locus = tmp[0][0]
                        
                    interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
                    if len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
                        #qualities[interset_tmp[0][1]] = 0
                        qualities= 1 
                        if i== len(g.path.mapping)-2:
                            if interset_tmp[0][1] == 0 or interset_tmp[0][1] ==1:
                                 read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
                        else:
                            next_edge1 = tuple((int(g.path.mapping[i+1].position.name), int(g.path.mapping[i+2].position.name)))
                            next_edge2 = tuple((int(g.path.mapping[i+2].position.name), int(g.path.mapping[i+1].position.name)))

                            if next_edge1 not in reverse_mapping and next_edge2 not in reverse_mapping:
                                if interset_tmp[0][1] == 0 or interset_tmp[0][1] ==1:
                                    read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)    

                        locus= interset_tmp[0][0]
                        

            readset.add(read)

    readset1=ReadSet()
    tmp_duplicated=set()
    for read in readset:
        if read.sort() ==1:
            duplicated = duplicated +1
            tmp=[]
            for variant in read:
                tmp.append(variant.position)
            #print("duplicated variant")
            x = [item for item, count in collections.Counter(tmp).items() if count > 1]
            for a in x:
                tmp_duplicated.add(a)
            continue
        else:
            tmp=[]
            for variant in read:
                tmp.append(variant.position)
            #print("duplicated variant")
            x = [item for item, count in collections.Counter(tmp).items() if count > 1]
            if len(x)>0:
                continue
            if len(read) >=4:
                tmp =[]
                for variant in read:
                   tmp.append(variant.position)
                flag=0
                for i, x in enumerate(tmp):
                    if i > 0:
                        #print(int(x - tmp[i - 1]))
                        if int(x - tmp[i - 1])>20:
                            flag=1
                            break
                if flag ==0:    
                    #print(read)
                    readset1.add(read)
            #if len(read) >=5:
             #   readset1.add(read)
    #print("length of duplicated bubbles")
    #print(tmp_duplicated)
    #print(len(list(tmp_duplicated)))


    readset1.sort()
    
    #print("duplicated")
    #print(duplicated)
    print("reads considered before read-selection")
    print(len(readset1))
    #print(readset1)
    return readset1, alleles_per_pos, locus_branch_mapping

def write_read_list(readset, bipartition, sample_components, numeric_sample_ids, output_file):
    """
    Write a list of reads that has been used for phasing to given file object.
    readset -- core.ReadSet object with reads to be written
    bipartition -- bipartition of reads, i.e. iterable with one entry from {0,1} for each read in readset
    sample_components -- a dictionary that maps each sample to its connected components

            Each component in turn is a dict that maps each variant position to a
            component, where a component is identified by the position of its
            left-most variant

    numeric_sample_ids -- core.NumericSampleIds object mapping sample names to numeric ids as stored in each read
    output_file -- file object to write to
    """
    assert len(readset) == len(bipartition)
    f = open(output_file, 'w')
    #numeric_id_to_name = numeric_sample_ids.inverse_mapping()
    for read, haplotype in zip(readset, bipartition):
        #sample = numeric_id_to_name[read.sample_id]
        sample = read.sample_id
        #print(read.sample_id)
        components = sample_components[read.sample_id]
        phaseset = components[read[0].position] + 1
        #f = open(output_file, 'w')
        print(read.name, read.source_id, sample, str(output_file)+str(phaseset), haplotype, len(read), read[0].position+1, read[-1].position+1, file=f)

def reverse_complement(seq):
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
    return "".join([seq_dict[base] for base in reversed(seq)])

    
# assumption that ped has m,f,c format and the GAMs are also in m,f,c.
# phase_input_files is a list of paths to gams files for m, f and c.

# def run_phaseg(locus_file, phase_input_files, use_ped_samples=False, read_list_filename=None):
def run_phaseg(read_list_filename,  use_ped_samples, locus_file, phase_input_files):
    #print("PIF", phase_input_files)
    recombrate=1.26
    #needs to be 5 for trio, around 15 for individual
    max_coverage = 9
    all_heterozygous = False
    distrust_genotypes = True

    readsets = dict()
    total_readsets = dict()
    for sample in [0,1,2]:
        # with timers('read_bam'):
        #print("i am here 0")
        readset, alleles_per_pos, locus_branch_mapping = vg_reader(locus_file, phase_input_files[sample], sample)
        #print("i am here 1")
        #print(readset)
        #print("i am here 1")
        total_readsets[sample] = readset
        #print(readset)
        selected_indices = readselection(readset, max_coverage)
       # print("i am here 1.2")
        selected_reads = readset.subset(selected_indices)
       # print(" i am here 1.5", len(selected_reads))
       # print("i am here 2 ")
       # print(selected_reads)
            #selected_reads = select_reads(readset, max_coverage_per_sample, preferred_source_ids = vcf_source_ids), check here for vcf_source_ids
        readsets[sample] = selected_reads
       # print("i am here 3 ")

    # Merge reads into one ReadSet (note that each Read object
    # knows the sample it originated from).
    all_reads = ReadSet()
    for sample, readset in readsets.items():
        for read in readset:
            #print("read yay", #read)
            #read.sort()
            assert read.is_sorted(), "Add a read.sort() here"
            all_reads.add(read)
        #print('test')
    #print('exited loop')
    print("allreads", len(all_reads))
    all_reads.sort()
    accessible_positions = sorted(all_reads.get_positions())
    total_reads = ReadSet()
    for sample, readset in total_readsets.items():
         for read in readset:
             assert read.is_sorted(), "Add a read.sort() here"
             total_reads.add(read)

    # # TEST STARTS HERE!!!!
   # print('hi')
    all_heterozygous = False

    recombcost = [100] * len(accessible_positions) 
    pedigree = Pedigree(NumericSampleIds())
    #genotype_likelihoods = []
    
    #for pos in accessible_positions:
     #   possible_genotypes = 3
      #  genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
    #print(genotype_likelihoods)
    genotype_likelihoods = [None if all_heterozygous else PhredGenotypeLikelihoods(0,0,0)] * len(accessible_positions)
    pedigree.add_individual('mother', [1] * len(accessible_positions), genotype_likelihoods) # all genotypes heterozygous
    pedigree.add_individual('father', [1] * len(accessible_positions), genotype_likelihoods) # all genotypes heterozygous
    pedigree.add_individual('child', [1] * len(accessible_positions), genotype_likelihoods) # all genotypes heterozygous
    pedigree.add_relationship('mother', 'father', 'child') 
    #print("ALLREADS", all_reads)
   # print ('hey again', pedigree)
    #distrust_genotypes = False
    print(all_reads)
    dp_table = PedigreeDPTable(all_reads, recombcost, pedigree, distrust_genotypes = not all_heterozygous)
    superreads_list, transmission_vector = dp_table.get_super_reads()
    #print('i am here')
    for s in superreads_list:
        print(s)
    master_block= None
    # #TEST ENDS HERE
    # If we distrusted genotypes, we need to re-determine which sites are homo-/heterzygous after phasing
    if distrust_genotypes:
        hom_in_any_sample = set()
        heterozygous_positions_by_sample = {}
        heterozygous_gts = frozenset({(0, 1), (1, 0)})
        homozygous_gts = frozenset({(0, 0), (1, 1)})
        #print(family_trios['1'][0][0])
        #fam_list = [family_trios['1'][0][0], family_trios['1'][0][1], family_trios['1'][0][2]]
        for sample, sample_superreads in zip([0,1,2], superreads_list):
            hets = set()
            for v1, v2 in zip(*sample_superreads):
                assert v1.position == v2.position
                if v1.position not in accessible_positions:
                    continue
                gt = (v1.allele, v2.allele)
                if gt in heterozygous_gts:
                    hets.add(v1.position)
                elif gt in homozygous_gts:
                    hom_in_any_sample.add(v1.position)
            heterozygous_positions_by_sample[sample] = hets
            master_block = sorted(hom_in_any_sample)
    # TODO: update master block code, for now, take None.
    overall_components = find_components(accessible_positions, all_reads, master_block, heterozygous_positions_by_sample)
    n_phased_blocks = len(set(overall_components.values()))
    components = defaultdict()
    print('homozygous in an sample master block', locus_file, master_block)
    # Superreads in superreads_list are in the same order as individuals were added to the pedigree
    f = open(str(locus_file)+".allreads", 'w')
    for sample, sample_superreads in zip([0,1,2], superreads_list):
        #superreads[sample] = sample_superreads
        # assert len(sample_superreads) == 2
        # assert sample_superreads[0].sample_id == sample_superreads[1].sample_id == numeric_sample_ids[sample]
        # # identical for all samples
        components[sample] = overall_components
        haplotag(sample_superreads, total_readsets[sample], components[sample], accessible_positions, locus_file, 1, f)

    
    if read_list_filename:
        write_read_list(all_reads, dp_table.get_optimal_partitioning(), components, {0:0, 1:1, 2:2}, read_list_filename)
        #haplotag(superreads_list[1], total_readsets[1], components[1], accessible_positions, locus_file, 1)
        #for sample, sample_superreads in zip(fam_list, superreads_list):
            #print("HEYHEY")
            #accessible_positions = sorted(readsets[sample].get_positions())
            #print(sample_superreads)
        #haplotag(sample_superreads, total_readsets[sample], components[sample], accessible_positions, locus_file, 1)
# def add_arguments(parser):
#   arg = parser.add_argument
#   # Positional arguments
#   arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
#   arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
#       help='BAM, CRAM or VCF file(s) with phase information, either through '
#           'sequencing reads (BAM/CRAM) or through phased blocks (VCF)')
#   arg('--output-read-list', metavar='FILE', default=None, dest='read_list_filename',
#       help='Write reads that have been used for phasing to FILE.')
#   arg('--use-ped-samples', dest='use_ped_samples',
#       action='store_true', default=False,
#       help='Only work on samples mentioned in the provided PED file.')

def add_arguments(parser):
   arg = parser.add_argument
   # Positional arguments
   arg('read_list_filename', metavar='FILE', help='Write reads that have been used for phasing to FILE.')
   arg('use_ped_samples', metavar='PED', help='Only work on samples mentioned in the provided PED file.')
   arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
   arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
       help='BAM, CRAM or VCF file(s) with phase information, either through '
           'sequencing reads (BAM/CRAM) or through phased blocks (VCF)')

def main(args):
    run_phaseg(**vars(args))
