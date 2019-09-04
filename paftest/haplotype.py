import sys
from collections import defaultdict, OrderedDict
import argparse
from multiprocessing import Pool
from cigarStringParser2 import snpDetector
from Bio import SeqIO
import networkx as nx
import subprocess
import clustering
from stringGraph import StringGraph, Node, Overlap, NodeCoords, Edge
import copy, operator

class Read(object):
    """docstring for Read"""
    def __init__(self, name, length, tgt_start, tgt_end, snp_coord, del_region):
        self.name = name
        self.length = length
        self.tgt_start = tgt_start
        self.tgt_end = tgt_end
        self.pos_var = {}
        self.pos_ref = {}
        for snp in snp_coord:
            self.pos_var[snp[0]] = snp[3]
            self.pos_ref[snp[0]] = snp[2]
        self.ori_var = copy.deepcopy(self.pos_var)
        self.deletionPos = set()
        for region in del_region:
            for pos in range(region[0], region[1]):
                    self.deletionPos.add(pos)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

class pileup():
    """
    Store the variant call from all query reads for one target read. And 
    apply filtering to get the true variants.
    All the coordinates included are coordinates on the target sequence.
    """
    def __init__(self, target_length, targetName):
        self.targetName = targetName
        self.target_length = target_length
        self.read_set = {}
        self.tgt_cover = {}
        self.posToVarread = defaultdict(set)
        self.posRef = dict()
        self.pos_var_count = defaultdict(dict)
        self.previous_aln = 0

    def add_read(self, read):
        assert isinstance(read, Read)
        assert read.tgt_start >= self.previous_aln, 'Input paf should be sorted by target overlapping start position'
        self.previous_aln = read.tgt_start

        if read.name in self.read_set and len(read.pos_ref) >= len(self.read_set[read.name].pos_var):
            return False
        else:
            
            self.tgt_cover[read.name] = [(read.tgt_start, 1), (read.tgt_end, -1)]
            for pos, base in read.pos_ref.items():
                if pos in self.posRef:
                    assert self.posRef[pos] == base
                else:
                    self.posRef[pos] = base
            if read.name in self.read_set:
                for pos in self.read_set[read.name].pos_var:
                    self.posToVarread[pos].remove(self.read_set[read.name])
            for pos in read.pos_var:
                self.posToVarread[pos].add(read)
            self.read_set[read.name] = read
            log(read, read.pos_var)
            return True

    def depth(self, pos):
        c = 0
        tc = []
        for r, t in self.tgt_cover.items():
            tc.extend(t)
        tc = sorted(tc)
        for i in range(len(tc)):
            if tc[i][0] < pos:
                c += tc[i][1]
            if tc[i][0] == pos:
                while i < len(tc) and tc[i][0] == pos:
                    c += 1
                    i += 1
                break
            if tc[i][0] > pos:
                break
        for readName, read in self.read_set.items():
            if pos in read.deletionPos:
                c -= 1
        return c

    def getAltVar(self):
        for pos in range(self.target_length):
            if pos in self.posToVarread:
                for read in self.posToVarread[pos]:
                    allele = read.pos_var[pos]
                    try:
                        self.pos_var_count[pos][allele] += 1
                    except KeyError:
                        self.pos_var_count[pos][allele] = 1

    def removeError(self):
        '''Filter variants by allele coverage and allele fraction'''
        tmp = defaultdict(dict)
        badPos = set()

        for pos, var_count in self.pos_var_count.items():
            Bad = 0
            good = 0
            altN = 0
            for b, n in var_count.items():
                log('looking at', pos, b, n, self.depth(pos))
                if (n / self.depth(pos) > 0.2 and n / self.depth(pos) < 0.8) and (n > 3 and self.depth(pos) - n > 3):
                    log('good')
                    good += 1
                altN += n
            refN = self.depth(pos) - altN
            if (refN / self.depth(pos) > 0.2 and refN / self.depth(pos) < 0.8) and (refN > 3 and self.depth(pos) - refN > 3):
                log('ref', refN, 'good')
                good += 1
                
            if good < 2 or self.depth(pos) <= 5:
                log('removing position', pos, 'count:', altN, 'depth:', self.depth(pos))
                del self.posToVarread[pos]
                badPos.add(pos)
            else:
                tmp[pos] = var_count.copy()
            '''
            if altN / self.depth(pos) <= 0.1 or altN / self.depth(pos) >= 0.9 or altN < 3 or self.depth(pos) - altN < 3: # THRESHOLD: allele freq to remove error
                Bad += 1
                log('removing position', pos, 'count:', altN, 'depth:', self.depth(pos))
            if Bad > 0:
                del self.posToVarread[pos]
                badPos.add(pos)
            else:
                tmp[pos] = var_count.copy()
            '''
        self.pos_var_count = tmp.copy()
        for readName, read in self.read_set.items():
            for pos in badPos:
                if pos in read.pos_var:
                    del read.pos_var[pos]

    def formSmallClusters(self):
        self.clusters = defaultdict(set)
        self.ovlp = {}
        for readName, read in self.read_set.items():
            log(read, read.pos_var, self.tgt_cover[readName][0][0], self.tgt_cover[readName][1][0])
            ordered_var = []
            nvar = 0
            for pos in range(self.target_length):
                if pos in self.pos_var_count:
                    if pos < self.tgt_cover[readName][0][0]:
                        continue
                    if pos > self.tgt_cover[readName][1][0]:
                        break
                    if pos in read.pos_var:
                        nvar += 1
                        ordered_var.append(str(pos)+'_'+read.pos_var[pos]+'_ALT')
                    else:
                        ordered_var.append(str(pos)+'_'+self.posRef[pos]+'_REF')
            self.ovlp[readName] = len(ordered_var)
            ordered_var.insert(0, nvar)
            self.clusters[tuple(ordered_var)].add(read)
        self.smallClusters = []
        for i in self.clusters:
            c = Cluster(i, self.clusters[i])
            self.smallClusters.append(c)

    def getFinalClusters(self):
        self.smallClusters = sorted(self.smallClusters, key=operator.attrgetter('nvar'))
        self.largeClusters = []
        dashed = []
        log('Basic small clusters:')
        for c in self.smallClusters:
            log(c)
            log(c.members)
        for c in self.smallClusters:
            if c.nvar > 10:
                break
            if len(c.pos_allele) == 0:
                dashed.append(c)
                continue
            if self.largeClusters == []:
                self.largeClusters.append(c)
                continue
            merged = False
            for C in self.largeClusters:
                if C.distance(c) == 0:
                    C.merge(c)
                    merged = True
                    break
                elif len(c.pos_allele) > 3 and len(C.pos_allele) > 3 and C.distance(c) <= 1:
                    C.merge(c)
                    merged = True
                    break
            if not merged:
                self.largeClusters.append(c)

        log('Final clusters')
        for c in self.largeClusters:
            log(c)
            log(c.members)

        C1 = {}
        if len(self.largeClusters) > 0:
            self.largeClusters[0].selfCheck()
            for read in self.largeClusters[0].members:
                if len(read.ori_var)/(read.tgt_end - read.tgt_start) < 0.001:
                    C1[read.name] = self.ovlp[read.name]
                else:
                    log('Discarded', read, 'because originally too many snps', len(read.ori_var), 'with overlapping length', read.tgt_end - read.tgt_start)
        dash = {}
        for c in dashed:
            for read in c.members:
                if len(read.ori_var)/(read.tgt_end - read.tgt_start) < 0.001: 
                    dash[read.name] = 0
                else:
                    log('Discarded', read, 'because originally too many snps', len(read.ori_var), 'with overlapping length', read.tgt_end - read.tgt_start)

        return C1, dash

    def getTrueHaps(self, outputCov):

        self.getAltVar()
        self.removeError()
        if outputCov:
            log('All allele position and coverage for target %s:'%self.targetName)
            log(self.pos_var_count)
            log('REF base at positions')
            log(self.posRef)
        
        self.formSmallClusters()
        C1, dash = self.getFinalClusters()
     
        return C1, dash

class Cluster(object):

    def __init__(self, ordered_var, members):
        self.nvar = ordered_var[0]
        self.pos_allele = {}
        self.pos_var = {}
        for i in ordered_var[1:]:
            pos = int(i.split('_')[0])
            allele = i.split('_')[1]
            self.pos_allele[pos] = allele
            if i.split('_')[2] == 'ALT':
                self.pos_var[pos] = allele
        #assert self.nvar == len(self.pos_allele)
        self.members = set(members)
        self.diff = {}
        self.diff_reads = defaultdict(set)
    
    def merge(self, cluster):
        assert isinstance(cluster, Cluster)
        for pos, allele in cluster.pos_allele.items():
            if pos not in self.pos_allele:
                self.pos_allele[pos] = allele
            else:
                if allele != self.pos_allele[pos]:
                    log('Allowed', pos, allele, 'when merging')
                    log(self)
                    log(cluster)
                    self.diff[pos] = allele
                    self.diff_reads[pos] = self.diff_reads[pos].union(cluster.members)
        # TODO nvar update??? Dont think is necessary for now...
        self.members = self.members.union(cluster.members)

    def distance(self, cluster):
        assert isinstance(cluster, Cluster)
        d = 0
        for pos, allele in cluster.pos_allele.items():
            if pos in self.pos_allele:
                if self.pos_allele[pos] != allele:
                    d += 1
            if pos in self.diff:
                if self.diff[pos] != allele:
                    d += 1
        return d

    def inClusterDepth(self, pos):
        depth = 0
        tgt_cover = []
        for read in self.members:
            tgt_cover.append([read.tgt_start, 1])
            tgt_cover.append([read.tgt_end, -1])
        tgt_cover = sorted(tgt_cover)
        for i in range(len(tgt_cover)):
            if tgt_cover[i][0] < pos:
                depth += tgt_cover[i][1]
            if tgt_cover[i][0] == pos:
                while i < len(tgt_cover) and tgt_cover[i][0] == pos:
                    depth += 1
                    i += 1
                break
            if tgt_cover[i][0] > pos:
                break
        for read in self.members:
            if pos in read.deletionPos:
                depth -= 1
        return depth

    def selfCheck(self):
        toRemove = set()
        for pos, readSet in self.diff_reads.items():
            if len(readSet) / self.inClusterDepth(pos) > 0.1 or len(readSet) / self.inClusterDepth(pos) < 0.9:
                toRemove = toRemove.union(readSet)
        
        self.members = self.members.difference(toRemove)

    def __str__(self):
        return 'Cluster with '+str(self.nvar)+' variants; with alleles:'+str(self.pos_allele)+'; with %d members'%(len(self.members))+'; with allowed difference at: '+str(self.diff)+'; SNP: '+str(self.pos_var)      

def align(fasta, t):
    cmd = 'minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t %d %s %s | sort -k8n -k6'%(t, fasta, fasta)
    sys.stderr.write('minimap2 cmd: %s\n'%cmd)
    sys.stderr.write('Aligning...\n')
    paf_run = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE, stdout = subprocess.PIPE)

    com = paf_run.communicate()
    returnCode = paf_run.returncode
    stdout = com[0].decode()
    stderr = com[1].decode()
    if returnCode != 0:
        sys.stderr.write('Error: minimap2 reuturn code: %d. Stderr follows: \n' % returnCode)
        sys.stderr.write(stderr)

    pafout = open(fasta.split('/')[-1]+'.paf', 'w')
    pafout.write(stdout)
    log(stderr)
    pafout.close()
    paflines = stdout.split('\n')[:-1]

    return paflines

def split_jobs(paflines):
    
    blocks = defaultdict(list)
    for line in paflines:
        if line == '':
            continue
        readname = line.split('\t')[5]
        blocks[readname].append(line)
    return blocks

def log(*content):
    out = []
    for i in content:
        out.append(str(i))
    sys.stderr.write(' '.join(out)+'\n')

def hasClipping(target_start, target_end, target_length, query_start, query_end, query_length, query_dir):
    if target_start > 20 and target_length - target_end > 20 and query_start > 20 and query_length - query_end > 20:
        return True

def variantCall(targetName, paflines, outputCov):
    query_record = dict()
    target_length = int(paflines[0].split('\t')[6])
    caller = pileup(target_length, targetName)
    target = Read(targetName, target_length, 0, target_length, [], [])
    caller.add_read(target)
    for line in paflines:
        tokens = line.split('\t')
        query_name = tokens[0]
        query_length = int(tokens[1])
        query_start = int(tokens[2])
        query_end = int(tokens[3])
        query_dir = tokens[4]
        target_name = tokens[5]
        target_start = int(tokens[7])
        target_end = int(tokens[8])
        aln_block_length = int(tokens[10])
        for field in tokens[10:]:
            if field[:5] == 'cs:Z:':
                cigar2 = field[5:]
            if field[:5] == 'cg:Z:':
                cigar1 = field[5:]
        target_coord = target_start
        if (query_end - query_start) / query_length <= 0.15: # THRESHOLD: minimum ovlp length to consider a query for variant calling
            continue
        if query_dir == '+':
            query_coord = 0
        else:
            query_coord = query_end
        snp_coord, del_region = snpDetector(cigar2, target_coord, query_coord, query_dir)
        #if len(snp_coord)/aln_block_length >= 0.04: # THRESHOLD to remove reads coming from other regions
        #    continue
        if hasClipping(target_start, target_end, target_length, query_start, query_end, query_length, query_dir):
            log('skipped', query_name, 'for clipping')
            continue
        
        read = Read(query_name, query_length, target_start, target_end, snp_coord, del_region)

        if caller.add_read(read):
            query_record[query_name] = [0, target_start, target_end, 
                                    target_length, query_start, query_end, 
                                    query_length, query_dir, target_end - target_start, cigar1]

    C1, dash = caller.getTrueHaps(outputCov)
    Llines = []
    contained = set()
    for query_name, query_info in query_record.items():
        # First check prefix-suffix.

        # Set a small window (20bp) for suffix prefix starting/end pos
        
        if query_info[4] <= 20 and query_info[6] - query_info[5] > 20 \
           and query_info[1] > 20 and query_info[3] - query_info[2] <= 20 \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] <= 20 and query_info[5] == query_info[6] \
           and query_info[1] > 20 and query_info[3] - query_info[2] <= 20 \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] <= 20 and query_info[6] - query_info[5] > 20 \
           and query_info[1] == 0 and query_info[3] - query_info[2] <= 20 \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        

        elif query_info[1] <= 20 and query_info[3] - query_info[2] > 20\
             and query_info[4] > 20 and query_info[6] - query_info[5] <= 20 \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] <= 20 and query_info[2] == query_info[3] \
             and query_info[4] > 20 and query_info[6] - query_info[5] <= 20 \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] <= 20 and query_info[3] - query_info[2] > 20 \
             and query_info[4] == 0 and query_info[6] - query_info[5] <= 20 \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]

        elif query_info[1] <= 20 and query_info[3] - query_info[2] > 20 \
             and query_info[4] <= 20 and query_info[6] - query_info[5] > 20 \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #    '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] <= 20 and query_info[2] == query_info[3] \
             and query_info[4] <= 20 and query_info[6] - query_info[5] > 20 \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] <= 20 and query_info[3] - query_info[2] > 20 \
             and query_info[4] <= 20 and query_info[5] == query_info[6] \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #    '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]

        elif query_info[4] > 20 and query_info[6] - query_info[5] <= 20 \
             and query_info[1] > 20 and query_info[3] - query_info[2] <= 20 \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] == 0 and query_info[6] - query_info[5] <= 20 \
             and query_info[1] > 20 and query_info[3] - query_info[2] <= 20 \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] > 20 and query_info[6] - query_info[5] <= 20 \
             and query_info[1] == 0 and query_info[3] - query_info[2] <= 20 \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]

        else:
            if query_info[1] <= 20 and query_info[3] - query_info[2] <= 20 and (query_name in C1 or query_name in dash):
                contained.add(target_name)
            elif query_info[4] <= 20 and query_info[6] - query_info[5] <= 20 and (query_name in C1 or query_name in dash):
                contained.add(query_name)
            continue

        if query_name in C1:
            L.append(C1[query_name])
            L.append('GREEN')
            Llines.append(L)
        elif query_name in dash:
            L.append(dash[query_name])
            L.append('DASH')
            Llines.append(L)
    return targetName, Llines, contained, caller.pos_var_count

def selectOvlp(Llines):
    node2node = defaultdict(set)
    for L in Llines:
        # Here src & tgt is defined as the 'left' read in L line and 'right' read in L line.
        # tgt here doesn't mean target in paf
        #log(L)

        src = L[0]
        tgt = L[5]
        srcdir = L[4]
        tgtdir = L[9]
        ovlp = L[10]
        srclen = L[1]
        tgtlen = L[6]
        nvar = L[-2]

        if srcdir == '+':
            node2node[src].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, tgt, 'E', tuple(L)))
        else:
            node2node[src].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, tgt, 'B', tuple(L)))
        if tgtdir == '+':
            node2node[tgt].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, src, 'B', tuple(L)))
        else:
            node2node[tgt].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, src, 'E', tuple(L)))

    Llinesout = []
    seenPair = set()
    for src, info in node2node.items():
        if len(info) == 1:
            continue
        info = sorted(list(info), reverse = True)
        maxovlp = info[0][0]
        minovlp = maxovlp * 0.5 # THRESHOLD, for selecting the best ovlpped connections
        good = []
        potential = []
        bad = []
        dash = []
        for conn in info:
            log(conn)
            if conn[3] >= 4 and ((conn[1] >= 0.5 and conn[2] >= 0.50) or conn[0] > 4000): # Before was 6000, for ragoo based, adjusted to 4000
                good.append(conn)
            elif conn[3] > 0 and ((conn[1] >= 0.50 and conn[2] >= 0.50) or conn[0] > 4000):
                potential.append(conn)
            elif conn[3] == 0 and ((conn[1] >= 0.4 and conn[2] >= 0.4) or conn[0] > 4000):
                dash.append(conn)
                 
            else:
                bad.append(conn)
        log('good', good)
        log('potential', potential)
        log('dash', dash)
        connected = []
        nB, nE, t = 0, 0, 0

        while ( nB < 4 or nE < 4 ) and t < len(good):
            if good[t][5] == 'E':

                if nE < 4:
                    nE += 1
                    connected.append(good[t])
            if good[t][5] == 'B':

                if nB < 4:
                    nB += 1
                    connected.append(good[t])
            t += 1

        t = 0
        while ( nB < 5 or nE < 5 ) and t < len(potential):
            if potential[t][5] == 'E':
                if nE < 5:
                    nE += 1
                    connected.append(potential[t])
            if potential[t][5] == 'B':
                if nB < 5:
                    nB += 1
                    connected.append(potential[t])
            t += 1

        t = 0
        while ( nB < 5 or nE < 5 ) and t < len(dash):
            if dash[t][5] == 'E':
                if nE < 5:
                    nE += 1
                    connected.append(dash[t])
            if dash[t][5] == 'B':
                if nB < 5:
                    nB += 1
                    connected.append(dash[t])
            t += 1   
        
        #if nE == 0 or nB == 0:
        #    continue
        for conn in connected:
            if (src, conn[4]) not in seenPair and (conn[4], src) not in seenPair:
                if conn[-2] == 'E':
                    Llinesout.append(list(conn[-1]))
                elif conn[-2] == 'B':
                    Llinesout.append(list(conn[-1]))
                elif conn[3] == 0:
                    Llinesout.append(list(conn[-1]))
                seenPair.add((src, conn[4]))

    return Llinesout


def transitive_reduction(Llines):
    
    nodeSet = dict()
    existingPair = set()
    dirdic = {1: True, -1: False}
    record = {}
    G = StringGraph()
    for line in Llines:
        if line[-1] == 'GREEN' or line[-1] == 'ONETWO' or line[-1] == 'DASH':
            
            srcName, srcLength, srcStart, srcEnd, srcDir = line[:5]
            tgtName, tgtLength, tgtStart, tgtEnd, tgtDir = line[5:10]

            if srcDir == '+':
                srcDir = 1
            else:
                srcDir = -1
            if tgtDir == '+':
                tgtDir = 1
            else:
                tgtDir = -1
            
            pair = (srcName, srcDir, tgtName, tgtDir)
            reverse_pair = (tgtName, -1 * tgtDir, srcName, -1 * srcDir)
            if reverse_pair not in existingPair and pair not in existingPair:
                record[(srcName, tgtName)] = line
                existingPair.add(pair)
                if srcName not in nodeSet:
                    nodeSet[srcName] = Node(srcName, srcLength)
                    G.addNode(nodeSet[srcName])
                if tgtName not in nodeSet:
                    nodeSet[tgtName] = Node(tgtName, tgtLength)
                    G.addNode(nodeSet[tgtName])
                G.addOverlap(Overlap(NodeCoords(nodeSet[srcName], srcStart, srcEnd, dirdic[srcDir]), 
                                     NodeCoords(nodeSet[tgtName], tgtStart, tgtEnd, dirdic[tgtDir])))
    G.transitiveReduce()
    final_L = []
    reducted = set()
    for e in G.edges():
        srcName = e.src.id
        tgtName = e.target.id
        if (srcName, tgtName) in record:
            L = record[(srcName, tgtName)]
        elif (tgtName, srcName) in record:
            L = record[(tgtName, srcName)]
        else:
            log('transitive reduction generates an inexisting pair???', srcName, tgtName)
        if (srcName, tgtName) not in reducted and (tgtName, srcName) not in reducted:
            final_L.append(L)
            reducted.add((srcName, tgtName))

    return final_L
    '''
    walk = G.findUnipaths()
    for path in walk:
        for e in path.edges:
            final_L.append(edgeObj2Lline(e))
    '''

def removeTip(graph):
    nodeB = defaultdict(set)
    nodeE = defaultdict(set)
    outL = []
    for L in graph:
        src = L[0]
        tgt = L[5]
        srcdir = L[4]
        tgtdir = L[9]
        ovlp = L[10]
        srclen = L[1]
        tgtlen = L[6]
        nvar = L[-2]

        # 'E' stands for connecting to the original end of 'key' read
        # 'B' stands for connecting to the original beginning of 'key' read
        if srcdir == '+':
            nodeE[src].add((ovlp, tgt, tuple(L)))
        else:
            nodeB[src].add((ovlp, tgt, tuple(L)))
        if tgtdir == '+':
            nodeB[tgt].add((ovlp, src, tuple(L)))
        else:
            nodeE[tgt].add((ovlp, src, tuple(L)))

    tips = set()
    for node in nodeB:
        if node not in nodeE:
            tips.add(node)
    for node in nodeE:
        if node not in nodeB:
            tips.add(node)
    out = []
    log('Single tip reads:')
    log(tips)
    for L in graph:
        if L[0] in tips or L[5] in tips:
            continue
        out.append(L)
    return out

def oneConn(graph):
    nodeB = defaultdict(set)
    nodeE = defaultdict(set)
    outL = []
    for L in graph:
        src = L[0]
        tgt = L[5]
        srcdir = L[4]
        tgtdir = L[9]
        ovlp = L[10]
        srclen = L[1]
        tgtlen = L[6]
        nvar = L[-2]
        # 'E' stands for connecting to the original end of 'key' read
        # 'B' stands for connecting to the original beginning of 'key' read
        if srcdir == '+':
            nodeE[src].add((ovlp, tgt, tuple(L)))
        else:
            nodeB[src].add((ovlp, tgt, tuple(L)))
        if tgtdir == '+':
            nodeB[tgt].add((ovlp, src, tuple(L)))
        else:
            nodeE[tgt].add((ovlp, src, tuple(L)))
    
    for node, conns in nodeB.items():
        conns = list(conns)
        if len(conns) == 1:
            outL.append(conns[0][2])
        else:
            conns = sorted(conns)
            outL.append(conns[-1][2])
    for node, conns in nodeE.items():
        conns = list(conns)
        if len(conns) == 1:
            outL.append(conns[0][2])
        else:
            conns = sorted(conns)
            outL.append(conns[-1][2])
    
    return outL

def gfaToNX(Llines):
    g = nx.Graph()
    for line in Llines:
        if line[-1] == 'GREEN' or line[-1] == 'DASH':
            n1 = line[0]
            n2 = line[5]
            g.add_node(n1)
            g.add_node(n2)
            g.add_edge(n1, n2)
    return g    

def draw2(fasta, Llines, out, haveDoneTR):

    consideredReads = set()
    existingPair = set()
    Lliness = copy.deepcopy(Llines)
    for Lline in Lliness:
        if haveDoneTR:
            consideredReads.add(Lline[0])
            consideredReads.add(Lline[2])
            Lline.insert(0, 'L')
            Lline[-2], Lline[-3], Lline[-4] = str(Lline[-2]), str(Lline[-3]), str(Lline[-4])
            Lline[-5] = str(Lline[-5])+'M'
            content = '\t'.join(Lline)+'\n'
        else:
            consideredReads.add(Lline[0])
            consideredReads.add(Lline[5])
            content = 'L\t' + '\t'.join([Lline[0], Lline[4], Lline[5], Lline[9], str(Lline[10])+'M', str(Lline[1]), str(Lline[6]), str(Lline[-2]), Lline[-1]]) + '\n'
        
        out.write(content)
    for record in SeqIO.parse(fasta, 'fasta'):
        if str(record.id) not in consideredReads:
            continue
        out.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
    out.close()

def removeWrongEdge(Llines, ground_truths, out):
    p1 = set(open(ground_truths[0], 'r').read().split('\n')[:-1])
    p2 = set(open(ground_truths[1], 'r').read().split('\n')[:-1])
    filtered = []
    for line in Llines:
        if line[-1] != 'GREEN':
            filtered.append(line)
        else:
            if line[0] in p1 and line[5] in p2:
                line[-1] = 'WRONG'
                
            if line[0] in p2 and line[5] in p1:
                line[-1] = 'WRONG'
                
            filtered.append(line)
    return filtered

def outputVariant(out, calls):
    if out != None:
        varout = open(out, 'w')
        varout.write('all_var = {')
        for tn, var in calls.items():
            varout.write('"'+tn+'": '+' '.join(str(var).split()[2:])[:-1]+',\n')
        varout.write('}')
        varout.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hardcodeTest', action = 'store_true', 
                         help = 'Run a hardcoded test, all other input might stil be required for convenience in writing the scripts, but they won\'t be touched at all.')
    parser.add_argument('-p', '--paf', metavar = 'PAF', type = str, 
                         help = 'PAF file, should be consistent with the input FASTA file; if not given, we will run minimap2 on the input fasta and generate one.')
    parser.add_argument('-o', '--output', metavar = 'GFA', type = str, required = True, 
                         help = 'Output to FILE, default to stdout. ')
    parser.add_argument('fasta', metavar = 'FASTA', type = str,
                        help = 'Input all read records in one FASTA file.')
    parser.add_argument('-t', '--threads', metavar = 'INT', type = int, default = 4, 
                         help = 'Maximum number of threads to use. [4]')
    parser.add_argument('--ground-truth', metavar = 'FILE', nargs = 2, 
                         help = '2 ground truth file, each has read names in lines.')
    parser.add_argument('-l', '--partition-list', metavar = 'FILE', type = str, 
                        help = 'sprcify to output the clustering of reads. reads in one continuous chain will be printed in one line dilimited by space')
    parser.add_argument('--het', metavar = 'FLOAT', type = float, default = 0.1,
                         help = 'Heterozygosity rate, set as a filter for abnormal alignment. [0.1]')
    parser.add_argument('-c', '--log-coverage', action = 'store_true', 
                         help = 'Output the depth of SNP calls in the log')
    parser.add_argument('-v', '--output-variant', metavar = '.PY', type = str, 
                         help = 'Output variant calls to a Python script to import.')
    parser.add_argument('-V', '--input-variant', metavar = '.PY', type = str, 
                         help = 'Import variants calls from python script')

    args = parser.parse_args()

    if args.hardcodeTest:
        runHCtest()
        exit()

    if args.paf == None:
        paflines = align(args.fasta, args.threads)
    else:
        log('Using given paf file rather than align again')
        paflines = open(args.paf, 'r').read().split('\n')[:-1]
    blocks = split_jobs(paflines)
    
    if args.input_variant != None:
        log('Sorry input_variant function is not finished yet...')
        exit()
        getattr(__import__(args.input_variant, fromlist=['all_var']), 'all_var')
        #from args.input_variant import all_var

    p = Pool(args.threads)

    processes = []
    for target, paflines in blocks.items():
        processes.append(p.apply_async(variantCall, (target, paflines, args.log_coverage)))

    p.close()
    p.join()

    Llines = []
    contained = set()
    calls = {}
    for proc in processes:
        targetName, Llines_block, contained_block, calls_block = proc.get()
        Llines.extend(Llines_block)
        contained = contained.union(contained_block)
        calls[targetName] = calls_block

    co = open('.'.join(args.output.split('.')[:-1]) + '.contained.reads', 'w')
    for i in contained:
        co.write(i+'\n')
    co.close()

    outputVariant(args.output_variant, calls)

    trueConn = []
    seenConn = dict()
    count = 0
    if len(blocks) > 0:
        for L in Llines:
            if L[0] in contained or L[5] in contained:
                if L[0] in contained:
                    log(L[0], 'is contained')
                else:
                    log(L[5], 'is contained')
                log('skipped edge for containing', L[0], L[5])
                continue
            if (L[0], L[5]) not in seenConn and (L[5], L[0]) not in seenConn:
                #log('L not seen before')
                seenConn[(L[0], L[5])] = L # ovlp length and color
            else:
                if (L[0], L[5]) in seenConn:
                    if L[-2] != seenConn[(L[0], L[5])][-2]:
                        #log('True variant numbers differ', L[-2], seenConn[(L[0], L[5])][-2])
                        count += 1
                        if L[-2] < seenConn[(L[0], L[5])][-2]:
                            trueConn.append(L)
                        if L[-2] > seenConn[(L[0], L[5])][-2]:
                            trueConn.append(seenConn[(L[0], L[5])])
                        del seenConn[(L[0], L[5])]
                        continue
                    if L[-3] > seenConn[(L[0], L[5])][-3]:
                        if L[-1] == 'GREEN' or L[-1] == 'DASH' or L[-1] == 'ONETWO':
                            trueConn.append(L)
                    elif L[-3] < seenConn[(L[0], L[5])][-3]:
                        if seenConn[(L[0], L[5])][-1] == 'GREEN' or seenConn[(L[0], L[5])][-1] == 'DASH' or seenConn[(L[0], L[5])][-1] == 'ONETWO':
                            trueConn.append(seenConn[(L[0], L[5])])
                    del seenConn[(L[0], L[5])]
                elif (L[5], L[0]) in seenConn:
                    if L[-2] != seenConn[(L[5], L[0])][-2]:
                        #log('True variant numbers differ', L[-2], seenConn[(L[5], L[0])][-2])
                        count += 1
                        if L[-2] < seenConn[(L[5], L[0])][-2]:
                            trueConn.append(L)
                        if L[-2] > seenConn[(L[5], L[0])][-2]:
                            trueConn.append(seenConn[(L[5], L[0])])
                        del seenConn[(L[5], L[0])]
                        continue
                    if L[-3] > seenConn[(L[5], L[0])][-3]:
                        if L[-1] == 'GREEN' or L[-1] == 'DASH' or L[-1] == 'ONETWO':
                            trueConn.append(L)
                    elif L[-3] < seenConn[(L[5], L[0])][-3]:
                        if seenConn[(L[5], L[0])][-1] == 'GREEN' or seenConn[(L[5], L[0])][-1] == 'DASH' or seenConn[(L[5], L[0])][-1] == 'ONETWO':
                            trueConn.append(seenConn[(L[5], L[0])])
                    del seenConn[(L[5], L[0])]
        for pair, L in seenConn.items():
            #log(pair, 'happens only once', L)
            if L[-1] == 'GREEN' or L[-1] == 'DASH' or L[-1] == 'ONETWO':
                trueConn.append(L)
        Llines = trueConn
    #log(count, 'cases have different number of truevars')

    obefore = open('.'.join(args.output.split('.')[:-1]) + '.beforePostProc.gfa', 'w')
    draw2(args.fasta, Llines, obefore, False)

    #out = open(args.output, 'w')
    if args.ground_truth != None:
        log('Trying to remove wrong connection basing on provided ground truth', args.ground_truth)
        Llines = removeWrongEdge(Llines, args.ground_truth, out)
    #draw2(args.fasta, Llines, out, False)
 
    Llines = selectOvlp(Llines)
    #redName = '.'.join(args.output.split('.')[:-1]) + '.proc.gfa'
    #out3 = open(redName, 'w')
    #draw2(args.fasta, Llines, out3, False)

    Llines = transitive_reduction(Llines)
    Llines = removeTip(Llines)
    Llines = oneConn(Llines)

    redName = '.'.join(args.output.split('.')[:-1]) + '.reducted.gfa'
    out2 = open(redName, 'w')
    draw2(args.fasta, Llines, out2, False)
    #exit()
    if args.partition_list != None:
        nxGraph = gfaToNX(Llines)
        subgraphs = clustering.divide(nxGraph, args.partition_list)
    #subgraphs = clustering.divide(nxGraph)
    #print(len(subgraphs))

if __name__ == '__main__':
    main()
