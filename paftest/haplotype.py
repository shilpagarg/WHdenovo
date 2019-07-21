import sys
from collections import defaultdict, OrderedDict
import argparse
from multiprocessing import Pool
from cigarStringParser2 import snpDetector
from Bio import SeqIO
import networkx as nx
import subprocess
import clustering

class pileup():
    """
    Store the variant call from all query reads for one target read. And 
    apply filtering to get the true variants.
    All the coordinates included are coordinates on the target sequence.
    """
    def __init__(self, target_length, targetName):
        self.targetName = targetName
        self.target_length = target_length
        self.read_var = {self.targetName: {}}
        #self.read_var = defaultdict(dict)
        self.read_order = [self.targetName]
        self.mode = 'normal'
        self.tgt_cover = {self.targetName: [(0, 1), (target_length, -1)]}
        self.posToVarread = defaultdict(set)
        self.posRef = dict()
        self.pos_var_count = defaultdict(dict)
        self.maxCov = 0
        self.copyNum = 1
        self.trueHapReads = set()
        self.badPos = set()
        self.previous_aln = 0
        self.removed = set()
        self.read_length = {self.targetName: target_length}
        self.read_del = dict()

    def add_snp(self, query_name, target_start, target_end, snp_coord, query_length, del_region):
        assert target_start >= self.previous_aln, 'Input paf should be sorted by target overlapping position'
        self.previous_aln = target_start
        if query_name in self.read_var and len(snp_coord) < len(self.read_var[query_name]):
            log('Updating', query_name, 'records', [i[0] for i in snp_coord])
            for pos in self.read_var[query_name]:
                self.posToVarread[pos].remove(query_name) 
            self.read_var[query_name] = dict()
            self.read_del[query_name] = []
            for r in del_region:
                for pos in range(r[0], r[1]):
                    self.read_del[query_name].append(pos)
            for snp in snp_coord:
                self.read_var[query_name][snp[0]] = snp[3]
                self.posToVarread[snp[0]].add(query_name)
                if snp[0] in self.posRef:
                    assert self.posRef[snp[0]] == snp[2], "REF base inconsistent, %s vs %s"%(self.posRef[snp[0]], snp[2])
                self.posRef[snp[0]] = snp[2]
            self.tgt_cover[query_name] = [(target_start, 1), (target_end, -1)]
            return True
        elif query_name not in self.read_var:
            self.read_var[query_name] = dict()
            log(query_name, 'first here adding', [i[0] for i in snp_coord])
            for snp in snp_coord:
                self.read_var[query_name][snp[0]] = snp[3]
                self.posToVarread[snp[0]].add(query_name)
                if snp[0] in self.posRef:
                    assert self.posRef[snp[0]] == snp[2], "REF base inconsistent, %s vs %s"%(self.posRef[snp[0]], snp[2])
                self.posRef[snp[0]] = snp[2]
            self.read_order.append(query_name)
            self.read_del[query_name] = []
            for r in del_region:
                for pos in range(r[0], r[1]):
                    self.read_del[query_name].append(pos)
            self.read_length[query_name] = query_length
            self.tgt_cover[query_name] = [(target_start, 1), (target_end, -1)]
            return True
        else:
            log(query_name, 'not added???')
            log(len(snp_coord), self.read_var[query_name].keys())
            return False

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
        for read, del_pos in self.read_del.items():
            if pos in del_pos:
                c -= 1
        return c

    def tgtAvgCov(self):
        total_depth = 0
        for i in range(self.target_length):
            total_depth += self.depth(i)
        return total_depth / self.target_length

    def getAltVar(self):
        self.pos_var_count = defaultdict(dict)
        for pos in range(self.target_length):
            try:
                reads = self.posToVarread[pos]
            except KeyError:
                continue
            for query_name in reads:
                allele = self.read_var[query_name][pos]
                try:
                    self.pos_var_count[pos][allele] += 1
                except KeyError:
                    self.pos_var_count[pos][allele] = 1

    def removeError(self):
        '''Second filter variants by allele coverage and allele fraction'''
        tmp = defaultdict(dict)
        for pos, var_count in self.pos_var_count.items():
            Bad = 0
            altN = 0
            for b, n in var_count.items():
                altN += n
                #if n / self.depth(pos) <= 0.06 or n / self.depth(pos) >= 0.94 or n <= 2 or self.depth(pos) - n <= 2: #//changed
                #    Bad += 1
            if altN / self.depth(pos) <= 0.2 or altN / self.depth(pos) >= 0.8 or altN <= 3 or self.depth(pos) - altN <= 3: # THRESHOLD: allele freq to remove error
                Bad += 1
            if Bad == 0:
                tmp[pos] = var_count.copy()
            else:
                log('removing var at pos', pos, 'depth', self.depth(pos), var_count)
                self.badPos.add(pos)
        self.pos_var_count = tmp.copy()
        for read in self.read_var:
            for pos in self.badPos:
                if pos in self.read_var[read]:
                    del self.read_var[read][pos]

    def getTargetHap(self, Ta = 3, Tr = 2): # changed, Ta = 4 before
        green = dict()
        red = dict()
        culprit = dict()
        dash = dict()
        for read, pos_var in self.read_var.items():
            d = len(pos_var)
            ovlp = 0
            for pos in self.pos_var_count:
                if pos >= self.tgt_cover[read][0][0] and pos <= self.tgt_cover[read][1][0]:
                    ovlp += 1
            assert d <= ovlp, 'Differences are more than true variants'
            if ovlp == 0:
                log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'No true var in ovlp, dash')
                dash[read] = ovlp
                continue
            #if ovlp == 1:
            #    if d == 0:
            #        green.append(read)
            #    else:
            #        red.append(read)
            #    log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, '1 true var in ovlp')
            #    continue
            a = ovlp - d
            if a >= 200:
                continue
            if d > 10:
                continue
            if ovlp >= Ta:
                if d != 0 and (a/d) >= Tr:
                    green[read]=ovlp
                    log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'disagree', d, 'agree', a, 'Considered')
                elif d == 0:
                    green[read]=ovlp
                    log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'disagree', d, 'agree', a, 'Considered')
                else:
                    red[read]=ovlp
                    log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'disagree', d, 'agree', a, 'Discarded')
            else:
                #if d == 0:
                #    green[read]=ovlp
                #else:
                #    red[read] = ovlp
                culprit[read]=ovlp
                log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'disagree', d, 'agree', a, 'Discarded')
        return green, red, dash, culprit

    def hamming(self, read1, read2):
        ovL = min(self.tgt_cover[read1][1][0], self.tgt_cover[read2][1][0]) - max(self.tgt_cover[read1][0][0], self.tgt_cover[read2][0][0])
        ovV = 0
        d = 0
        for pos in range(max(self.tgt_cover[read1][0][0], self.tgt_cover[read2][0][0]), min(self.tgt_cover[read1][1][0], self.tgt_cover[read2][1][0])):
            if pos in self.pos_var_count:
                ovV += 1
                if pos in self.read_var[read1] and pos in self.read_var[read2]:
                    if self.read_var[read1][pos] != self.read_var[read2][pos]:
                        d += 1
                elif pos in self.read_var[read1] and pos not in self.read_var[read2]:
                    d += 1
                elif pos not in self.read_var[read1] and pos in self.read_var[read2]:
                    d += 1
        return d, ovL, ovV

    def cluster(self, maxD = 1, minOLR = 0):
        clusters = []
        culprit = dict()
        dash = dict()
        minOL = self.target_length * minOLR
        read_ov_dict = {self.targetName: len(self.pos_var_count)}
        for read, pos_var in self.read_var.items():
            if clusters == []:
                # first add the target
                clusters.append([read])
                continue
            d = len(pos_var)
            ovlp = 0
            for pos in self.pos_var_count:
                if pos >= self.tgt_cover[read][0][0] and pos <= self.tgt_cover[read][1][0]:
                    ovlp += 1
            read_ov_dict[read] = ovlp
            assert d <= ovlp, 'Differences are more than true variants'
            if ovlp == 0:
                log(read, 'target:', self.targetName, 'from', self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0], '(%d)'%(self.tgt_cover[read][1][0] - self.tgt_cover[read][0][0]), pos_var, 'No true var in ovlp, dash')
                dash[read] = ovlp
                continue
            #if ovlp < 3:
            #    culprit[read] = ovlp
            #    continue
            clustered = False
            for clust in clusters:
                matchMember = 0
                disagree = 0
                for read2 in clust:
                    d, ovL, ovV = self.hamming(read, read2) # Distance, ov_length, ov_Var
                    log('Distance between', read, read2, d, ovL, ovV, self.read_var[read], self.read_var[read2])
                    if ovlp == 1 or ovlp == 2:
                        if ovV >= 1 and d <= 0 and ovL > minOL:
                            matchMember += 1
                        elif d > 0:
                            disagree += 1
                    else:
                        if ovV >= 3 and d <= maxD and ovL > minOL: # Having at least one read overlapping with at least 3 truevar
                            matchMember += 1
                        elif d > maxD:
                            disagree += 1
                if matchMember > 0 and disagree == 0:
                    clust.append(read)
                    clustered = True
                    break
            if not clustered:
                clusters.append([read])

        #log('Final clusters for target', self.targetName, clusters)
        log(len(clusters), 'clusters')
        
        log('one or two:', culprit)
        #picker = []
        two = []
        for i in range(len(clusters)):
            log('cluster', i, 'for target', self.targetName, self.clusterStats(clusters[i]), len(clusters[i]))
            log(clusters[i])
            start, end, var, het = self.clusterStats(clusters[i])
            if ovlp == 2:
                two.append((len(var), i))
        if len(two) > 1:
            two = sorted(two)
            C1 = clusteres[two[-1][1]]
            log('In the case of only 2 var. Picked cluster', two[-1][1], 'for target', self.targetName, 'with', len(C1), 'reads inside')
        else:
            C1 = clusters[0]
            log('finally chose cluster 0', 'for target', self.targetName, 'with', len(C1), 'reads inside')
        C1dict = dict()
        for read in C1:
            C1dict[read] = read_ov_dict[read]
            log(read, 'total overlapping var', C1dict[read])
        
        # TODO then output the valid cluster(s?)

        return C1dict, dash, culprit

    def clusterStats(self, cluster):
        start = self.tgt_cover[cluster[0]][0][0]
        end = self.tgt_cover[cluster[0]][1][0]
        var = set(self.read_var[cluster[0]])
        het = 0
        #self.tgt_cover[read][0][0], 'to', self.tgt_cover[read][1][0]
        for i in range(1, len(cluster)):
            start = min(start, self.tgt_cover[cluster[i]][0][0])
            end = max(end, self.tgt_cover[cluster[i]][1][0])
            var = var.union(set(self.read_var[cluster[i]]))
        het = len(var) / (end - start)
        return start, end, var, het

    def isPS(self, read1, read2):
        '''is Prefix-Suffix'''
        if self.tgt_cover[read1][1][0] > self.tgt_cover[read2][1][0]:
            if self.tgt_cover[read1][1][0] - self.read_length[read1] > self.tgt_cover[read2][1][0] - self.read_length[read2]:
                return True
        if self.tgt_cover[read1][1][0] < self.tgt_cover[read2][1][0]:
            if self.tgt_cover[read1][1][0] - self.read_length[read1] < self.tgt_cover[read2][1][0] - self.read_length[read2]:
                return True

    def RGalgo(self, green, red, culprit):
        # Connections within C1 (green)
        Green = nx.Graph()
        for i in green:
            for j in green:
                if i == j:
                    continue
                if self.isPS(i, j):
                    Green.add_edge(i, j)

        # Connections within C2 (red)
        for i in red:
            for j in red:
                if i == j:
                    continue
                if self.isPS(i, j):
                    Green.add_edge(i, j)
        
        # Red edges across C1 and C2
        Red = nx.Graph()
        for i in green:
            for j in red:
                if self.isPS(i, j):
                    Red.add_edge(i,j)

        log(Green.degree)
        log(Red.degree)
        log()

        for read in culprit:
            Green2 = Green.copy()
            Red2 = Red.copy()
            for node in green:
                if self.isPS(read, node):
                    Green2.add_edge(read, node)
            for node in red:
                if self.isPS(read, node):
                    Red2.add_edge(read, node)
            log('Tried to put', read, 'in C1')
            log(Green2.degree)
            log(Red2.degree)
            log()
        return green, red

    def getTrueHaps(self, outputCov):
        self.getAltVar()
        self.removeError()
        if outputCov:
            log('All allele position and coverage for target %s:'%self.targetName)
            log(self.pos_var_count)
            log('REF base at positions')
            log(self.posRef)
        log('Total var pos after clean up for target', self.targetName, len(self.pos_var_count.keys()), self.pos_var_count.keys())
        C1, dash, culprit = self.cluster()
        
        return C1, dash, culprit

        #green, red = self.RGalgo(C1, C2, culprit)
        #dash = [] # TODO HARDCODED
        #return [], [], [] # TODO HARDCODED
        #return green, red, dash

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

def variantCall(targetName, paflines, outputCov):
    #Pos_varCount = defaultdict(dict)
    #Pos_Ref = dict()
    # keys are positions on target, 
    # values are dictionaries where keys are alleles values are their count
    query_record = dict()
    target_length = int(paflines[0].split('\t')[6])
    caller = pileup(target_length, targetName)
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
        #cigar1 = tokens[21][5:]
        #cigar2 = tokens[22][5:]
        target_coord = target_start
        log('Query', query_name, 'overlap length', query_end - query_start, 'read length', query_length, (query_end - query_start) / query_length)
        if (query_end - query_start) / query_length <= 0.15: # THRESHOLD: minimum ovlp length to consider a query for variant calling
            log('Query', query_name, 'not considered for overlapping', query_end - query_start, query_length, (query_end - query_start) / query_length)
            continue
        if query_dir == '+':
            query_coord = 0
        else:
            query_coord = query_end
        snp_coord, del_region = snpDetector(cigar2, target_coord, query_coord, query_dir)
        log('HHII Query', query_name, snp_coord, del_region)
        #if len(caller.read_var[query_name]) == 0 or len(snp_coord) < len(caller.read_var[query_name]):
        #    log('Read name:', query_name, snp_coord)
        if caller.add_snp(query_name, target_start, target_end, snp_coord, query_length, del_region):
            query_record[query_name] = [0, target_start, target_end, 
                                    target_length, query_start, query_end, 
                                    query_length, query_dir, target_end - target_start, cigar1]
    #green, red, dash = caller.getTrueHaps(outputCov)
    C1, dash, culprit = caller.getTrueHaps(outputCov)
    Llines = []
    for query_name, query_info in query_record.items():
        # First check prefix-suffix.
        if query_info[4] == 0 and query_info[5] < query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] == query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] < query_info[6] \
           and query_info[1] == 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]

        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] <= query_info[3] \
             and query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] >= 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]

        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] == 0 and query_info[5] < query_info[6] \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #    '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] == query_info[3] \
             and query_info[4] == 0 and query_info[5] < query_info[6] \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #     '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] == 0 and query_info[5] == query_info[6] \
             and query_info[7] == '-':
            # query-suffix '-' <==> prefix-target
            #L = [query_name, query_info[7], target_name, 
            #    '+', str(query_info[9])]
            L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_info[8], str(query_info[9])]

        elif query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[1] > 0 and query_info[2] == query_info[3] \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] == query_info[6] \
             and query_info[1] > 0 and query_info[2] == query_info[3] \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]
        elif query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[7] == '-':
            # target-suffix  <==> prefix-query '-'
            #L = [target_name, '+', query_name, 
            #     query_info[7], str(query_info[9])]
            L = [target_name, query_info[3], query_info[1], query_info[2], '+', 
                 query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
                 query_info[8], str(query_info[9])]

        else:
            continue
            #L = [query_name, query_info[6], query_info[4], query_info[5], query_info[7], 
            #     target_name, query_info[3], query_info[1], query_info[2], '+', 
            #     query_info[8], str(query_info[9])]
            #log('what is this???', L)
        
        if query_name in C1:
            L.append(C1[query_name])
            L.append('GREEN')
            Llines.append(L)
        elif query_name in dash:
            L.append(dash[query_name])
            L.append('DASH')
            Llines.append(L)
        elif query_name in culprit:
            L.append(culprit[query_name])
            L.append('ONETWO')
            Llines.append(L)
    for i in Llines:
        log(i)
    return Llines
    '''
        if query_name in green:
            L.append(green[query_name])
            L.append('GREEN')
        elif query_name in dash:
            L.append(dash[query_name])
            L.append('DASH')
        elif query_name in red:
            L.append(red[query_name])
            L.append('RED')
        else:
            continue

        Llines.append(L)
    return Llines
    '''
def gfaToNX(Llines):
    g = nx.Graph()
    for line in Llines:
        if line[-1] != 'GREEN':
            continue
        n1 = line[0]
        n2 = line[5]
        g.add_node(n1)
        g.add_node(n2)
        g.add_edge(n1, n2)

    return g

def draw2(fasta, Llines, out, haveDoneTR):

    consideredReads = set()
    existingPair = set()
    for Lline in Llines:
        if haveDoneTR:
            consideredReads.add(Lline.split('\t')[1])
            consideredReads.add(Lline.split('\t')[3])
            content = Lline+'\n'
        else:
            consideredReads.add(Lline[0])
            consideredReads.add(Lline[5])
            content = 'L\t' + '\t'.join([Lline[0], Lline[4], Lline[5], Lline[9], str(Lline[10])+'M', str(Lline[-2]), Lline[-1]]) + '\n'
        
        out.write(content)
    for record in SeqIO.parse(fasta, 'fasta'):
        if str(record.id) not in consideredReads:
            continue
        out.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
    out.close()

def edgeObj2Lline(edge):
    if not isinstance(edge, Edge):
        raise TypeError('Check input object type, I need a stringGraph.Edge object')
    else:
        srcName = edge.src.id
        tgtName = edge.target.id
        if edge.srcEnd == 'E':
            srcDir = '+'
        else:
            srcDir = '-'
        if edge.targetEnd == 'B':
            tgtDir = '+'
        else:
            tgtDir = '-'

        return 'L\t%s\t%s\t%s\t%s\t%s'%(srcName, srcDir, tgtName, tgtDir, str(edge.overlap.getLength())+'M')

def transitive_reduction(Llines):
    
    nodeSet = dict()
    existingPair = set()

    G = StringGraph()
    for line in Llines:
        if line[-1] == 'RED':
            continue
        srcName, srcLength, srcStart, srcEnd, srcDir = line[:5]
        tgtName, tgtLength, tgtStart, tgtEnd, tgtDir = line[5:10]
        if srcDir == '+':
            srcDir = True
        else:
            srcDir = False
        if tgtDir == '+':
            tgtDir = True
        else:
            tgtDir = False
        
        pair = (srcName, srcDir, tgtName, tgtDir)
        if srcDir and tgtDir:
            reverse_pair = (tgtName, tgtDir, srcName, srcDir)
        else:
            reverse_pair = (tgtName, srcDir, srcName, tgtDir)
        if reverse_pair not in existingPair and pair not in existingPair:
            existingPair.add(pair)
            if srcName not in nodeSet:
                nodeSet[srcName] = Node(srcName, srcLength)
                G.addNode(nodeSet[srcName])
            if tgtName not in nodeSet:
                nodeSet[tgtName] = Node(tgtName, tgtLength)
                G.addNode(nodeSet[tgtName])
            G.addOverlap(Overlap(NodeCoords(nodeSet[srcName], srcStart, srcEnd, srcDir), 
                                 NodeCoords(nodeSet[tgtName], tgtStart, tgtEnd, tgtDir)))

    G.transitiveReduce()

    final_L = []
    walk = G.findUnipaths()
    for path in walk:
        for e in path.edges:
            final_L.append(edgeObj2Lline(e))

    return final_L

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

def main():
    
    parser = argparse.ArgumentParser()
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
    args = parser.parse_args()

    if args.paf == None:
        paflines = align(args.fasta, args.threads)
    else:
        log('Using given paf file rather than align again')
        paflines = open(args.paf, 'r').read().split('\n')[:-1]
    blocks = split_jobs(paflines)
    
    p = Pool(args.threads)

    processes = []
    for target, paflines in blocks.items():
        processes.append(p.apply_async(variantCall, (target, paflines, args.log_coverage)))

    p.close()
    p.join()

    Llines = []
    for proc in processes:
        Llines_block = proc.get()
        Llines.extend(Llines_block)
    log('Here', len(Llines), 'edges were made.')
    trueConn = []
    seenConn = dict()
    count = 0
    if len(blocks) > 1:
        log('Working on all Llines now')
        for L in Llines:
            if (L[0], L[5]) not in seenConn and (L[5], L[0]) not in seenConn:
                log('L not seen before')
                seenConn[(L[0], L[5])] = L # ovlp length and color
            else:
                if (L[0], L[5]) in seenConn:
                    if L[-2] != seenConn[(L[0], L[5])][-2]:
                        log('True variant numbers differ', L[-2], seenConn[(L[0], L[5])][-2])
                        count += 1
                        if L[-2] < seenConn[(L[0], L[5])][-2]:
                            trueConn.append(L)
                        if L[-2] > seenConn[(L[0], L[5])][-2]:
                            trueConn.append(seenConn[(L[0], L[5])])
                        del seenConn[(L[0], L[5])]
                        continue
                    if L[-1] == 'RED' and seenConn[(L[0], L[5])][-1] == 'RED':
                        log(L, 'both red')
                    if L[-1] == 'DASH' and seenConn[(L[0], L[5])][-1] == 'DASH':
                        log(L, 'both dash')
                    if L[-1] == 'DASH' and seenConn[(L[0], L[5])][-1] == 'RED':
                        log(L, 'one dash one red')
                    if L[-1] == 'RED' and seenConn[(L[0], L[5])][-1] == 'DASH':
                        log(L, 'one dash one red')
                    if L[-3] > seenConn[(L[0], L[5])][-3]:
                        if L[-1] == 'GREEN' or L[-1] == 'DASH':
                            trueConn.append(L)
                        if seenConn[(L[0], L[5])][-1] == 'GREEN' and L[-1] != 'GREEN':
                            log(L, 'bad connection dominates.')
                    elif L[-3] < seenConn[(L[0], L[5])][-3]:
                        if seenConn[(L[0], L[5])][-1] == 'GREEN' or seenConn[(L[0], L[5])][-1] == 'DASH':
                            trueConn.append(seenConn[(L[0], L[5])])
                        if L[-1] == 'GREEN' and seenConn[(L[0], L[5])][-1] != 'GREEN':
                            log(L, 'bad connection dominates.')
                    del seenConn[(L[0], L[5])]
                elif (L[5], L[0]) in seenConn:
                    if L[-2] != seenConn[(L[5], L[0])][-2]:
                        log('True variant numbers differ', L[-2], seenConn[(L[5], L[0])][-2])
                        count += 1
                        if L[-2] < seenConn[(L[5], L[0])][-2]:
                            trueConn.append(L)
                        if L[-2] > seenConn[(L[5], L[0])][-2]:
                            trueConn.append(seenConn[(L[5], L[0])])
                        del seenConn[(L[5], L[0])]
                        continue
                    if L[-1] == 'RED' and seenConn[(L[5], L[0])][-1] == 'RED':
                        log(L, 'both red')
                    if L[-1] == 'DASH' and seenConn[(L[5], L[0])][-1] == 'DASH':
                        log(L, 'both dash')
                    if L[-1] == 'DASH' and seenConn[(L[5], L[0])][-1] == 'RED':
                        log(L, 'one dash one red')
                    if L[-1] == 'RED' and seenConn[(L[5], L[0])][-1] == 'DASH':
                        log(L, 'one dash one red')
                    if L[-3] > seenConn[(L[5], L[0])][-3]:
                        if L[-1] == 'GREEN' or L[-1] == 'DASH':
                            trueConn.append(L)
                        if L[-1] != 'GREEN' and seenConn[(L[5], L[0])][-1] == 'GREEN':
                            log(L, 'bad connection dominates.')
                    elif L[-3] < seenConn[(L[5], L[0])][-3]:
                        if seenConn[(L[5], L[0])][-1] == 'GREEN' or seenConn[(L[5], L[0])][-1] == 'DASH':
                            trueConn.append(seenConn[(L[5], L[0])])
                        if seenConn[(L[5], L[0])][-1] != 'GREEN' and L[-1] == 'GREEN':
                            log(L, 'bad connection dominates.')
                    del seenConn[(L[5], L[0])]
        for pair, L in seenConn.items():
            log(pair, 'happens only once', L)
            if L[-1] == 'GREEN' and L[-1] == 'DASH':
                trueConn.append(L)
        Llines = trueConn
    log(count, 'cases have diferent number of truevars')

    out = open(args.output, 'w')
    if args.ground_truth != None:
        log('Trying to remove wrong connection basing on provided ground truth', args.ground_truth)
        Llines = removeWrongEdge(Llines, args.ground_truth, out)

    draw2(args.fasta, Llines, out, False)
    #log('gfa output')
    #final_L = transitive_reduction(Llines)
    #log(final_L)
    #draw2(args.fasta, final_L, 'reducted.gfa', True)
    #exit()
    if args.partition_list != None:
        nxGraph = gfaToNX(Llines)
        log(nxGraph)
        subgraphs = clustering.divide(nxGraph, args.partition_list)
    #log('nx graph constructed')
    #log('start clustering')
    #subgraphs = clustering.divide(nxGraph)
    #print(len(subgraphs))

if __name__ == '__main__':
    main()
