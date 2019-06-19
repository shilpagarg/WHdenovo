import sys
from collections import defaultdict, OrderedDict
import argparse
from multiprocessing import Pool
from cigarStringParser2 import snpDetector
from Bio import SeqIO
import networkx as nx
import subprocess


class pileup():
    """
    Store the variant call from all query reads for one target read. And 
    apply filtering to get the true variants.
    All the coordinates included are coordinates on the target sequence.
    """
    def __init__(self, target_length, avgCov, targetName):
        self.targetName = targetName
        self.target_length = target_length
        self.avgCov = avgCov
        self.read_var = defaultdict(dict)
        self.read_order = []
        self.mode = 'normal'
        self.tgt_cover = dict()
        self.posToVarread = defaultdict(set)
        self.posRef = dict()
        self.pos_var_count = defaultdict(dict)
        self.maxCov = 0
        self.copyNum = 1
        self.trueHapReads = set()
        self.badPos = set()
        self.previous_aln = 0

    def add_snp(self, query_name, target_start, target_end, snp_coord):
        assert target_start >= self.previous_aln, 'Input paf should be sorted by target overlapping position'
        for snp in snp_coord:
            self.read_var[query_name][snp[0]] = snp[2]
            self.posToVarread[snp[0]].add(query_name)
            if snp[0] in self.posRef:
                assert self.posRef[snp[0]] == snp[1]
            self.posRef[snp[0]] = snp[1]
        self.read_order.append(query_name)
        self.tgt_cover[query_name] = [(target_start, 1), (target_end, -1)]
        self.previous_aln = target_start

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
        return c + 1

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

    def addRefAllele(self):
        for pos, Vars in self.pos_var_count.items():
            pos_cov = self.depth(pos)
            ref_cov = pos_cov
            for b, n in Vars.items():
                ref_cov -= n
            self.pos_var_count[pos][self.posRef[pos]] = ref_cov
        
    def removeError(self):
        tmp = defaultdict(dict)
        for pos, var_count in self.pos_var_count.items():
            Bad = 0
            for b, n in var_count.items():
                if n / self.depth(pos) <= 0.2 or n / self.depth(pos) >= 0.8:# or n < 5 or self.depth(pos) - n < 5:
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

    def cluster(self, T = 3):
        patterns = []
        clusters = []
        for read in self.read_order:
            p = pattern(self.read_var[read], self.tgt_cover[read][0][0], self.tgt_cover[read][1][0])
            if not patterns:
                patterns.append(p)
                clusters.append([read])
            else:
                picker = defaultdict(list)
                log('checking', read)
                for i in range(len(patterns)):
                    d = patterns[i].distance(p)
                    if d < T:
                        picker[d].append(i)
                if picker:
                    mind = min(picker)
                    idx = picker[mind]
                    clusters[idx[0]].append(read)
                    patterns[idx[0]].update(p)
                else:
                    patterns.append(p)
                    clusters.append([read])
            #log(read, clusters)
        picker = []
        for i in range(len(patterns)):
            if len(clusters[i]) <= 5:
                continue
            nVar = len(patterns[i].pos_var)
            picker.append((nVar, i))
        picker = sorted(picker)
        if len(picker) < 2:
            log('WARNING: LESS THEN 2 HAPLOTYPES.')
            return [], []
        hap1 = clusters[picker[0][1]]
        hap2 = clusters[picker[1][1]]

        log(self.targetName, hap1, hap2, 'nPatters', len(picker))

        return hap1, hap2

    def getTrueHaps(self):
        self.getAltVar()
        self.removeError()
        hap1, hap2 = self.cluster()

        return hap1, hap2

class pattern():
    '''
    for representing the het pattern of one haplotype
    '''
    def __init__(self, pos_var, start, end):
        assert type(pos_var) == dict
        self.pos_var = pos_var
        self.start = start
        self.end = end
    
    def update(self, p):
        self.start = min(self.start, p.start)
        self.end = max(self.end, p.end)
        for pos in p.pos_var:
            if pos in self.pos_var:
                if self.pos_var[pos] != p.pos_var[pos]:
                    del self.pos_var[pos]
            else:
                self.pos_var[pos] = p.pos_var[pos]
    
    def distance(self, p):
        log('start distance')
        log('self', self.pos_var)
        log('query', p.pos_var)
        ov1 = {}
        d = 0
        for pos, var in self.pos_var.items():
            if pos >= p.start and pos <= p.end:
                ov1[pos] = var
        log('ov1', ov1)
        ov1_pos = sorted(list(ov1.keys()))
        ov2 = {}
        for pos, var in p.pos_var.items():
            if pos >= self.start and pos <= self.end:
                ov2[pos] = var
        ov2_pos = sorted(list(ov2.keys()))
        log('ov2', ov2)
        i = 0
        while ov2_pos and i < len(ov1_pos):
            pos = ov1_pos[i]
            if pos == ov2_pos[0]:
                if ov2[pos] != ov1[pos]:
                    d += 1
                ov2_pos.pop(0)
                i += 1
            elif pos > ov2_pos[0]:
                ov2_pos.pop(0)
                d += 1
            elif pos < ov2_pos[0]:
                i += 1
                d += 1
        d += len(ov2_pos) + len(ov1_pos[i:])
        log('distance', d, '\n')
        
        return d

def align(fasta, t):
    cmd = 'minimap2 -c -x asm20 --cs -DP --no-long-join -n500 -t %d %s %s | sort -k8n -k6'%(t, fasta, fasta)
    sys.stderr.write('minimap2 cmd: %s\n'%cmd)
    sys.stderr.write('Aligning...\n')
    paf_run = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    #paf_run = subprocess.Popen(['minimap2', '-c', '--cs', '-x', 'asm20', 
    #                            '-DP', '--no-long-join', '-n500', '-t', str(t),
    #                            fasta, fasta, '|', 'sort', '-k8n', '-k6'], shell = True, stderr = subprocess.PIPE, 
    #                           stdout = subprocess.PIPE)
    com = paf_run.communicate()
    returnCode = paf_run.returncode
    stdout = com[0].decode()
    stderr = com[1].decode()
    if returnCode != 0:
        sys.stderr.write('Error: minimap2 reuturn code: %d. Stderr follows: \n' % returnCode)
        sys.stderr.write(stderr)
    # minimap2 -c --cs -x asm20 -DP --no-long-join -n500 48841_mhc_50k_x/pacbioccs_child.fasta 48841_mhc_50k_x/pacbioccs_child.fasta > 48841_mhc_50k_x/all.paf
    pafout = open(fasta.split('/')[-1]+'.paf', 'w')
    pafout.write(stdout)
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

def variantCall(targetName, paflines, avgCov):
    #Pos_varCount = defaultdict(dict)
    #Pos_Ref = dict()
    # keys are positions on target, 
    # values are dictionaries where keys are alleles values are their count
    query_record = dict()
    target_length = int(paflines[0].split('\t')[6])
    caller = pileup(target_length, avgCov, targetName)
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
        cigar1 = tokens[21][5:]
        cigar2 = tokens[22][5:]
        target_coord = target_start
        if query_dir == '+':
            query_coord = 0
        else:
            query_coord = query_end
        snp_coord = snpDetector(cigar2, target_coord, query_coord, query_dir)
        if query_name == 'S1_271_child1':
            log('DEBUG S1_271_child1', snp_coord)
        caller.add_snp(query_name, target_start, target_end, snp_coord)
        query_record[query_name] = [0, target_start, target_end, 
                                    target_length, query_start, query_end, 
                                    query_length, query_dir, aln_block_length, cigar1]
    hap1, hap2 = caller.getTrueHaps()

    Llines = []
    for query_name, query_info in query_record.items():
        # First check prefix-suffix.
        if query_info[8] < 4000:
            continue
        '''
        if query_info[1] == 0 and query_info[2] == query_info[3] \
           and query_info[4] == 0 and query_info[5] == query_info[6]:
            # Complete overlap
            L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        '''
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

        if query_name in hap1:
            L.append('GREEN')
        elif query_name in hap2:
            L.append('RED')
        else:
            continue

        '''
        match_true, snp_true, nVar = 0, 0, 0
        query_var = query_info[0]

        #if targetName == 'S1_33_child2':# and query_name == 'S1_61_child2':
        #    log('wierd target! S1_33_child2', Pos_varCount)
        #    #log('wierd query! S1_61_child2', query_info[0])
            
        for pos, Vars in Pos_varCount.items():
            if pos <= query_info[1] or pos >= query_info[2]:
                continue
            try:
                qv = query_var[pos]
                snp_true += 1
                #if target_name == 'S1_4053_child2':
                #    log(target_name, query_name, pos, Pos_Ref[pos], query_var[pos], Vars, 'RED')
                nVar += 1
            except KeyError:
                match_true += 1
                #if target_name == 'S1_4053_child2':
                #    log(target_name, query_name, pos, Pos_Ref[pos], 'Match', Vars, 'GREEN')
                nVar += 1
            
        if nVar == 0:
            #log(target_name, query_name, match_true, snp_true, nVar, 'DASH')
            L.append('DASH')
        elif match_true / nVar == 1:
            #log(target_name, query_name, match_true, snp_true, nVar, 'GREEN')
            L.append('GREEN')
        elif snp_true / nVar > 0:
            #log(target_name, query_name, match_true, snp_true, nVar, 'RED')
            L.append('RED')
        else:
            #log(target_name, query_name, match_true, snp_true, nVar, 'somethingWrong')
            L.append('ERROR')
        '''
        Llines.append(L)
    return Llines

def gfaToNX(Llines):
    g = nx.Graph()
    for line in Llines:
        if line[-1] != 'GREEN':
            continue
        n1 = line[0]
        n2 = line[2]
        g.add_node(n1)
        g.add_node(n2)
        g.add_edge(n1, n2)

    return g

def draw2(fasta, Llines, output, haveDoneTR):

    if output != None:
        graph = open(output, 'w')

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
            content = 'L\t' + '\t'.join([Lline[0], Lline[4], Lline[5], Lline[9], str(Lline[10])+'M', Lline[12]]) + '\n'
        
        if output != None:
            graph.write(content)
        else:
            sys.stdout.write(content)
    for record in SeqIO.parse(fasta, 'fasta'):
        if str(record.id) not in consideredReads:
            continue
        if output != None:
            graph.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
        else:
            sys.stdout.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
    if output != None:
        graph.close()

def draw(fasta, Llines, output):

    if output != None:
        graph = open(output, 'w')

    consideredReads = set()
    existingPair = set()
    for Lline in Llines:

        consideredReads.add(Lline[0])
        consideredReads.add(Lline[2])

        content = 'L\t' + '\t'.join(Lline) + '\n'
        
        if output != None:
            graph.write(content)
        else:
            sys.stdout.write(content)

    for record in SeqIO.parse(fasta, 'fasta'):
        if str(record.id) not in consideredReads:
            continue
        if output != None:
            graph.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
        else:
            sys.stdout.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')

    if output != None:
        graph.close()

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

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--paf', metavar = 'PAF', type = str, 
                         help = 'PAF file, should be consistent with the input FASTA file; if not given, we will run minimap2 on the input fasta and generate one.')
    parser.add_argument('-o', '--output', metavar = 'GFA', type = str, 
                         help = 'Output to FILE, default to stdout. ')
    parser.add_argument('fasta', metavar = 'FASTA', type = str,
                        help = 'Input all read records in one FASTA file.')
    parser.add_argument('-t', '--threads', metavar = 'INT', type = int, default = 4, 
                         help = 'Maximum number of threads to use. [4]')
    parser.add_argument('--het', metavar = 'FLOAT', type = float, default = 0.1,
                         help = 'Heterozygosity rate, set as a filter for abnormal alignment. [0.1]')
    parser.add_argument('-c', '--avg-coverage', metavar = 'INT', type = int, default = 40, 
                         help = 'average coverage, set as a filter for abnormal alignment. [40]')
    args = parser.parse_args()

    if args.paf == None:
        paflines = align(args.fasta, args.threads)
    else:
        paflines = open(args.paf, 'r').read().split('\n')[:-1]

    blocks = split_jobs(paflines)
    
    p = Pool(args.threads)

    processes = []
    for target, paflines in blocks.items():
        processes.append(p.apply_async(variantCall, (target, paflines, args.avg_coverage)))

    p.close()
    p.join()

    Llines = []
    for proc in processes:
        Llines_block = proc.get()
        Llines.extend(Llines_block)
    draw2(args.fasta, Llines, args.output, False)
    #log('gfa output')
    final_L = transitive_reduction(Llines)
    log(final_L)
    draw2(args.fasta, final_L, 'reducted.gfa', True)
    exit()
    nxGraph = gfaToNX(Llines)
    log('nx graph constructed')
    log('start clustering')
    subgraphs = clustering.divide(nxGraph)
    print(len(subgraphs))

if __name__ == '__main__':
    main()
