import sys
from collections import defaultdict
import argparse
from multiprocessing import Pool
from cigarStringParser2 import snpDetector
from Bio import SeqIO
import subprocess

def align(fasta, t):
    
    paf_run = subprocess.Popen(['minimap2', '-c', '--cs', '-x', 'asm20', 
                                '-DP', '--no-long-join', '-n500', '-t', str(t),
                                fasta, fasta], stderr = subprocess.PIPE, 
                               stdout = subprocess.PIPE)
    
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
    
    paflines = sorted(paflines)
    blocks = defaultdict(list)
    for line in paflines:
        if line == '':
            continue
        readname = line.split('\t')[5]
        blocks[readname].append(line)
    return blocks

def coverage(tgt_cover, pos):
    c = 0
    for i in range(len(tgt_cover)):
        if tgt_cover[i][0] < pos:
            c += tgt_cover[i][1]
        if tgt_cover[i][0] == pos:
            while tgt_cover[i][0] == pos and i < len(tgt_cover):
                c += 1
                i += 1
            break
        if tgt_cover[i][0] > pos:
            break
    return c + 1 # This 1 is for the target itself

def isHet(b1, b2):
    if b1 == 'a' and b2 == 't':
        return False
    elif b1 == 't' and b2 == 'a':
        return False
    elif b1 == 'c' and b2 == 'g':
        return False
    elif b1 == 'g' and b2 == 'c':
        return False
    else:
        return True

def log(*content):
    out = []
    for i in content:
        out.append(str(i))
    sys.stderr.write(' '.join(out)+'\n')

def variantCall(targetName, paflines):
    #if targetName != 'S1_1_child1':
    #    return None
    Pos_varCount = defaultdict(dict)
    Pos_Ref = dict()
    # keys are positions on target, 
    # values are dictionaries where keys are alleles values are their count
    query_record = dict()
    tgt_cover = []
    if targetName == 'S1_47_child2' or targetName == 'S1_309_child1':
        log(targetName,'block!')
        log(len(paflines))
    for line in paflines:
        tokens = line.split('\t')
        query_name = tokens[0]
        query_length = int(tokens[1])
        query_start = int(tokens[2])
        query_end = int(tokens[3])
        query_dir = tokens[4]
        target_name = tokens[5]
        target_length = int(tokens[6])
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
        #print(target_name, query_name, aln_block_length, len(snp_coord)/aln_block_length, snp_coord)
        if targetName == 'S1_33_child1' and query_name == 'S1_553_child2':
            log('wierd aln 553_2', snp_coord)
        if snp_coord != None:
            tgt_cover.append((target_start, 1))
            tgt_cover.append((target_end, -1))
        else:
            continue
        query_var = dict()
        # Keys are SNP position on the target, values are alleles the queries have.

        for pos_snp in snp_coord:
            try:
                Pos_varCount[pos_snp[0]][pos_snp[2]] += 1
            except KeyError:
                Pos_varCount[pos_snp[0]][pos_snp[2]] = 1
            query_var[pos_snp[0]] = pos_snp[2]
            Pos_Ref[pos_snp[0]]= pos_snp[1]
        query_record[query_name] = (query_var, target_start, target_end, 
                                    target_length, query_start, query_end, 
                                    query_length, query_dir, aln_block_length, cigar1)
    
    tgt_cover = sorted(tgt_cover)
    if targetName == 'S1_33_child1':
        log('wierd target ori vars', Pos_varCount)
    Pos_varCounttmp = defaultdict()
    for pos, Vars in Pos_varCount.items():
        pos_cov = coverage(tgt_cover, pos)
        ref_cov = pos_cov
        for b, n in Vars.items():
            ref_cov -= n
        try:
            Pos_varCount[pos][Pos_Ref[pos]] += ref_cov
        except KeyError:
            Pos_varCount[pos][Pos_Ref[pos]] = ref_cov
        alleles = list(Vars.keys())
        true_var = 0
        if targetName == 'S1_33_child2':
            log('before filter', Pos_varCount)
        #if pos == 413:
        #    log('wierd pos 413', pos_cov, Vars)
        for b, n in Vars.items():
            if n / pos_cov >= 0.4 and pos_cov >= 5:
                true_var += 1
        #if pos == 413:
        #    log('wierd pos 413', pos_cov, true_var, Vars)
        if true_var == 2:
            Pos_varCounttmp[pos] = Vars

        '''
        if len(alleles) == 2:
            if alleles[0] == 'a' and ( alleles[1] == 'c' or alleles[1] == 'g'):
                Pos_varCounttmp[pos] = Vars
            elif alleles[0] == 't' and ( alleles[1] == 'c' or alleles[1] == 'g'):
                Pos_varCounttmp[pos] = Vars
            elif alleles[0] == 'c' and ( alleles[1] == 'a' or alleles[1] == 't'):
                Pos_varCounttmp[pos] = Vars
            elif alleles[0] == 'g' and ( alleles[1] == 'a' or alleles[1] == 't'):
                Pos_varCounttmp[pos] = Vars
        '''
    Pos_varCount = Pos_varCounttmp.copy()
    #if targetName == 'S1_33_child1':
    #    log('wierd target 33_1 filtered var', Pos_varCount)
    Llines = []
    for query_name, query_info in query_record.items():
        # First check prefix-suffix.
        if query_info[8] < 2000:
            continue
        
        if query_info[1] == 0 and query_info[2] == query_info[3] \
           and query_info[4] == 0 and query_info[5] == query_info[6]:
            # Complete overlap
            L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        
        elif query_info[4] == 0 and query_info[5] < query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] == query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] < query_info[6] \
           and query_info[1] == 0 and query_info[2] == query_info[3] \
           and query_info[7] == '+':
            # target-suffix <==> prefix-query '+'
            L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]


        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] <= query_info[3] \
             and query_info[4] > 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] >= 0 and query_info[5] == query_info[6] \
             and query_info[7] == '+':
            # query-suffix '+' <==> prefix-target
            L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]

        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] == 0 and query_info[5] < query_info[6] \
             and query_info[7] == '-':
             # query-suffix '-' <==> prefix-target
             L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] == query_info[3] \
             and query_info[4] == 0 and query_info[5] < query_info[6] \
             and query_info[7] == '-':
             # query-suffix '-' <==> prefix-target
             L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]
        elif query_info[1] == 0 and query_info[2] < query_info[3] \
             and query_info[4] == 0 and query_info[5] == query_info[6] \
             and query_info[7] == '-':
             # query-suffix '-' <==> prefix-target
             L = [query_name, query_info[7], target_name, 
                 '+', str(query_info[9])]

        elif query_info[4] > 0 and query_info[5] == query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '-':
           # target-suffix  <==> prefix-query '-'
           L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        elif query_info[4] == 0 and query_info[5] == query_info[6] \
           and query_info[1] > 0 and query_info[2] == query_info[3] \
           and query_info[7] == '-':
           # target-suffix  <==> prefix-query '-'
           L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]
        elif query_info[4] > 0 and query_info[5] == query_info[6] \
           and query_info[1] == 0 and query_info[2] < query_info[3] \
           and query_info[7] == '-':
           # target-suffix  <==> prefix-query '-'
           L = [target_name, '+', query_name, 
                 query_info[7], str(query_info[9])]

        else:
            continue

        match_true, snp_true, nVar = 0, 0, 0
        query_var = query_info[0]

        if targetName == 'S1_33_child2':# and query_name == 'S1_61_child2':
            log('wierd target! S1_33_child2', Pos_varCount)
            #log('wierd query! S1_61_child2', query_info[0])
            
        for pos, Vars in Pos_varCount.items():
            if pos <= query_info[1] or pos >= query_info[2]:
                continue
            pos_cov = coverage(tgt_cover, pos)
            try:
                #if isHet(Pos_Ref[pos], query_var[pos]):
                #if Vars[Pos_Ref[pos]] / pos_cov >= 0.4 \
                #   and Vars[query_var[pos]] / pos_cov >= 0.4:
                qv = query_var[pos]
                snp_true += 1
                log(target_name, query_name, pos, Pos_Ref[pos], query_var[pos], Vars, 'RED')
                nVar += 1
                #else:
                    # Do something for different strand SNP?.
                #    log(target_name, query_name, pos, Pos_Ref[pos], query_var[pos], Vars, 'STRAND?')
                #    pass
                #
            except KeyError:
                #if Vars[Pos_Ref[pos]] / pos_cov > 0.4 \
                #   and Vars[Pos_Ref[pos]] / pos_cov < 0.6:
                match_true += 1
                log(target_name, query_name, pos, Pos_Ref[pos], 'Match', Vars, 'GREEN')
                nVar += 1
        if target_name == 'S1_40_child2' and query_name == 'S1_27_child1':
            log('caonima???', snp_true, match_true, nVar, Pos_varCount)
            
        if nVar == 0:
            log(target_name, query_name, match_true, snp_true, nVar, 'DASH')
            L.append('DASH')
        elif match_true / nVar > 0.5:
            log(target_name, query_name, match_true, snp_true, nVar, 'GREEN')
            L.append('GREEN')
        elif snp_true / nVar > 0.5:
            log(target_name, query_name, match_true, snp_true, nVar, 'RED')
            L.append('RED')
        else:
            log(target_name, query_name, match_true, snp_true, nVar, 'somethingWrong')
            L.append('ERROR')
        Llines.append(L)
    return Llines

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
    '''
    for read in consideredReads:
        if output != None:
            graph.write('L\t%s\t+\t+%s\t-\t10000M'%(read,read))
            
        else:
            sys.stdout.write('L\t%s\t+\t+%s\t-\t10000M'(read,read))
    '''
    for record in SeqIO.parse(fasta, 'fasta'):
        if str(record.id) not in consideredReads:
            continue
        if output != None:
            graph.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')
        else:
            sys.stdout.write('S\t' + str(record.id) + '\t' + str(record.seq) + '\n')

    if output != None:
        graph.close()


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--paf', metavar = 'PAF', type = str, 
                         help = 'PAF file, if not given, we will run minimap2 on the input fasta')
    parser.add_argument('-o', '--output', metavar = 'GFA', type = str, 
                         help = 'Output to FILE, default to stdout. ')
    parser.add_argument('-f', '--fasta', metavar = 'FASTA', type = str, required = True,
                        help = 'FASTA file which corresponds to the PAF file.')
    parser.add_argument('-t', '--threads', metavar = 'INT', type = int, default = 4, 
                         help = 'Maximum number of threads to use. [4]')
    args = parser.parse_args()

    if args.paf == None:
        paflines = align(args.fasta, args.threads)
    else:
        paflines = open(args.paf, 'r').read().split('\n')[:-1]

    blocks = split_jobs(paflines)
    
    p = Pool(args.threads)

    processes = []
    for target, paflines in blocks.items():
        processes.append(p.apply_async(variantCall, (target, paflines)))

    p.close()
    p.join()


    Llines = []
    for proc in processes:
        Llines_block = proc.get()
        Llines.extend(Llines_block)

    draw(args.fasta, Llines, args.output)

if __name__ == '__main__':
    main()
