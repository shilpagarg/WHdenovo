#import sys

def reverse_complement(b):
    d = {'a':'t', 't':'a', 'g':'c', 'c':'g'}
    if type(b) == str:
        return d[b]
    if type(b) == dict:
        tmp = dict()
        for base, n in b.items():
            tmp[d[base]] = n
        return tmp

def snpDetector(cigar, target_coord, query_coord, query_dir):
    snp_coord = []
    cigar = list(cigar)
    #target_coord -= 1
    matchRegion = [] # [(target0, target1, query0, query1)]
    while cigar != []:
        if cigar[0] == ':':
            # Matches
            match = []
            cigar.pop(0)
            for i in cigar:
                if i.isdigit():
                    match.append(i)
                else:
                    break
            matchLength = int(''.join(match))
            target_prev = target_coord
            query_prev = query_coord
            target_coord += matchLength
            if query_dir == '+':
                query_coord += matchLength
                matchRegion.append((target_prev, target_coord, query_prev, query_coord))
                #print(target_prev, target_coord, query_prev, query_coord)
            else:
                query_coord -= matchLength
                matchRegion.append((target_prev, target_coord, query_coord, query_prev))
                #print(target_prev, target_coord, query_coord, query_prev)
            for i in range(len(match)):
                cigar.pop(0)

        elif cigar[0] == '+':
            # Insertion in query
            queryInsertion = []
            cigar.pop(0)
            for i in cigar:
                if i.isalpha():
                    queryInsertion.append(i)
                else:
                    break
            query_prev = query_coord
            #if len(queryInsertion) > 100:
            #    return None
            if query_dir == '+':
                query_coord += len(queryInsertion)
            else:
                query_coord -= len(queryInsertion)
            #print(target_coord, '-', ''.join(queryInsertion), query_prev)
            for i in queryInsertion:
                cigar.pop(0)
        elif cigar[0] == '-':
            # Deletion in query
            queryDeletion = []
            cigar.pop(0)
            for i in cigar:
                if i.isalpha():
                    queryDeletion.append(i)
                else:
                    break
            #if len(queryDeletion) > 100:
            #    return None
            #print(target_coord, ''.join(queryDeletion), '-', query_coord)
            target_coord += len(queryDeletion)
            
            for i in queryDeletion:
                cigar.pop(0)
        elif cigar[0] == '*':
            # SNP !!!
            snp = []
            cigar.pop(0)
            for i in cigar:
                if i.isalpha():
                    snp.append(i)
                else:
                    break
            
            target_prev = target_coord
            query_prev = query_coord
            target_coord += 1
            if query_dir == '+':
                query_coord += 1
            else:
                query_coord -= 1
            #print(target_prev, snp[0], snp[1], query_prev)
            if query_dir == '+':
                snp_coord.append((target_prev, snp[0], snp[1]))
            else:
                snp_coord.append((target_prev, snp[0], snp[1]))
            for i in snp:
                cigar.pop(0)
        else:
            break
    return snp_coord
'''
paf = open(sys.argv[1], 'r').readlines()

for pafLine in paf:
    pafLine = pafLine.strip()
    tokens = pafLine.split('\t')
    query_name = tokens[0]
    query_length = int(tokens[1])
    query_start = int(tokens[2])
    query_end = int(tokens[3])
    query_dir = tokens[4]
    target_name = tokens[5]
    target_length = int(tokens[6])
    target_start = int(tokens[7])
    target_end = int(tokens[8])
    target_dir = '+'
    cigar = tokens[22][5:]
    if query_dir == '+':
        snp_coord = snpDetector(cigar, target_start, query_start, query_dir)
    else:
        snp_coord = snpDetector(cigar, target_start, query_end, query_dir)
    print(snp_coord)

'''
