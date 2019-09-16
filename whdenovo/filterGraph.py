from collections import defaultdict
import sys
def flt(GFAin, GFAout):
    degrees = defaultdict(set)
    allNodes = set()
    gfa = open(GFAin, 'r').read().split('\n')[:-1]
    for line in gfa:
        if line[0] == 'S':
            allNodes.add(line.split('\t')[1])
        if line[0] == 'L':
            tokens = line.split('\t')
            n1 = tokens[1]
            n2 = tokens[3]
            degrees[n1].add(n2)
            degrees[n2].add(n1)

    wrong = set()
    for n in allNodes:
        if len(degrees[n]) > 70 or len(degrees[n]) == 0:
            #print(n, len(degrees[n]))
            wrong.add(n)
    o = open(GFAout, 'w')
    for line in gfa:
        if line[0] == 'S':
            if line.split('\t')[1] in wrong:
                continue
        if line[0] == 'L':
            tokens = line.split('\t')
            n1 = tokens[1]
            n2 = tokens[3]
            if n1 in wrong or n2 in wrong:
                continue
        o.write(line + '\n')

def main():
    flt(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
