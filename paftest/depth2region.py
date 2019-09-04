import argparse

def readDepth(depth):
    with open(depth, 'r') as DF:
        for line in DF.readlines():
            line = line.rstrip()
            tokens = line.split('\t')
            chrm = tokens[0]
            pos = int(tokens[1])
            dpth = int(tokens[2])
            yield (chrm, pos, dpth)

def outputRegion(region, form, minl):
    if region[1] == -1 or region[2] == -1:
        return False
    if region[2] - region[1] < minl:
        return False

    if form == 'BED':
        print(region[0]+'\t'+str(region[1])+'\t'+str(region[2]))
    elif form == 'region':
        print(region[0]+':'+str(region[1])+'-'+str(region[2]))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth', type = str, metavar = 'FILE', help = 'Output file from "samtools depth".', required = True)
    parser.add_argument('-m', '--mincov', type = int, metavar = 'INT', default = 20, help = 'Minimum depth [20] to consider as high-coverage.')
    parser.add_argument('-min-length', type = int, metavar = 'INT', default = 10000, help = 'Minimum region length [10000] to consider')
    parser.add_argument('-w', '--window', type = int, default = 5, help = 'While going through each position, allow [5] bp that has coverage lower than threshold. (i.e. A short drop in coverage')
    parser.add_argument('-f', '--format', type = str, metavar = '[BED|region]', default = 'BED', help = 'Output format. BED refers to "chr\\tfrom\\tto" for each line; region refers to "chr:from-to" for each line')
    args = parser.parse_args()
    
    lastChr = ''
    lastPos = -1
    regionStart = -1
    regionEnd = -1
    lastEnd = -1
    inregion = False
    newChr = False
    for i in readDepth(args.depth):
        
        currentChr = i[0]
        currentPos = i[1]

        if currentChr != lastChr:
            newChr = True
            if lastChr != '' and inregion:
                outputRegion((lastChr, regionStart, lastPos), args.format, args.min_length)
        else:
            newChr = False
            
        if i[2] >= args.mincov:

            if not inregion:
                inregion = True
                if newChr:
                    regionStart = i[1]
                elif currentPos - lastEnd >= args.window:
                    outputRegion((currentChr, regionStart, lastEnd), args.format, args.min_length)
                    regionStart = i[1]
            else:
                if currentPos - lastPos >= args.window:
                    outputRegion((currentChr, regionStart, lastPos), args.format, args.min_length)
                    regionStart = i[1]
        else:
            if inregion:
                lastEnd = lastPos
            inregion = False

        lastPos = currentPos
        lastChr = currentChr

    if lastEnd < regionStart:
        outputRegion((currentChr, regionStart, currentPos), args.format, args.min_length)
    else:
        outputRegion((currentChr, regionStart, lastEnd), args.format, args.min_length)

if __name__ == '__main__':
    main()
