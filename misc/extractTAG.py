import sys
total_count = 0
tag_count = 0
while True:
    read = sys.stdin.readline()
    read = read.strip()
    fields = read.split('\t')
    readname = fields[0]
    HP = ''
    PS = ''
    for i in fields:
        if i[:2] == 'HP':
            HP = i
        if i[:2] == 'PS':
            PS = i
    total_count += 1
    if HP != '' and PS != '':
        tag_count += 1
    sys.stdout.write(readname+'\t'+HP+'\t'+PS+'\n')
    if not read:
        break
sys.stderr.write('Total reads #: '+str(total_count)+'\n'+'Has both tag #: '+str(tag_count)+'\n')
