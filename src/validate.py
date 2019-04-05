#! /usr/bin/env python3
import sys, os
import re
import subprocess
from subprocess import PIPE

'''
Usage: python validator.py temp_PID_TIME/bc1 pacbio/fasta
'''

path = os.path.abspath(sys.argv[1])
pacbioFA = sys.argv[2]
fileList = os.listdir(path)
readslist = []
allreadslist = []
for i in fileList:
    if re.match(r'.*_contig.*.trans.reads', i):
        readslist.append(i)
    elif re.match(r'.*_contig.*.trans.allreads', i):
        allreadslist.append(i)
subprocess.call('cat %s/*.reads > %s/all.reads'%(path, path), shell = True)
subprocess.call('cat %s/*.allreads > %s/all.allreads'%(path, path), shell = True)


print('=== For .reads ===')
subprocess.call('python src/compute_read_partitioning_test.py %s/all.reads'%path, shell = True)
print('=== For .allreads ===')
subprocess.call('python src/compute_read_partitioning_test_allreads.py %s/all.allreads'%path, shell = True)

subprocess.call("cat %s/all.allreads | cut -d' ' -f1 | sort -u | uniq > %s/all.tmpreads"%(path, path), shell = True)
subprocess.call('python src/get_unpartitioned.py %s/all.tmpreads %s'%(path, pacbioFA), shell = True)
