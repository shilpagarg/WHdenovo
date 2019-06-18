#! /usr/bin/env python3
import sys, os
import re
import subprocess
from subprocess import PIPE

'''
This code is for validating test runs on simulated data
Usage: python validator.py outputPath/bc1/ pacbio/child.fasta
'''

path = os.path.abspath(sys.argv[1])
pacbioFA = sys.argv[2]
fileList = os.listdir(path)
readslist = []
allreadslist = []
for i in fileList:
    elif re.match(r'.*_block_*.allreads', i):
        allreadslist.append(i)
subprocess.call('cat %s/*.allreads > %s/all.allreads'%(path, path), shell = True)
subprocess.call('python src/validate_simulate.py %s/all.allreads'%path, shell = True)
subprocess.call("cat %s/all.allreads | cut -d' ' -f1 | sort -u | uniq > %s/all.tmpreads"%(path, path), shell = True)
subprocess.call('python src/get_unpartitioned.py %s/all.tmpreads %s'%(path, pacbioFA), shell = True)
