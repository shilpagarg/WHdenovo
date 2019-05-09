#! /usr/bin/env python3

import subprocess

def call_variants(fastqs, outputPath):
	# fastqs is a list of file names.
	fastqcmd = ' '.join(fastqs)
	subprocess.call('cat ' + fastqcmd + ' > %s/all.fasta' % outputPath, shell = True)
	subprocess.call('minimap2 -c --cs -x asm20 -DP --no-long-join -n500 %s/all.fasta %s/all.fasta > %s/all.paf' % (outputPath, outputPath, outputPath), shell = True)
	pass