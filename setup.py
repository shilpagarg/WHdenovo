#!/usr/bin/env python

import sys
import os.path
from setuptools import setup
import subprocess

if sys.version_info < (3, 4):
	sys.stdout.write("At least Python 3.4 is required.\n")
	sys.exit(1)

#install_whtrio = 'cd trioasm/whatshap_trioasm/ && python setup.py install'
#install_whtrio_p = subprocess.call(install_whtrio, shell = True)
#if install_whtrio_p != 0:
#	sys.exit(1)

#install_whind = 'cd trioasm/whatshap/ && python setup.py install'
#install_whind_p = subprocess.call(install_whind, shell = True)
#if install_whtrio_p != 0:
#	sys.exit(1)

setup(
	name = 'whdenovo',
	version = 1.0,
	author = 'WHdenovo authors',
	author_email = ['shilpa.garg2k7@gmail.com'], 
	url = 'https://github.com/shilpagarg/WHdenovo',
	description = 'A graph-based approach to diploid assembly for single samples and trios (In principle, should work for duos).',
	entry_points={'console_scripts': ['whdenovo = whdenovo.__main__:main']},
	install_requires = [
		'xopen',
		'networkx',
		'pystream-protobuf'
	],
	extras_require = {
		'dev': ['Cython', 'pytest', 'sphinx', 'sphinx_issues'],
	},
	python_requires = '>=3.4',
	packages = ['whdenovo'],
	classifiers = [
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Programming Language :: Cython",
		"Programming Language :: Python",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
