#!/usr/bin/env python

import sys
import os.path
from setuptools import setup
import subprocess
import versioneer
from distutils.command.build_ext import build_ext as _build_ext


cmdclass = versioneer.get_cmdclass()

class build_ext(cmdclass.get('build_ext', _build_ext)):
	def run(self):
		print('Start to compile whatshap_trio')
		install_whtrio = 'cd trioasm/whatshap_trioasm/ && python setup.py build_ext -i'
		install_whtrio_p = subprocess.call(install_whtrio, shell = True)
		if install_whtrio_p != 0:
			print('whatshap_trio compile failed')
			sys.exit(1)

if sys.version_info < (3, 4):
        print("At least Python 3.4 is required.")
        sys.exit(1)
cmdclass['build_ext'] = build_ext

setup(
	name = 'whdenovo',
	version = 1.0,
	author = 'WHdenovo authors',
	author_email = ['shilpa.garg2k7@gmail.com'], 
	url = 'https://github.com/shilpagarg/WHdenovo',
	description = 'A graph-based approach to diploid assembly for single samples and trios (In principle, should work for duos).',
	entry_points={'console_scripts': ['whdenovo = whdenovo.__main__:main']},
	cmdclass = cmdclass,
	install_requires = [
		'xopen',
		'networkx',
		'pystream-protobuf',
		'Biopython'
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
