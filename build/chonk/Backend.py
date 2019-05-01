#!/usr/bin/env python3
import sys, os
"""
backend functions
"""

def check_path(f):
	if not os.path.isfile(f):
		sys.stderr.write('FATAL ERROR: {} not found\n'.format(f))
		sys.exit(1)

def check_chrom(chrom,flag):
	if 'chr' in chrom and flag==False: chrom=chrom.replace('chr','')
	if 'chr' not in chrom and flag==True: chrom='chr'+chrom
	return chrom