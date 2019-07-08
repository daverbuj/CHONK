#!/usr/bin/env python3
import chonk.Backend as Backend
import os,pysam,sys

class Bam(object):
	def __init__(self,f):
		# check if file exists		
		Backend.check_path(f)
		self.bam_path = os.path.abspath(f)
		# bam attribute is a pysam AlignmentFile object
		self.bam = pysam.AlignmentFile(f,'r')
		
		self.chrom_flag = True # reference contigs begin with "chr"
		if 'chr' not in self.bam.references[0]: self.chrom_flag=False

		# get sample name
		if self.bam.header.get('RG') == None: 
			sys.stderr.write('FATAL ERROR: bam file is missing read group (@RG) entry in the header\n')
			sys.exit(1)

		self.sample = tuple(sorted(set( x['SM'] for x in self.bam.header['RG'] )))

		if len(self.sample) == 0: 
			sys.stderr.write('FATAL ERROR: bam file read group is missing sample name (@RG SM:sample_name) in the header\n')
			sys.exit(1)

		if len(self.sample) > 1: 
			sys.stderr.write('WARNING: bam file read group contains more than one sample name entry... Using {} as the sample name\n'.format(self.sample[0]))

		self.sample = self.sample[0]
