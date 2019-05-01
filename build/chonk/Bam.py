#!/usr/bin/env python3
import chonk.Backend as Backend
import pysam

class Bam():
	def __init__(self,f):
		Backend.check_path(f)
		self.bam = pysam.AlignmentFile(f,'r')
		
		self.chrom_flag = True # reference contigs begin with "chr"
		if 'chr' not in self.bam.references[0]: self.chrom_flag=False
