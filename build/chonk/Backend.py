#!/usr/bin/env python3
import math, os, sys 
"""
backend functions
"""

def alignment_midpoint(Aln):
	# return the mapped midpoint position of an alignment
	return int(((Aln.reference_end - Aln.reference_start) / 2) + Aln.reference_start)

def check_path(f):
	# simply check if the file exists, if not exit
	if not os.path.isfile(f):
		sys.stderr.write('FATAL ERROR: {} not found\n'.format(f))
		sys.exit(1)

def check_chrom(contig,flag):
	# formats the contig name according to the format in the BAM file
	# returns a formatted contig name
	if 'chr' in contig and flag==False: contig=contig.replace('chr','')
	if 'chr' not in contig and flag==True: contig='chr'+contig
	return contig

def get_gc_bin(gc_content):
	# takes the gc content as a decimal as input (example 0.3454)
	# returns the bin number with respect to gc content

	# this conditional because of the header 
	try: gc_content = round(float(gc_content) * 100)
	except: return -9

	#gc bin encompasses the percentages (gc_bin*4-3)--(gc_bin*4). example: gc_bin 4 == 13%-16%
	gc_bin = 0
	if gc_content % 4 == 0 and gc_content != 0:
		gc_bin = int(gc_content / 4)
	else:
		gc_bin = int( (gc_content + (4-(gc_content % 4))) / 4 )
	return gc_bin

def reporter(s):
	space, cols = 4, 80
	buff = ' ' * space
	ln = '-' * cols+buff
	sys.stderr.write('\n'+ln+'\n'+buff+s+buff+'\n'+ln+'\n\n')

def tokenize_user_contigs(r):
	return tuple(r.split(','))

def tuple2string(t):
	return ','.join(map(str,(t)))

def avg(x):
	x_sum = 0
	x_count = 0
	for num in x:
		if num == None: continue
		x_sum += x
		x_count += 1
	return (x_sum / x_count)


class Welford(object):
	""" 
	    Implements Welford's algorithm for computing a running mean
	    and standard deviation as described at: 
		http://www.johndcook.com/standard_deviation.html
		https://gist.github.com/alexalemi/2151722
	"""

	def __init__(self,lst=None):
		self.k = 0
		self.M = 0
		self.S = 0
		
		self.__call__(lst)
	
	def update(self,x):
		if x is None: return
		self.k += 1
		newM = self.M + (x - self.M)*1./self.k
		newS = self.S + (x - self.M)*(x - newM)
		self.M, self.S = newM, newS

	def consume(self,lst):
		lst = iter(lst)
		for x in lst:
			self.update(x)
	
	def __call__(self,x):
		if hasattr(x,"__iter__"):
			self.consume(x)
		else:
			self.update(x)
			
	@property
	def mean(self):
		return self.M

	@property
	def std(self):
		if self.k==1:
			return 0
		return math.sqrt(self.S/(self.k-1))