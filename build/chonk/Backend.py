#!/usr/bin/env python3
import math, os, sys, subprocess as sp
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

def get_gc_perc(fasta, window_size, contig, start, end):
	gc = 0
	win_size = 0
	region = '{}:{}-{}'.format(contig, start, end)
	# get the sequence of the region selected
	fasta_cmd = 'samtools faidx {} {}'.format(fasta, region)
	sequence_raw = sp.check_output(fasta_cmd, shell = True)
	sequence_raw = sequence_raw.decode("utf-8")
	region_fasta = '>'+region
	sequence = sequence_raw.replace(region_fasta,'').strip()
	# count up the number of g's and c's in the sequence
	for base in sequence:
		win_size += 1
		if base == 'G' or base == 'C':
			gc += 1
	gc_perc = gc / window_size
	return gc_perc

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

# find mean of a list excluding Nones
def mean(lst):
	lst_sum = 0
	lst_count = 0
	for num in lst:
		if num == None: continue
		lst_sum += num
		lst_count += 1
	return (lst_sum / lst_count)

# find median of a list excluding Nones
def median(lst):
	newLst = []
	for num in lst:
		if num != None:
			newLst.append(num)
	sortedLst = sorted(newLst)
	lstLen = len(newLst)
	index = (lstLen - 1) // 2
	if (lstLen % 2):
		return sortedLst[index]
	else:
		return (sortedLst[index] + sortedLst[index + 1])/2.0

# deprecated since recursion was close to causing a stack overflow in some cases
'''def create_windows(window_size, start, end, windows = []):
	# returns a list of tuples with start and end values for each window in the region
	if (end - start) >= window_size:
		win_end = start + (window_size-1)
		win = (start, win_end)
		windows.append(win)
		return create_windows(window_size, win_end +1, end, windows)
	else:
		return windows'''

def create_windows(window_size, start, end):
    # 0-based I believe
    windows = []
    tmp_start = start
    tmp_end = start + (window_size-1)
    for x in range( int( (end-start)/window_size ) ):
        windows.append( (tmp_start, tmp_end) )
        tmp_start = tmp_end + 1
        tmp_end = tmp_start + (window_size-1)
    return windows

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