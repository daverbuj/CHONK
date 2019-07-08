#!/usr/bin/env python3
import re

### functions that use cigarstrings to 
#   determine the query start position
#   and the right most aligned position
def get_cigartuples(cigar=None):
	# reverse-engineered pysam.cigartuples() 
	CODE2CIGAR= ['M','I','D','N','S','H','P','=','X','B']
	CIGAR2CODE = dict([ord(y), x] for x, y in enumerate(CODE2CIGAR))
	CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
	parts = CIGAR_REGEX.findall(cigar)
	return [(CIGAR2CODE[ord(y)], int(x)) for x,y in parts]

def get_qpos(cigartuple=None):
	# returns the first position on the query
	# the first aligned portion relative to the query
	# 1-base
	qInd=0
	qStart=0
	for (flg,leng) in cigartuple:
		if flg == 0 or flg == 1 or flg==4 or flg==7 or flg==8: qInd+=leng
		# this is the first mapped postion
		if (flg==0 or flg==7) and qStart==0: 
			qStart=qInd-leng+1
			break
	return qStart

def get_rpos(cigartuple=None,lpos=None):
	# return the right mapped position relative to the reference
	rInd=0
	for (flg,leng) in cigartuple:
		if flg == 0 or flg == 2 or flg==3 or flg==7 or flg==8: rInd+=leng
	rpos=lpos + rInd - 1
	return rpos

##################

class Alignment():
	"""
	aligned segment. A split-read has at least 2 aligned segments
	"""
	def __init__(self):
		self.chrom=None
		self.lpos=None # left alignment position (rel. reference)
		self.rpos=None # right alignment position (rel. reference)
		self.strand='+' # forward by default
		self.mapq=None 
		self.qpos=None # first mapped position (rel. query)
	
	def primary_alignment(self,Aln=None):
		# populate the Alignment attributes for the 
		# primary alignment, takes full advantage
		# of pysam functions.

		self.chrom = Aln.reference_name
		self.lpos = Aln.reference_start #Aln.reference_start is in 0-base
		self.rpos = Aln.reference_end
		if Aln.is_reverse: self.strand = '-'
		self.mapq = Aln.mapping_quality
		self.qpos = Aln.query_alignment_start+1 # 1-base position
		
	def secondary_alignment(self,sa_list=None):
		# populate the Alignment attributes for the
		# secondary alignment lpos is 0-based
	
		self.chrom = sa_list[0]
		self.lpos = int(sa_list[1])-1
		self.strand = sa_list[2]
		self.mapq = int(sa_list[4])

		cigartuple = get_cigartuples(sa_list[3])
		# determine the right aligned position
		self.rpos = get_rpos(cigartuple,self.lpos)
		# determine the first aligned position on query
		self.qpos = get_qpos(cigartuple)


