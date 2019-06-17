#!/usr/bin/env python3
from chonk.Alignment import Alignment
import operator, pysam

# class of split-aln

def get_split(Aln=None):
    # returns a list of secondary alignments
    s_alns=[]
    if Aln.has_tag("SA"):
        s_alns = Aln.get_tag("SA").split(";")
    return s_alns # len(s_alns) == number of secondary alignments


################################
# main split-read function
def splitread(Aln,s_alns,chrom):
	# list of Alignment objects which are each alignment for a read
	alns = []
	breaks=[]

	# Alignment attributes for the primary alignment 
	paln = Alignment()
	paln.primary_alignment(Aln)
	alns.append(paln)
	
	# Alignment attributes for each secondary alignment
	for sec_aln in s_alns:
		#sec_aln is a string
		if sec_aln == '': continue
		sa_list = sec_aln.split(',')
		if chrom != sa_list[0]: continue 
		saln = Alignment()
		saln.secondary_alignment(sa_list)
		alns.append(saln)

	if len(alns)<2: return breaks 
	
	# sort the split alignments by query position
	key = operator.attrgetter("qpos")
	alns.sort(key=key,reverse=False)
	
	# infer breakpoints by comparing the mapped positions
	# in pairs of alignments, starting with the first aln
	# with respect to the query position
	for i in range(0,len(alns)-1):
		left = alns[i] # "left" aln on query
		right = alns[i+1] # "right" aln on query

		# skip if not on same strand or chromosome or if mapped positions overlap
		if (left.chrom != right.chrom
			or  list(set(range(left.lpos,left.rpos+1)) & set(range(right.lpos,right.rpos+1)))
			): 
				continue
						
		if left.strand == right.strand:
		# checking if deletion or duplication:
			breakpoint_start, breakpoint_end = left.rpos, right.lpos
			svtype = 'DEL'
			if left.lpos > right.rpos:
				breakpoint_start, breakpoint_end = right.lpos, left.rpos
				svtype = 'DUP'

		elif left.strand != right.strand:
			svtype = 'INV'
			breakpoint_start, breakpoint_end = left.rpos, right.rpos
			if Aln.mate_is_reverse:
				breakpoint_start, breakpoint_end = left.lpos, right.lpos

		breaks.append( (chrom,breakpoint_start,breakpoint_end,svtype,'SR') )

	return breaks