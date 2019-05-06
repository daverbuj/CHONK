import sys
import pysam
import re

# class of split-aln
class Alignment():
	def __init__(self):
		self.chrom=None
		self.lpos=None # left alignment position (rel. reference)
		self.rpos=None # right alignment position (rel. reference)
		self.strand='+' # forward by default
		self.mapq=None 
		self.qpos=None # first mapped position (rel. query)
		#debugging
		self.cigar=None

def get_split(Aln=None):
    # returns a list of secondary alignments
    s_alns=[]
    if Aln.has_tag("SA"):
        s_alns = Aln.get_tag("SA").split(";")
    return s_alns # len(s_alns) == number of secondary alignments

def primary_aln(Aln=None,Alignment=None):
	# populate the Alignment attributes for the 
	# primary alignment, takes full advantage
	# of pysam functions.

	Alignment.chrom = Aln.reference_name
	Alignment.lpos = Aln.reference_start+1 #Aln.reference_start is in 0-base
	Alignment.rpos = Aln.reference_end
	if Aln.is_reverse: Alignment.strand = '-'
	Alignment.mapq = Aln.mapping_quality
	Alignment.qpos = Aln.query_alignment_start+1 # 1-base position
	
	#debugging
	Alignment.cigar= Aln.cigarstring  
	return Alignment

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

def secondary_alignment(sa_list=None,Alignment=None):
	# populate the Alignment attributes for the
	# secondary alignment lpos is 1-based

	#['1', '1829833', '-', '20S30M50S', '60', '0;']
    #[CHROM, LEFT-POS, STRAND, CIGAR, MAPQ, ? ] 
	
	Alignment.chrom = sa_list[0]
	Alignment.lpos = int(sa_list[1])
	Alignment.strand = sa_list[2]
	Alignment.mapq = int(sa_list[4])

	#debugging
	Alignment.cigar=sa_list[3]

	cigartuple = get_cigartuples(sa_list[3])
	# determine the right aligned position
	Alignment.rpos = get_rpos(cigartuple,Alignment.lpos)
	# determine the first alignmed position on query
	Alignment.qpos = get_qpos(cigartuple)
	return Alignment

def sort_qpos(sub_li): 
    l = len(sub_li) 
    for i in range(0, l): 
        for j in range(0, l-i-1): 
            if (sub_li[j][3] > sub_li[j + 1][3]): 
                tempo = sub_li[j] 
                sub_li[j]= sub_li[j + 1] 
                sub_li[j + 1]= tempo 
    return

###################################################################

bam_fh = sys.argv[1]
chrom = sys.argv[2]

bam = pysam.AlignmentFile(bam_fh,'r')

processed_reads = [] #list of reads that we looked at already

for Aln in bam.fetch(reference=chrom):
	# process split-reads
	s_alns = get_split(Aln)
	if len(s_alns)>0:
		# these reads have splits

		# get the read name
		read='1'
		if Aln.is_read2: read='2' # is_read2 should be called "is_second_in_pair"
		read_name = Aln.query_name+"-"+read

		# skip if we looked at the read before
		if read_name in processed_reads: continue

		processed_reads.append(read_name)

		alns = [] # list of Alignment objects which are each alignment for a read

		# load the Alignment attributes for the primary alignment 
		alns.append( primary_aln(Aln, Alignment()) )

		# load the Alignment attributes for each secondary alignment
		for sec_aln in s_alns:
			#sec_aln is a string
			if sec_aln == '': continue
			sa_list = sec_aln.split(',')
			if chrom != sa_list[0]: continue 

			alns.append( secondary_alignment(sa_list, Alignment()) )

		if len(alns)<2: continue 
		# test

		#print('-'*40)
		#for alignment in alns:
		#	print(read_name,alignment.lpos,alignment.rpos,alignment.qpos)
		#print('-'*40)

		# Adding all alignments to a list
		nl_alns = []
		for alignment in alns:
			nl_alns.append([read_name,alignment.lpos,alignment.rpos,alignment.qpos, 
							alignment.strand,alignment.chrom])
		# sorting read lists by qpos
		sort_qpos(nl_alns)
		
		for i in range(0,len(nl_alns)-1):
			left = nl_alns[i]
			right = nl_alns[i+1]
			# left and right are on the same strand
			if left[4] == right[4]:
				# left and right are on the same chr
				if left[5] == right[5]:
					# positions do not overlap
					left_range = set(range(left[1],left[2]+1))
					right_range = set(range(right[1],right[2]+1))
					if not list(left_range & right_range):
						left_l = left[1]
						left_r = left[2]
						right_l = right[1]
						right_r = right[2]
						# checking if deletion or duplication:
						if left_r < right_l:
							sv = 'deletion'
							breakpoint_start = left_r
							breakpoint_end = right_l
						elif left_l > right_r:
							sv = 'duplication'
							breakpoint_start = right_l
							breakpoint_end = left_r
						#print(read_name, sv, chrom+':'+str(breakpoint_start)+'-'+str(breakpoint_end), 
						#	alignment.mapq) #which mapq is showing?

						# Printing for .bed file
						print(chrom,'\t',breakpoint_start,'\t',breakpoint_end,'\t', sv,'\t',read_name, sep='')
