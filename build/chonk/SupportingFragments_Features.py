#!/usr/bin/env python3
import pysam, numpy as np
from chonk.Alignment import Alignment
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import chonk.Backend as Backend


#########################

class Fragment():

	def __init__(self):
		self.supp = False

		# only one is true if supp
		self.split = False
		self.disc = False
		self.clip = False

		self.mapq = []
		self.baseq = []
		self.mapq_mean = 0
		self.baseq_mean = 0

		self.n_forward=0 
		self.n_reverse=0 
		self.n_aln = 0
		self.qname = None

	def update(self, features): # features == (su_caller, (l_strand, l_mapq), (r_strand, r_mapq), basequals, queryname)
		su_caller, left, right, basequals, queryname = features[0], features[1], features[2], features[3], features[4]
		l_strand, l_mapq, r_strand, r_mapq = left[0], left[1], right[0], right[1]
		# Update / occupy the existing Fragment object with the appropriate feature information

		# Ensure that SR > DPE > CLIP > NonSU, and that nothing less than the existing su_caller can overwrite feature data
		if self.split == True and su_caller != 'SR':
			return
		if self.disc == True:
			if su_caller == 'DPE': print('SUPP FRAG, line 37 - DPE read already found')
			if su_caller == 'CLIP' or su_caller == None:
				return
		if self.clip == True:
			if su_caller == 'CLIP': print('SUPP FRAG, line 41 - CLIP read already found')
			if su_caller == None:
				return

		# appends to Fragment if another SR found with same query name
		if su_caller == 'SR' and self.split == True:
			self.mapq.extend([l_mapq, r_mapq])
			self.baseq.extend(basequals)
			self.n_aln += 2
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# writes / overwrites to Fragment.
		elif su_caller == 'SR':
			self.supp, self.split, self.disc, self.clip = True, True, False, False
			self.n_forward, self.n_reverse = 0,0
			self.baseq = []
			self.mapq = [l_mapq, r_mapq]
			self.n_aln = 2
			self.qname = queryname
			self.baseq.extend(basequals)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# writes / overwrites to Fragment if DPE and the existing is clip or None
		elif su_caller == 'DPE':
			self.supp, self.split, self.disc, self.clip = True, False, True, False
			self.n_forward, self.n_reverse = 0,0
			self.baseq = []
			self.mapq = [l_mapq, r_mapq]
			self.n_aln = 2
			self.qname = queryname
			self.baseq.extend(basequals)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# only writes / overwrites to Fragment if blank or currently consumed by NonSU
		elif su_caller == 'CLIP':
			self.supp, self.split, self.disc, self.clip = True, False, False, True
			self.n_forward, self.n_reverse = 0,0
			self.baseq = []
			self.mapq = [l_mapq]
			self.n_aln = 1
			self.qname = queryname
			self.baseq.extend(basequals)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# NonSU calls
		elif su_caller == None:
			self.baseq.extend(basequals)
			self.qname = queryname
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			self.n_aln += 1

	def update_baseq(self, qualities):
		self.baseq.extend(qualities)

#################################

def supporting_frags(AlnFile, Meta, contig, start, end, cistart, ciend, svtype):
	"""
	main function for SV supporting fragment feature extraction
	"""

	# import pertinent metadata
	tlen_mean = Meta.tlen[contig]
	tlen_std = Meta.tlen_std[contig]
	read_len = Meta.read_len[contig]

	# determine left and right CI windows
	l_window_start, l_window_end, r_window_start, r_window_end = create_flanking_windows(start, end, cistart, ciend, tlen_mean, read_len)

	# extract SV supporting fragment/read features from both windows
	fragment_data = window_features(AlnFile, Meta, contig, l_window_start, l_window_end, r_window_start, r_window_end, svtype)

	# analyze the fragments to determine features for each sv
	fragment_features = analyze_fragments(fragment_data)

	return fragment_features

#################################

def create_flanking_windows(start, end, ci_start, ci_end, tlen_mean, read_len):
	# determine left and right CI windows
	if ci_start == '.' or ci_start == '0,0':
		ci_start = '{},{}'.format((-tlen_mean)/2, tlen_mean/2)
	if ci_end == '.' or ci_end == '0,0':
		ci_end = '{},{}'.format((-tlen_mean)/2, tlen_mean/2)
	l_window_start = int(start) + int(float(ci_start.split(',')[0]))
	l_window_end = int(start) + int(float(ci_start.split(',')[1]))
	r_window_start = int(end) + int(float(ci_end.split(',')[0]))
	r_window_end = int(end) + int(float(ci_end.split(',')[1]))

	# if there is overlap in the windows, then shorten the windows to prevent overlap
	if l_window_end > r_window_start:
		l_window_start = int(start) - (read_len / 2)
		l_window_end = int(start) + (read_len / 2)
		r_window_start = int(end) - (read_len / 2)
		r_window_end = int(end) + (read_len / 2)

	# if there is still overlap, shorten windows more
	if l_window_end > r_window_start:
		l_window_start = int(start) - 20
		l_window_end = int(start) + 20
		r_window_start = int(end) - 20
		r_window_end = int(end) + 20
	return (l_window_start, l_window_end, r_window_start, r_window_end)

#################################

def window_features(AlnFile, Meta, contig, lwin_start, lwin_end, rwin_start, rwin_end, truesvtype):
	fragments = {}

	## look through all the reads in the left window. Checking for DPE, SR, CLIP, or Non-SU
	for Aln in AlnFile.bam.fetch(reference = contig, start = int(lwin_start), end = int(lwin_end)):
		# initialize each read as a Fragment object. Will end up being either supporting (SR, DPE, or clipped), or non-supporting
		Frag = Fragment()
		fragments[Aln.query_name] = Frag

		# SU - split-read
		s_alns = Splitread.get_split(Aln)
		sr_features = splitread(Aln, s_alns, contig, rwin_start, rwin_end, truesvtype)
		if sr_features:
			fragments[Aln.query_name].update(sr_features)
			continue

		# SU - discordant paired end	
		dpe_features = dpe(Aln, Meta, contig, rwin_start, rwin_end, truesvtype)
		if dpe_features:
			fragments[Aln.query_name].update(dpe_features)
			continue

		# SU - clipped reads
		left_clip_features = clipread(Aln, truesvtype)
		if left_clip_features:
			fragments[Aln.query_name].update(left_clip_features)
			continue

		# NonSU
		else:
			strand = '+'
			if Aln.is_reverse: strand = '-'
			left_nonsu_features = (None, (strand, None), (None, None), Aln.query_alignment_qualities, Aln.query_name)
			fragments[Aln.query_name].update(left_nonsu_features)

	## look through all the reads in the right window. Update if clips found, update base qualities, and add nonSUs
	for Aln in AlnFile.bam.fetch(reference = contig, start = int(rwin_start), end = int(rwin_end)):

		# if fragment hasnt been looked at then it could be a nonSU or a clip. check for both and update
		if not Aln.query_name in fragments:
			Frag = Fragment()
			fragments[Aln.query_name] = Frag

			# check if clipped read
			right_clip_features = clipread(Aln, truesvtype)
			if right_clip_features:
				fragments[Aln.query_name].update(right_clip_features)
				continue

			# check if NonSU
			else:
				strand = '+'
				if Aln.is_reverse: strand = '-'
				right_nonsu_features = (None, (strand, None), (None, None), Aln.query_alignment_qualities, Aln.query_name)
				fragments[Aln.query_name].update(right_nonsu_features)

		# if fragment has already been looked at, then update the Fragment to include basequality of Aln in the right window
		else:
			fragments[Aln.query_name].update_baseq(Aln.query_alignment_qualities)

	return fragments

##############

def splitread(Aln, s_alns, chrom, owins, owine, truesvtype):
	# opposite window start (owins), opposite window end (owine)
	features = ()

	# Left is aligment in left window.
	Left = Alignment()
	Left.primary_alignment(Aln)

	# look through alignments to determine the Right, which whill be in the right window.
	Right = Alignment()
	for sec_aln in s_alns:
		#sec_aln is a string
		if sec_aln == '': continue
		sa_list = sec_aln.split(',')
		if chrom != sa_list[0]: continue 
		Saln = Alignment()
		Saln.secondary_alignment(sa_list)
		
		# check to see which secondary alignment is the Right alignment
		if ((Saln.rpos > owins and Saln.lpos < owine) # is in the right window
			and not list(set(range(Left.lpos,Left.rpos+1)) & set(range(Saln.lpos,Saln.rpos+1))) # Left and Right dont overlap
			):
			Right = Saln

			# decide what type of SV this split is detecting
			if Left.strand == Right.strand:
			# checking if deletion or duplication:
				svtype = 'DEL'
				if Left.qpos > Right.qpos: # Danny: remember how we thought that the qpos might be buggy and was affecting the SR INV's? That wouldn't affect this since were only looking at left window, right?
					svtype = 'DUP'

			elif Left.strand != Right.strand:
				svtype = 'INV'

			# if this split read detects the same svtype as the one in the breakpoint file, this is a supporting split-read fragment
			if svtype == truesvtype:
				features = ('SR', (Left.strand, Left.mapq), (Right.strand, Right.mapq), Aln.query_alignment_qualities, Aln.query_name)

	return features

#################

def dpe(Aln, Meta, chrom, owins, owine, truesvtype):
	features = ()

	# get relevant metadata info for dpe calling
	tlen_mean = Meta.tlen[chrom]
	tlen_std = Meta.tlen_std[chrom]

	# cases in which dpe returns as empty.
	if Aln.mate_is_unmapped: return features
	try: 
		Aln.reference_end < Aln.next_reference_start
	except: return features

	# checking that Aln is proper for dpe
	if (Aln.is_paired # Read is Paired
		and Aln.next_reference_name == Aln.reference_name # Same chr
		and Aln.next_reference_start < owine # Mate is in the opposite window
		and (Aln.next_reference_start + Aln.query_length) > owins
		and Aln.reference_end < Aln.next_reference_start # Mates don't overlap
		):

		# checking if dpe
		dist = abs(Aln.template_length)
		if dist > (float(tlen_mean) + (3.5 * float(tlen_std))):

			# check if DEL or DUP (opposite strands)
			if Aln.is_reverse != Aln.mate_is_reverse:
				svtype = 'DEL'
				l_strand, r_strand = '+', '-'

				# checking if DUP
				if Aln.is_reverse:
					svtype = 'DUP'
					l_strand, r_strand = '-', '+'

			# check if INV (same strand)
			elif Aln.is_reverse == Aln.mate_is_reverse:
				svtype = 'INV'
				l_strand, r_strand = '+', '+'
				if Aln.is_reverse:
					l_strand, r_strand = '-', '-'

			# mapping quality of left and right alignment
			l_mapq = Aln.mapping_quality
			try:
				r_mapq = Aln.get_tag('MQ')
			except:
				r_mapq = 0 # should this be zero or same as left mapq? When no tag for mate's mapq. Tried investigating why this exception would occur

			# if this split read detects the same svtype as the one in the breakpoint file, this is a supporting split-read fragment
			if svtype == truesvtype:
				features = ('DPE', (l_strand, l_mapq), (r_strand, r_mapq), Aln.query_alignment_qualities, Aln.query_name)

	return features

################################

def clipread(Aln, truesvtype, left_window = True):
	features = ()

	# if Aln is split, don't consider it clipped
	if Aln.has_tag('SA'): return features

	## check if read is clipped
	left_clip, right_clip = clipped(Aln)
	if not left_clip and not right_clip:
		return features

	# check if clip is informative for a DEL
	if truesvtype == 'DEL':
		if ((left_window and right_clip)
			or (not left_window and left_clip)
			):
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), (None, None), Aln.query_alignment_qualities, Aln.query_name)

	# check if clip is informative for a DUP
	elif truesvtype == 'DUP':
		if ((left_window and left_clip)
			or (not left_window and right_clip)
			):
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), (None, None), Aln.query_alignment_qualities, Aln.query_name)

	#check if clip is informative for INV
	elif truesvtype == 'INV':
		if left_clip or right_clip:
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), (None, None), Aln.query_alignment_qualities, Aln.query_name)

	return features

def clipped(Aln):
	left_clip = False
	right_clip = False

	# if the Aln is not considered a split read
	if not Aln.has_tag("SA"):
		# check if read is left clipped
		cigar = Aln.cigartuples

		# some cigar strings == None, which causes problems
		try:
			if ((cigar[0][0] == 4 or cigar[0][0] == 5)
				and cigar[0][1] > 20
				):
				left_clip = True
			# check if read is right clipped
			if ((cigar[-1][0] == 4 or cigar[-1][0] == 5)
				and cigar[-1][1] > 20
				):
				right_clip = True
		except: return left_clip, right_clip

	return left_clip, right_clip

##################################################

def analyze_fragments(fragments):
	## initialize feature-relevant variables for analysis
	# sf_ratio
	n_sf = 0
	n_nonsf = 0
	# split_ratio
	n_split = 0
	n_nonsplit = 0
	# disc_ratio
	n_disc = 0
	n_nondisc = 0
	# clip_ratio
	n_clip = 0
	n_nonclip = 0
	# mapq
	sf_mapq = []
	nonsf_mapq = []
	#baseq
	sf_baseq = []
	nonsf_baseq = []


	# compile relevant information 
	for frag_name in fragments:
		Frag = fragments[frag_name]

		if Frag.supp: 
			# count of supporting fragments
			n_sf += 1
			# compile list of all mapq and baseq for SFs
			sf_mapq.extend(Frag.mapq)
			sf_baseq.extend(Frag.baseq)

		else: 
			# count of non supporting fragments
			n_nonsf += 1
			# compile list of all mapq and baseq for nonSFs
			nonsf_mapq.extend(Frag.mapq)
			nonsf_baseq.extend(Frag.baseq)

		# fraction of split reads to non split reads (split_ratio)
		if Frag.split: n_split += 1
		else: n_nonsplit += 1

		# fraction of DPE to non-SF (disc_ratio)
		if Frag.disc: n_disc += 1
		else: n_nondisc += 1

		# fraction of reads with clips (clip_ratio)
		if Frag.clip: n_clip += 1
		else: n_nonclip += 1

	# determine ratios
	sf_ratio = n_sf / n_nonsf
	split_ratio = n_split / n_nonsplit
	disc_ratio = n_disc / n_nondisc
	clip_ratio = n_clip / n_nonclip

	# determine mean and median mapqs and baseqs
	sf_mapq_mean = Backend.mean(sf_mapq)
	sf_mapq_median = Backend.median(sf_mapq)
	nonsf_baseq_mean = Backend.mean(sf_baseq)
	nonsf_baseq_median = Backend.median(sf_baseq)

	# return all supporting fragment related features
	return (sf_ratio, split_ratio, disc_ratio, clip_ratio, sf_mapq_mean, sf_mapq_median, nonsf_baseq_mean, nonsf_baseq_median)


	