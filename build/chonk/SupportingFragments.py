#!/usr/bin/env python3
import pysam, os, csv, json
from chonk.Alignment import Alignment
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE


#########################

class Fragment():

	def __init__(self):
		self.supp = False

		# only one is true if supp
		self.split = False
		self.disc = False
		self.clip = False

		self.mapq = None
		self.baseq = []

		# left_clip / right_clip has to be in proper orientation
		self.left_clip = 0
		self.right_clip = 0 

		self.n_forward=0 
		self.n_reverse=0 
		self.n_aln = 0 
		self.qname = None

	def update(self, features): # features == (svtype, (L.strand, L.mapq), (R.strand, R.mapq), baseq, queryname)
		svtype, left, right, basequals, queryname = features[0], features[1], features[2], features[3], features[4]
		l_strand, l_mapq, r_strand, r_mapq = left[0], left[1], right[0], right[1]
		# Update / occupy the existing Fragment object with the appropriate feature information

		# Ensure that SR > DPE > CLIP > NonSU, and that nothing less than the existing svtype can overwrite feature data
		if self.split == True:
			if svtype == 'SR': print('SUPP FRAG, line 34 - SR read already found') # if so, fix, and also check line 51
			return
		if self.disc == True:
			if svtype == 'DPE': print('SUPP FRAG, line 37 - DPE read already found')
			if svtype == 'CLIP' or svtype == None:
				return
		if self.clip == True:
			if svtype == 'CLIP': print('SUPP FRAG, line 41 - CLIP read already found')
			if svtype == None:
				return

		# writes / overwrites to Fragment if SR, always
		if svtype == 'SR':
			self.qname = queryname
			self.supp, self.split, self.disc, self.clip = True, True, False, False
			self.left_clip, self.right_clip = 0, 0
			self.mapq = (l_mapq + r_mapq) / 2
			self.n_aln = 2
			for qual in basequals:
				self.baseq.append(qual)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# writes / overwrites to Fragment if DPE and the existing is clip or None
		elif svtype == 'DPE':
			self.qname = queryname
			self.supp, self.split, self.disc, self.clip = True, False, True, False
			self.left_clip, self.right_clip = 0,0
			self.mapq = (l_mapq + r_mapq) / 2
			self.n_aln = 2
			for qual in basequals:
				self.baseq.append(qual)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1

		# only writes / overwrites to Fragment if not supporting
		elif svtype == 'CLIP':
			self.qname = queryname
			self.supp, self.split, self.disc, self.clip = True, False, False, True
			self.left_clip, self.right_clip = 0,0
			self.mapq = (l_mapq + r_mapq) / 2
			self.n_aln = 2
			for qual in basequals:
				self.baseq.append(qual)
			if l_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1
			if r_strand == '+':
				self.n_forward += 1
			else:
				self.n_reverse += 1







	def update_baseq(self, qualities)
		for qual in qualities:
			self.baseq.apend(qual)







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
	if cistart == '.' or cistart == '0,0':
		ci_start = '{},{}'.format((-med_mpd)/2, med_mpd/2)
	if ci_end == '.' or ci_end == '0,0':
		ci_end = '{},{}'.format((-med_mpd)/2, med_mpd/2)
	l_window_start = int(start) + int(float(cistart.split(',')[0]))
	l_window_end = int(start) + int(float(cistart.split(',')[1]))
	r_window_start = int(end) + int(float(ciend.split(',')[0]))
	r_window_end = int(end) + int(float(ciend.split(',')[1]))

	#if there is overlap in the windows, then shorten the windows to prevent overlap
	if l_window_end > r_window_start:
		l_window_start = int(start) - (read_len / 2)
		l_window_end = int(start) + (read_len / 2)
		r_window_start = int(end) - (read_len / 2)
		r_window_end = int(end) + (read_len / 2)

	if l_window_end > r_window_start:
		print('SUPP FRAG, line 66 - there is still overlap in the windows') # if so, shorten window more and try again

	# extract SV supporting fragment/read features from both windows
	su_features = window_features(AlnFile, Meta, contig, l_window_start, l_window_end, r_window_start, r_window_end, svtype)







#################################


def window_features(AlnFile, Meta, contig, lwin_start, lwin_end, rwin_start, rwin_end, truesvtype):
	fragments = {}

	## look through all the reads in the left window. Checking for DPE, SR, or Non-SU
	for Aln in AlnFile.fetch(reference = contig, start = int(lwin_start), end = int(lwin_end)):
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
			fragments[Aln.queryname].update(dpe_features)
			continue

		# SU - clipped reads
		clip_features = clipread(Aln, truesvtype)
		if clip_features:
			fragment[Aln.queryname].update(clip_features)
			continue

		# NonSU
		else:




		









	# SU - discordant paired end just left window
		if not sr_found:
			if dpe(Aln,meta,contig,rwin_start, rwin_end, truesvtype):
				svcaller = 'DPE'
				dpe_found, qname, svtype, l_read, r_read = dpe(Aln, meta, contig, rwin_start, rwin_end, truesvtype)
				l_strand, l_mapq = l_read[0], l_read[1]
				r_strand, r_mapq = r_read[0], r_read[1]
				if not qname in reads_dict['SU']:
					reads_dict['SU'][qname] = [[l_strand, l_mapq, svtype, svcaller], [r_strand, r_mapq, svtype, svcaller]]
				else:
					print('DPE problem')
	# nonSU - left window
		if not sr_found and not dpe_found: # then this Aln is not a supporting read
			qname, mapq = Aln.query_name, Aln.mapping_quality
			strand = '+'
			if Aln.is_reverse: strand = '-'
			svcaller = None
			svtype = None
			if not qname in reads_dict['nonSU']:
				reads_dict['nonSU'][qname] = [[strand, mapq, svtype, svcaller]]
			else:
				reads_dict['nonSU'][qname].append([strand, mapq, svtype, svcaller]) # this should happen very rarely i believe
	# nonSU - right window excluding SUs found in left window
	for Aln in AlnFile.fetch(reference = contig, start = int(rwin_start), end = int(rwin_end)):
		qname = Aln.query_name
		if not qname in reads_dict['SU']:
			mapq = Aln.mapping_quality
			strand = '+'
			if Aln.is_reverse: strand = '-'
			svcaller = None
			svtype = None
			if not qname in reads_dict['nonSU']:
				reads_dict['nonSU'][qname] = [[strand, mapq, svtype, svcaller]]
			else:
				reads_dict['nonSU'][qname].append([strand, mapq, svtype, svcaller]) # this should happen very rarely i believe

	# Calculate total number of SU and nonSU frags and reads in both CI windows of each SV
	su_frag_count = 0
	nonsu_frag_count = 0
	su_read_count = 0
	nonsu_read_count = 0
	su_pos_strand = 0
	su_neg_strand = 0
	nonsu_pos_strand = 0
	nonsu_neg_strand = 0
	del_count = 0
	dup_count = 0
	su_mapq_sum = 0
	nonsu_mapq_sum = 0
	svt = 'None'
	for frag in reads_dict['SU']:
		su_frag_count += 1
		for read in reads_dict['SU'][frag]:
			su_read_count += 1
			su_mapq_sum += read[1]
			if read[2] == 'DEL':
				del_count+=1
			elif read[2] == 'DUP':
				dup_count+=1
			if read[0] == '+':
				su_pos_strand += 1
			elif read[0] == '-':
				su_neg_strand += 1				
	for frag in reads_dict['nonSU']:
		nonsu_frag_count += 1
		for read in reads_dict['nonSU'][frag]:
			nonsu_read_count += 1
			nonsu_mapq_sum += read[1]
			if read[0] == '+':
				nonsu_pos_strand += 1
			elif read[0] == '-':
				nonsu_neg_strand += 1
	su_mapq_mean, nonsu_mapq_mean = None, None
	if su_mapq_sum > 0 and su_read_count > 0:
		su_mapq_mean = su_mapq_sum/su_read_count
	if nonsu_mapq_sum > 0 and nonsu_read_count > 0:
		nonsu_mapq_mean = nonsu_mapq_sum/nonsu_read_count
	if del_count > dup_count:
		svt = 'DEL'
	if del_count < dup_count:
		svt = 'DUP'
	features = (su_frag_count, nonsu_frag_count, su_read_count, nonsu_read_count, truesvtype, svt, su_pos_strand, su_neg_strand, nonsu_pos_strand, nonsu_neg_strand,
		su_mapq_mean, nonsu_mapq_mean)
	return features


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
		saln = Alignment()
		saln.secondary_alignment(sa_list)
		
		# check to see that Right alignment is correct
		if ((saln.rpos > owins and saln.lpos < owine) # is in the right window
			and Left.chrom == Right.chrom # on the same chromosome
			and not list(set(range(Left.lpos,Left.rpos+1)) & set(range(Right.lpos,Right.rpos+1)))
			):
			Right = saln

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
				r_mapq = Aln.mapping_quality # should this be zero or same as left mapq? 

			# if this split read detects the same svtype as the one in the breakpoint file, this is a supporting split-read fragment
			if svtype == truesvtype:
				features = ('DPE', (l_strand l_mapq), (r_strand, r_mapq), Aln.query_alignment_qualities, Aln.query_name)

	return features

################################

def clipread(Aln, truesvtype, left_window = True):
	features = ()

	## check if read is clipped
	left_clip, right_clip = clipped(Aln)
	if not left_clip and not right_clip:
		return features

	# check if clip is informative for a DEL
	if truesvtype == 'DEL':
		if (left_window and right_clip)
			or (not left_window and left_clip
			):
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), None, Aln.query_alignment_qualities, Aln.query_name)

	# check if clip is informative for a DUP
	elif truesvtype == 'DUP':
		if (left_window and left_clip)
			or (not left_window and right_clip
			):
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), None, Aln.query_alignment_qualities, Aln.query_name)

	#check if clip is informative for INV
	elif truesvtype == 'INV':
		if left_clip or right_clip:
			strand = '+'
			if Aln.is_reverse: strand = '-'
			features = ('CLIP', (strand, Aln.mapping_quality), None, Aln.query_alignment_qualities, Aln.query_name)	

	return features


def clipped(Aln):
	left_clip = False
	right_clip = False

	# if the Aln is not considered a split read
	if not Aln.has_tag("SA"):
		# check if read is left clipped
		cigar = Aln.cigartuples
		if ((cigar[0][0] == 4 or cigar[0][0] == 5)
			and cigar[0][1] > 20
			):
			left_clip = True
		# check if read is right clipped
		if ((cigar[-1][0] == 4 or cigar[-1][0] == 5)
			and cigar[-1][1] > 20
			):
			right_clip = True

	return left_clip, right_clip

	