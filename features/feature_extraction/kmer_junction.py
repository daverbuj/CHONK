import csv, pysam, sys, os
import subprocess as sp
import pandas as pd
from chonk.Metadata import Metadata
from itertools import islice
from collections import defaultdict

def k_mers(sequence, k):
	it = iter(sequence)
	result = tuple(islice(it, k))
	if len(result) == k:
		yield "".join(result)
	for elem in it:
		result = result[1:] + (elem,)
		yield "".join(result)

def get_ref_alt_kmers(seq,mod_rlen,svtype, k):

	# isolate the sequences for just the refs
	l_ref_seq = seq[:mod_rlen*2]
	r_ref_seq = seq[-mod_rlen*2:]

	# create kmers for the refs
	l_ref_kmer, r_ref_kmer = set(), set()
	for lk in k_mers(l_ref_seq, k):
		l_ref_kmer.add(lk)
	for rk in k_mers(r_ref_seq, k):
		r_ref_kmer.add(rk)

	## create alt sequence and extract kmers based on svtype
	# if DEL
	if svtype == 'DEL':
		# isolate the sequences for the alt
		alt_seq = seq[:mod_rlen] + seq[-mod_rlen:]
		# create the kmers for the alt
		alt_kmer = set()
		for ak in k_mers(alt_seq, k):
			alt_kmer.add(ak)

		ref_alt_kmer = (l_ref_kmer,r_ref_kmer,alt_kmer)

		return ref_alt_kmer

	# if DUP
	elif svtype == 'DUP':
		# isolate the sequences for the alt
		alt_seq = seq[-mod_rlen*2:-mod_rlen] + seq[mod_rlen:mod_rlen*2]

		# create the kmers for the alt
		alt_kmer = set()
		for ak in k_mers(alt_seq, k):
			alt_kmer.add(ak)

		ref_alt_kmer = (l_ref_kmer, r_ref_kmer, alt_kmer)

		return ref_alt_kmer

	# if INV, have to account for the second sv junction
	elif svtype == 'INV':
		# isolate the sequence for the alts
		l_alt_seq = seq[:mod_rlen] + reverse_complement(seq[-mod_rlen*2:-mod_rlen])
		r_alt_seq = reverse_complement(seq[-mod_rlen:-mod_rlen*2]) + seq[-mod_rlen:]

		# create the kmers for the alts
		l_alt_kmer, r_alt_kmer = set(), set()
		for la in k_mers(l_alt_seq, k):
			l_alt_kmer.add(la)
		for ra in k_mers(r_alt_seq, k):
			r_alt_kmer.add(ra)

		ref_alt_kmer = (l_ref_kmer, r_ref_kmer, l_alt_kmer, r_alt_kmer)

		return ref_alt_kmer

def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	bases = list(seq)
	bases = [complement[base] for base in bases]
	n_seq = ''.join(bases)
	return n_seq[::-1]

def process_aln(AlnFile, chrom, start, end, K):
	kmers = set()
	for Aln in AlnFile.bam.fetch(chrom,start,end):
		# skip duplicates
		if Aln.is_duplicate and Aln.is_unmapped: continue
		if 'N' in Aln.query_sequence: continue
		if Aln.cigarstring==None: continue
		for k in k_mers(Aln.query_sequence,K):
			kmers.add(k)
	return kmers

def kp_process_aln(AlnFile, chrom, start, end, k, left_flank, svtype):
	"""
	this returns a set of kmers from left and right soft clipped regions
   	right now it's not correctly using the orientation of the clips (i.e. for deletions the left flank should use right clips) 
	Please fix this so that left flank -> right clip / right flank -> left clip (dels) ... etc 
	"""
	kmers = set()
	for Aln in AlnFile.bam.fetch(chrom,start,end):
		if Aln.cigarstring==None: continue
		if 'S' not in Aln.cigarstring: continue
		
		# if left side is soft clipped
		if Aln.cigartuples[0][0]==4:
			if (svtype=='DEL' and (not left_flank)) or (svtype=='DUP' and left_flank) or (svtype=='INV'):
				for lk in k_mers(Aln.query_sequence[0:Aln.query_alignment_start], k):
					kmers.add(lk)
		# if right side is soft clipped
		if Aln.cigartuples[-1][0]==4:
			if (svtype=='DEL' and left_flank) or (svtype=='DUP' and (not left_flank)) or (svtype=='INV'):
				for rk in k_mers(Aln.query_sequence[Aln.query_alignment_end:Aln.query_length], k):
					kmers.add(rk)

	return kmers

def get_ref_kmers(l_seq, r_seq, k, svtype):
	if svtype == 'INV':
		lo_seq, ro_seq = reverse_complement(l_seq), reverse_complement(r_seq)

	start_ref_kmers = set()
	end_ref_kmers = set()

	for sk in k_mers(l_seq, k):
		start_ref_kmers.add(sk)
	for ek in k_mers(r_seq, k):
		end_ref_kmers.add(ek)

	return start_ref_kmers, end_ref_kmers

def count_kmers(query,start_ref,end_ref):
	left,right=0,0
	
	left = len(query & start_ref)
	right = len(query & end_ref)

	return left, right 

def kp_ratios(start_kmers, end_kmers, suk_left, suk_right, euk_left, euk_right):
	start_ratio_left, start_ratio_right = 0.0, 0.0
	end_ratio_left, end_ratio_right = 0.0, 0.0

	if len(start_kmers)>0:
		start_ratio_left = suk_left / float(len(start_kmers))
		start_ratio_right = suk_right / float(len(start_kmers))
	if len(end_kmers)>0:
		end_ratio_left = euk_left / float(len(end_kmers))
		end_ratio_right = euk_right / float(len(end_kmers))

	return start_ratio_left, start_ratio_right, end_ratio_left, end_ratio_right

def compare_kmers(l_read_kmers, r_read_kmers, ref_alt):
	# initializing features we will end up writing out
	lwin_lref, lwin_rref, lwin_alt, rwin_lref, rwin_rref, rwin_alt  = 0,0,0,0,0,0
	lwin_lalt, lwin_ralt, rwin_lalt, rwin_ralt = 0,0,0,0

	# getting the left and right reference kmers
	try:
		lref_kmer, rref_kmer = ref_alt[0], ref_alt[1]
	except TypeError:
		print(ref_alt)

	if len(l_read_kmers) != 0:
		lwin_lref = len( lref_kmer & l_read_kmers ) / len(l_read_kmers)
		lwin_rref = len( rref_kmer & l_read_kmers ) / len(l_read_kmers)

	if len(r_read_kmers) != 0:
		rwin_lref = len( lref_kmer & r_read_kmers ) / len(r_read_kmers)
		rwin_rref = len( rref_kmer & r_read_kmers ) / len(r_read_kmers)

	# for DELS and DUPS
	if len(ref_alt) == 3:
		alt_kmer = ref_alt[2]
		if len(l_read_kmers) != 0:
			lwin_alt = len( alt_kmer & l_read_kmers ) / len(l_read_kmers)
		if len(r_read_kmers) != 0:
			rwin_alt = len( alt_kmer & r_read_kmers ) / len(r_read_kmers)

		return (lwin_lref, lwin_rref, lwin_alt, rwin_lref, rwin_rref, rwin_alt)

	# for INVs
	if len(ref_alt) == 4:
		lalt_kmer, ralt_kmer = ref_alt[2], ref_alt[3]
		if len(l_read_kmers) != 0:
			lwin_lalt = len( lalt_kmer & l_read_kmers ) / len(l_read_kmers)
			lwin_ralt = len( ralt_kmer & l_read_kmers ) / len(l_read_kmers)
		if len(r_read_kmers) != 0:
			rwin_lalt = len( lalt_kmer & r_read_kmers ) / len(r_read_kmers)
			rwin_ralt = len( ralt_kmer & r_read_kmers ) / len(r_read_kmers)

		return (lwin_lref, lwin_rref, lwin_lalt, lwin_ralt, rwin_lref, rwin_rref, rwin_lalt, rwin_ralt)

##########################################################################################################
# MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN #
##########################################################################################################

def kmer_features(chrom, start, end, AlnFile, sequence, svtype, mod_rlen, k, kp_k, ci_start, ci_end, lo_seq, ro_seq):

	### KMER JUNCTION (KJ)
	# get the reference and alt kmers from the sequence
	ref_alt_kmers = get_ref_alt_kmers(sequence, mod_rlen, svtype, k)

	# get the kmers from the reads in the bam file for the start and end regions of the SV
	l_start, l_end, r_start, r_end = start-mod_rlen, start, end, end+mod_rlen
	l_read_kmers = process_aln(AlnFile, str(chrom), l_start, l_end, k)
	r_read_kmers = process_aln(AlnFile, str(chrom), r_start, r_end, k)

	# determine ratio of matched kmers to all kmers for reads in each windows to left / right ref and alt
	kj_features = compare_kmers(l_read_kmers, r_read_kmers, ref_alt_kmers)

	### KMER PSEUDO-ALIGNMENT (KP)
	# START/LEFT alignment kmers for the soft clip, and side of clip
	left_flank = True
	cistart_start = int(ci_start[:ci_start.find(',')])
	cistart_end = int(ci_start[ci_start.find(',')+1:])
	start_kmers = kp_process_aln(AlnFile, chrom, int(start)+cistart_start, int(start)+cistart_end, kp_k, left_flank, svtype)
	# END/RIGHT alignment kmers for the soft clip, and side of clip
	left_flank = False
	ciend_start = int(ci_end[:ci_end.find(',')])
	ciend_end = int(ci_end[ci_end.find(',')+1:])
	end_kmers = kp_process_aln(AlnFile, chrom,int(end)+ciend_start ,int(end)+ciend_end, kp_k, left_flank,svtype)

	# get kmers for left and right reference flanks/ci_regions
	start_ref_kmers, end_ref_kmers = get_ref_kmers(lo_seq, ro_seq, kp_k, svtype)

	# find k-mers in soft-clipped sequence on one side on the reference on the other side
	suk_left, suk_right = count_kmers(start_kmers,start_ref_kmers,end_ref_kmers)
	euk_left, euk_right = count_kmers(end_kmers,start_ref_kmers,end_ref_kmers)

	kp_features = kp_ratios(start_kmers, end_kmers, suk_left, suk_right, euk_left, euk_right)

	return kj_features + kp_features