#!/usr/bin/env python3
from chonk.Bam import Bam
from chonk.Metadata import Metadata
import chonk.Backend as Backend
import chonk.Breakpoint as Breakpoint
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import chonk.SupportingFragments as SupportingFragments
import pysam, json, os, sys, csv, numpy as np


### Features needed
# * Coverage
# sv region DOC (/chrom doc) == sv_doc_fc       DONE
# left and right flank DOC (/left and right DOC) == left_flank_doc_fc , right_flank_doc_fc DONE
# sv region windowed gc (/null gc read count. and find the mean and STD)
# same as above but for the flanking regions


def features(Args):
	"""
	main function for SV feature extraction
	"""
	
	# metadata file
	Meta = Metadata()
	Meta.read_json(Args.m)

	# check if alignment file exists
	AlnFile = Bam(Meta.bam_path)

	# check that contigs selected match those in breakpoint file
	# if no contigs selected, use the contigs available in the breakpoint file
	contigs = Breakpoint.selected_contigs(Args.r, bp_contigs(Args.b), AlnFile)

	# output bed file
	output_file = Args.o + Meta.sample + '_' + '_'.join(contigs) + '_features.chonk'
	if contigs == Meta.contigs:
		output_file = Args.o + Meta.sample + '_features.chonk'
	bed_fh = open(output_file,'wt')
	bed_writer = csv.writer(bed_fh, delimiter='\t')
	bed_writer.writerow( ['#chrom', 'start', 'end', 'svtype'] ) # ADD THE OTHER FEATURES TO THIS

	# start extracting features from each sv in the breakpoint file
	sv_bed = open(Args.b)
	for sv in csv.reader(sv_bed, dialect='excel-tab'):
		# this is for validated SVs / training set
		try:
			contig, start, end, svtype, cistart, ciend, genome, genotype = sv[0], int(sv[1]), int(sv[2]), sv[3], sv[4], sv[5], sv[6], sv[7]
		# this is for unvalidated SVs
		except:
			contig, start, end, svtype, cistart, ciend = sv[0], int(sv[1]), int(sv[2]), sv[3], '.', '.', 

		## extract coverage features
		# sv region
		sv_doc_fc = doc_coverage(AlnFile, Meta, contig, start, end)
		sv_gc_mean, sv_gc_std = gcbin_coverage(AlnFile, Meta, contig, start, end)
		# left flank
		lf_doc_fc = doc_coverage(AlnFile, Meta, contig, start-1000, start)
		lf_gc_mean, lf_gc_std = gcbin_coverage(AlnFile, Meta, contig, start-1000, end)
		# right flank
		rf_doc_fc = doc_coverage(AlnFile, Meta, contig, end, end+1000)
		rf_gc_mean, rf_gc_std = gcbin_coverage(AlnFile, Meta, contig, end, end+1000)

		## extract supporting fragment features
		supporting_frag_features = SupportingFragments.supporting_frags(AlnFile, Meta, contig, start, end, cistart, ciend, svtype)


		


##################################################

def bp_contigs(bp_file):
	# find contigs used in breakpoint file
	bp_contigs = ()
	for row in csv.reader(bp_file, dialect="excel-tab"):
		if row[0] not in bp_contigs:
			bp_contigs += (row[0],)
	return bp_contigs

def doc_coverage(Alnfile, Meta, contig, start, end):
	meta_doc = Meta.doc[contig]

	# find the doc fold change for the sv region
	read_length_sum = 0
	region_size = end - start
	for Aln in AlnFile.bam.fetch(reference = contig, start = start, end = end):
		midpoint = Backend.alignment_midpoint(Aln)
		if start <= midpoint <= end:
			read_length_sum += Aln.reference_length
	if region_size != 0:
		doc = read_length_sum / region_size
	if meta_doc != 0:
		doc_fc = sv_doc / meta_doc

	return doc_fc

def gcbin_coverage(AlnFile, Meta, contig, start, end):
	fasta = Meta.fasta

	# choose window size based on the size of the sv
	window_lengths = Meta.window_lengths.sort()
	window_size = window_lengths[0]
	if (end - start) >= (2 * window_lengths[1]):
		window_size = window_lengths[1]

	# find the mean and STD binned coverage in the SV region in relation to GC content
	Region = Welford()

	read_count = 0
	sv_windows = create_windows(window_size, start, end)
	for window in sv_windows:
		win_start, win_end = window[0], window[1]

		# find the gc bin for the window
		gc_bin = Backend.get_gc_bin(get_gc_perc(fasta, window_size, contig, win_start, win_end))

		# count up the number of reads in the window
		for Aln in AlnFile.bam.fetch(reference = contig, start = win_start, end = win_end):
			midpoint = Backend.alignment_midpoint(Aln)
			if start <= midpoint <= end:
				read_count += 1

		# calculate the fold change of read count in window to the count of the null with similar gc content
		key = Backend.tuple2string( (contig, window_size, gc_bin) )
		gc_fc = read_count / Meta.gc_rc[key]

		# add the gc fold change value to list of all windows in the sv region
		Region.update(gc_fc)

	# return mean and STD for all the windows in the region
	return (Region.mean, Region.std)

def create_windows(window_size, start, end, windows = []):
	# returns a list of tuples with start and end values for each window in the region
	if (end - start) >= window_size:
		win_end = start + (window_size-1)
		win = (start, win_end)
		windows.append(win)
		return create_windows(window_size, win_end +1, end, windows)
	else:
		return windows

def get_gc_perc(fasta, window_size, contig, start, end):
	gc = 0
	region = '{}:{}-{}'.format(contig, start, end)
	# get the sequence of the region selected
	fasta_cmd = 'samtools faidx {} {}'.format(fasta, region)
	sequence_raw = sp.check_output(fa_cmd, shell = True)
	sequence_raw = sequence_raw.decode("utf-8")
	region_fasta = '>'+region
	sequence = sequence_raw.replace(region_fasta,'').strip()
	# count up the number of g's and c's in the sequence
	for base in sequence:
		if base == 'G' or base == 'C':
			gc += 1
	gc_perc = gc / window_size
	return gc_perc













