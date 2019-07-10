#!/usr/bin/env python3
from chonk.Bam import Bam
from chonk.Metadata import Metadata
from chonk.Backend import alignment_midpoint, Welford, create_windows, get_gc_bin, get_gc_perc, tuple2string
import chonk.Breakpoint as Breakpoint
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import pysam, csv

def bp_contigs(bp_file):
	# find contigs used in breakpoint file
	bp_contigs = ()
	bp_file_fh = open(bp_file, 'r')
	for row in csv.reader(bp_file_fh, dialect="excel-tab"):
		if row[0].startswith('#'): continue
		if row[0] not in bp_contigs:
			bp_contigs += (row[0],)
	return bp_contigs

def doc_coverage(AlnFile, Meta, contig, start, end):
	meta_doc = Meta.doc[contig]

	# find the doc fold change for the sv region
	read_length_sum = 0
	region_size = end - start
	for Aln in AlnFile.bam.fetch(reference = contig, start = start, end = end):
		# Some Aln.reference_start == None, which causes error
		try:
			midpoint = alignment_midpoint(Aln)
		except:
			continue
		if start <= midpoint <= end:
			read_length_sum += Aln.reference_length
	if region_size != 0:
		sv_doc = read_length_sum / region_size
	if meta_doc != 0:
		doc_fc = sv_doc / meta_doc

	return doc_fc

def gcbin_coverage(AlnFile, Meta, contig, start, end):
	fasta = Meta.fasta

	# choose window size based on the size of the region
	region_size = end-start
	window_size = Meta.window_lengths[0]
	if region_size >= (2 * Meta.window_lengths[1]):
		window_size = Meta.window_lengths[1]

	# find the mean and STD binned coverage in the SV region in relation to GC content
	Region = Welford()

	read_count = 0
	region_windows = create_windows(window_size, start, end)
	for window in region_windows:
		win_start, win_end = window[0], window[1]

		# find the gc bin for the window
		gc_bin = get_gc_bin(get_gc_perc(fasta, window_size, contig, win_start, win_end))

		# count up the number of reads in the window
		for Aln in AlnFile.bam.fetch(reference = contig, start = win_start, end = win_end):
			# Some Aln.reference_start == None, which causes error
			try:
				midpoint = alignment_midpoint(Aln)
			except:
				continue
			if start <= midpoint <= end:
				read_count += 1

		# calculate the fold change of read count in window to the count of the null with similar gc content
		key = tuple2string( (contig, window_size, gc_bin) )
		gc_fc = read_count / Meta.gc_rc[key]

		# add the gc fold change value to list of all windows in the sv region
		Region.update(gc_fc)

	# return mean and STD for all the windows in the region
	return (Region.mean, Region.std)