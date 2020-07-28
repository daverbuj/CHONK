from chonk.Bam import Bam
from chonk.Metadata import Metadata
from chonk.Backend import alignment_midpoint, Welford, get_gc_bin, tuple2string
import pysam, json, os, sys, csv

def create_windows(window_size, seq):
	windows = []
	tmp_start = 0
	tmp_end = tmp_start + window_size
	for x in range( int( len(seq)/window_size ) ):
		windows.append( (tmp_start, tmp_end) )
		tmp_start = tmp_end
		tmp_end = tmp_start + (window_size)
	return windows

def seq_gc_perc(seq):
	length = len(seq)
	gc = 0
	for base in seq:
		if base.upper() == 'G' or base.upper() == 'C':
			gc += 1
	return gc/length

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
		# skip duplicates
		if Aln.is_duplicate and Aln.is_unmapped: continue
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

def gcbin_coverage(AlnFile, Meta, contig, start, end, sequence):
	# choose window size based on the size of the region
	window_size = region_window_size(Meta.window_lengths, start, end)

	# find the mean and STD binned coverage in the SV region in relation to GC content
	Region = Welford()

	region_windows = create_windows(window_size, sequence)
	for window in region_windows:
		local_start, local_end = window[0], window[1]
		window_seq = sequence[local_start:local_end]

		# find the gc bin for the window
		gc_bin = get_gc_bin(seq_gc_perc(window_seq))

		# count up the number of reads in the window
		win_start, win_end = start+local_start, start+local_end    # 1 bp overlap in windows ((0,4),(4,8),(8,12)....)

		read_count = 0
		for Aln in AlnFile.bam.fetch(reference = contig, start = win_start, end = win_end):
			# skip duplicates
			if Aln.is_duplicate and Aln.is_unmapped: continue
			# Some Aln.reference_start == None, which causes error
			try:
				midpoint = alignment_midpoint(Aln)
			except:
				continue
			if start <= midpoint <= end:
				read_count += 1

		# calculate the fold change of read count in window to the count of the null with similar gc content
		key = tuple2string( (contig, window_size, gc_bin) )

		try:
			if int(Meta.gc_rc[key]) != 0:
				gc_fc = read_count / Meta.gc_rc[key]
			else:
				gc_fc = 0
		except KeyError:
			gc_fc = 0

		# add the gc fold change value to list of all windows in the sv region
		Region.update(gc_fc)

	# return mean and STD for all the windows in the region
	return (Region.mean, Region.std)

def region_window_size(window_lengths, start, end):
	region_size = end-start
	window_size = window_lengths[0]
	if region_size >= (2 * window_lengths[1]):
		window_size = window_lengths[1]
	return window_size

################################################################################################
#MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN  
################################################################################################

def cov_features(AlnFile, Meta, chrom, start, end, sv_seq, lf_seq, rf_seq):

	## extract coverage-based features
	# sv region
	sv_doc_fc = doc_coverage(AlnFile, Meta, chrom, start, end)
	sv_gc_mean, sv_gc_std = gcbin_coverage(AlnFile, Meta, chrom, start, end, sv_seq)
	# left flank
	lf_doc_fc = doc_coverage(AlnFile, Meta, chrom, start-1000, start)
	lf_gc_mean, lf_gc_std = gcbin_coverage(AlnFile, Meta, chrom, start-1000, start, lf_seq)
	# right flank
	rf_doc_fc = doc_coverage(AlnFile, Meta, chrom, end, end+1000)
	rf_gc_mean, rf_gc_std = gcbin_coverage(AlnFile, Meta, chrom, end, end+1000, rf_seq)

	return (sv_doc_fc, sv_gc_mean, sv_gc_std, lf_doc_fc, lf_gc_mean, lf_gc_std, rf_doc_fc, rf_gc_mean, rf_gc_std)