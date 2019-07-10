#!/usr/bin/env python3
from chonk.Bam import Bam
from chonk.Metadata import Metadata
import chonk.Coverage_Features as Coverage
import chonk.Breakpoint as Breakpoint
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import chonk.SupportingFragments_Features as SupportingFragments
import pysam, json, os, sys, csv


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
	contigs = Breakpoint.selected_contigs(Args.r, Coverage.bp_contigs(Args.b), AlnFile)

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
		# skip the first line (header)
		if sv[0].startswith('#'): continue

		# this is for validated SVs / training set
		try:
			contig, start, end, svtype, cistart, ciend, genome, genotype = sv[0], int(sv[1]), int(sv[2]), sv[3], sv[4], sv[5], sv[6], sv[7]
		# this is for unvalidated SVs
		except:
			contig, start, end, svtype, cistart, ciend = sv[0], int(sv[1]), int(sv[2]), sv[3], '.', '.'

		## extract coverage features
		# sv region
		sv_doc_fc = Coverage.doc_coverage(AlnFile, Meta, contig, start, end)
		sv_gc_mean, sv_gc_std = Coverage.gcbin_coverage(AlnFile, Meta, contig, start, end)
		# left flank
		lf_doc_fc = Coverage.doc_coverage(AlnFile, Meta, contig, start-1000, start)
		lf_gc_mean, lf_gc_std = Coverage.gcbin_coverage(AlnFile, Meta, contig, start-1000, end)
		# right flank
		rf_doc_fc = Coverage.doc_coverage(AlnFile, Meta, contig, end, end+1000)
		rf_gc_mean, rf_gc_std = Coverage.gcbin_coverage(AlnFile, Meta, contig, end, end+1000)

		## extract supporting fragment features
		supporting_fragment_features = SupportingFragments.supporting_frags(AlnFile, Meta, contig, start, end, cistart, ciend, svtype)
		sf_ratio, split_ratio, disc_ratio, clip_ratio, sf_mapq_mean, sf_mapq_median, nonsf_baseq_mean, nonsf_baseq_median = supporting_fragment_features

		'''
		print(sv_doc_fc, sv_gc_mean, sv_gc_std, lf_doc_fc, lf_gc_mean, lf_gc_std, rf_doc_fc, rf_gc_mean, rf_gc_std)
		print('~~~~~')
		print(supporting_fragment_features)
		print('++++++++++++++++++++++++++++')'''




