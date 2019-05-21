#!/usr/bin/env python3
import sys, pysam, os, csv

def discordantpe(Aln, meta, chrom):

	# convert the metadata tsv file to local variables
	with open(meta) as tsv:
		for metadata in csv.reader(tsv, dialect="excel-tab"):
			med_mpd = metadata[3]
			mad = metadata[4]

	# looking through bam 
	if (Aln.is_paired # Read is Paired
		and Aln.next_reference_name == Aln.reference_name # Same chr
		and  Aln.is_reverse != Aln.mate_is_reverse # Opposite strands
		and not Aln.is_reverse # Each mate/pair once; from positive strand
		):
		dist = abs(Aln.template_length)
		if dist > (float(med_mpd) + (3.5 * float(mad))):
			svtype = 'DEL'
			breakpoint_start = Aln.reference_end
			breakpoint_end = Aln.next_reference_start
			if Aln.reference_start > Aln.next_reference_start:
				svtype = 'DUP'
				breakpoint_start = Aln.next_reference_start
				breakpoint_end = Aln.reference_end
			return (chrom, breakpoint_start, breakpoint_end, svtype)
#		if dist < (float(med_mpd) - (3.5 * float(mad))):
#				sv = 'INV'
#				breakpoint_start = 
#				breakpoint_end = 
#				breaks.append([chrom, breakpoint_start, breakpoint_end, svtype])
