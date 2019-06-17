#!/usr/bin/env python3
import sys, pysam, os, json

def discordantpe(Aln, meta, chrom):

	# convert the metadata tsv file to a dict
	with open(meta) as json_file:
		data = json.load(json_file)
		chr_metadata = data[1]
		med_mpd = chr_metadata[chrom]['meta'][2]
		mad = chr_metadata[chrom]['meta'][3]

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
#		if dist < (float(med_mpd) - (3.5 * float(mad))):
#				sv = 'INS'
#				breakpoint_start = 
#				breakpoint_end = 
#				breaks.append([chrom, breakpoint_start, breakpoint_end, svtype])
				return (chrom, breakpoint_start, breakpoint_end, svtype, 'DPE')
			return (chrom, breakpoint_start, breakpoint_end, svtype, 'DPE')
	elif (Aln.is_paired # Read is Paired
		and Aln.next_reference_name == Aln.reference_name # Same chr
		and Aln.is_reverse == Aln.mate_is_reverse # Same strand
		and Aln.is_read1 # Each mate/pair once; from positive strand
		):
		dist = abs(Aln.template_length)
		if dist > (float(med_mpd) + (3.5 * float(mad))):
			svtype = 'INV'
			if not Aln.is_reverse and not Aln.mate_is_reverse: # if both reads are on the positive strand
				if Aln.reference_end < (Aln.next_reference_start + Aln.query_length):
					breakpoint_start = Aln.reference_end
					breakpoint_end = (Aln.next_reference_start + Aln.query_length)
				else:
					breakpoint_start = (Aln.next_reference_start + Aln.query_length)
					breakpoint_end = Aln.reference_end
				return (chrom, breakpoint_start, breakpoint_end, svtype, 'DPE')
			elif Aln.is_reverse and not Aln.mate_is_reverse: # if both reads are on the negative strand
				if Aln.reference_start < Aln.next_reference_start:
					breakpoint_start = Aln.reference_start
					breakpoint_end = Aln.next_reference_start
				else:
					breakpoint_start = Aln.next_reference_start
					breakpoint_end = Aln.reference_start
				return (chrom, breakpoint_start, breakpoint_end, svtype, 'DPE')

