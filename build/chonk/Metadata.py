#!/usr/bin/env python3
from chonk.Bam import Bam
import chonk.Backend as Backend
import pysam, sys, csv
import numpy as np

# Extracting average read length and depth of coverage across chromosome
def doc():
	read_lengths_raw = [] # Compiles read_length of all reads in the chromosome
	read_count = 0 # Counts up each read in the chromosome
	mpd = []
	

	return avg_length, chrom_doc


#################################################

# Need bamfile, cn2 file, output directory, chrom
def metadataExtraction(Args):

	# check if alignment file exists
	AlnFile = Bam(Args.i)
	# ensure the chromosome prefix matches the one in the aln file
	chrom = Backend.check_chrom(Args.r,AlnFile.chrom_flag)
	# output file
	out = open(Args.o,'w')
	out.write('#chrom\tavgLength\tdoc\tmpd\tmad\n')

	read_lengths_raw = []
	read_count = 0
	bp_count = 0
	mpd = []

	with open(Args.c) as tsv:
		for cn2_loc in csv.reader(tsv, dialect="excel-tab"): # CN2 LOCATIONS
			cn2_chrom = cn2_loc[0]
			if 'chr' in chrom and 'chr' not in cn2_chrom:
				cn2_chrom = 'chr' + cn2_chrom
			elif 'chr' not in chrom and 'chr' in cn2_chrom:
				cn2_chrom.replace('chr','')
			cn2_start = cn2_loc[1]
			cn2_end = cn2_loc[2]
			if chrom == cn2_chrom:
				bp_count += (int(cn2_end) - int(cn2_start))
				for Aln in AlnFile.bam.fetch(reference = chrom, start = int(cn2_start), end = int(cn2_end)): # READS IN CN2 LOCATIONS
					read_lengths_raw.append(Aln.reference_length)
					read_count += 1
					if (Aln.is_paired # Read is Paired
						and Aln.next_reference_name == Aln.reference_name # Same chr
						and  Aln.is_reverse != Aln.mate_is_reverse # Opposite strands
						and not Aln.has_tag('SA') # Is not a split read         DO WE NEED THIS FOR META EXTRACTION? ~~~~~~~~~~~~~~
						and not Aln.is_reverse # Each mate/pair once; from positive strand
						):
						dist = abs(Aln.template_length)
						mpd.append(dist)
		read_lengths = [x for x in read_lengths_raw if type(x) == int]
		avg_length = np.mean(read_lengths)       # BETTER TO DYNAMICALLY CALCULATE MEAN AS TO NOT TAKE UP too much MEMORY? ~~~~~~~~~~~~~~~~
		chr_doc = (avg_length * read_count) / bp_count
		med_mpd = np.median(mpd)
		abs_diff = []
		for x in mpd:
			abs_diff.append(abs(x - med_mpd))
		mad = np.median(abs_diff)
		metadata = [chrom, avg_length, chr_doc, med_mpd, mad]
		out.write('\t'.join(map(str,metadata))+'\n')
	out.close()
