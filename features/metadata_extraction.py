import pysam
import sys
import numpy as np
import subprocess as sp
import os
import csv

bam_directory = sys.argv[1] # Here we use the directory with the .bam files
if bam_directory[-1] != '/':
	bam_directory = bam_directory + '/' # Check that directory argument given ends with a forwardslash
cn2_bed = sys.argv[2] # cn2.bed file
output_directory = sys.argv[3] # Directory we want to write meta files to
if output_directory[-1] != '/':
	output_directory = output_directory + '/'
chrom = sys.argv[4]

	
# Finding the DOC for each chromosome
	# L (Average read length) == avg_length
	# N (Number of reads) == read_count
	# G (Number of basepairs) == bp_count
	# Chromosome doc = L X N / G
for filename_r in os.listdir(bam_directory): # GENOME
	if filename_r.endswith('.bam'): 
		bamfile = filename_r
		genome = bamfile[:7]
		bamfile_abspath = bam_directory + bamfile
		filename_abspath_w = output_directory + genome + '.meta.tsv'
		bam = pysam.AlignmentFile(bamfile_abspath,'r')
		with open(filename_abspath_w, 'w') as out_file: 
			tsv_writer = csv.writer(out_file, delimiter='\t')
			read_lengths_raw = [] # Compiles read_length of all reads in the chromosome
			read_count = 0 # Counts up each read in the chromosome
			mpd = []
			with open(cn2_bed) as tsv:
				for cn2_loc in csv.reader(tsv, dialect="excel-tab"): # CN2 LOCATIONS
					cn2_chrom = cn2_loc[0]
					cn2_start = cn2_loc[1]
					cn2_end = cn2_loc[2]
					if chrom == cn2_chrom:
						for Aln in bam.fetch(reference = chrom, start = int(cn2_start), end = int(cn2_end)): # READS IN CN2 LOCATIONS
							read_lengths_raw.append(Aln.reference_length)
							read_count += 1
							if (Aln.is_paired # Read is Paired
								and Aln.next_reference_name == Aln.reference_name # Same chr
								and  Aln.is_reverse != Aln.mate_is_reverse # Opposite strands
								and not Aln.has_tag('SA') # Is not a split read
								and Aln.is_read1 # Each mate/pair once
								):
								dist = Aln.template_length
								mpd.append(abs(dist))
				read_lengths = [x for x in read_lengths_raw if x is not None]
				avg_length = np.mean(read_lengths)
				# Finding number of bps in each chromosome using the cn2 bed file
				command = "awk '{ print $1,$3 - $2 }' "+cn2_bed+" | awk '$1 == "+chrom+"' | awk '{s+=$2} END {print s}'"
				bp_count_raw = sp.check_output(command,shell=True)
				bp_count = int(str(bp_count_raw)[2:-3])
				chr_doc = (avg_length * read_count) / bp_count
				med_mpd = np.median(mpd)
				abs_diff = []
				for x in mpd:
					abs_diff.append(abs(x - med_mpd))
				mad = np.median(abs_diff)
				tsv_writer.writerow([chrom, avg_length, chr_doc, med_mpd, mad])
