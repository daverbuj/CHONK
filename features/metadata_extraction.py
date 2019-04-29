import pysam
import sys
import numpy as np
import subprocess as sp
import os
import csv

directory = sys.argv[1] # Here we use the directory with the .bam files
if directory[-1] != '/':
	directory = directory + '/' # Check that directory argument given ends with a forwardslash
cn2_bed = sys.argv[2] # cn2.bed file
output_directory = sys.argv[3] # Directory we want to write meta files to

# Finding the DOC for each chromosome
	# L (Average read length) == avg_length
	# N (Number of reads) == read_count
	# G (Number of basepairs) == bp_count
	# Chromosome doc = L X N / G
for filename_r in os.listdir(directory): # GENOME
	if filename_r.endswith('.bam'): 
		bamfile = filename_r
		genome = bamfile[:7]
		bamfile_abspath = directory + bamfile
		filename_abspath_w = output_directory + genome + '.meta.tsv'
		bam = pysam.AlignmentFile(bamfile_abspath,'r')
		with open(filename_abspath_w, 'w') as out_file: 
			tsv_writer = csv.writer(out_file, delimiter='\t')
			for chrom in range(1,23): # CHROMSOME
				read_lengths_raw = [] # Compiles read_length of all reads in the chromosome
				read_count = 0 # Counts up each read in the chromosome
				with open(cn2_bed) as tsv:
					for cn2_loc in csv.reader(tsv, dialect="excel-tab"): # CN2 LOCATIONS
						cn2_chrom = cn2_loc[0]
						cn2_start = cn2_loc[1]
						cn2_end = cn2_loc[2]
						if str(chrom) == cn2_chrom:
							for Aln in bam.fetch(reference = str(chrom), start = int(cn2_start), end = int(cn2_end)): # loops through the alignments in the cn2 locations
								read_lengths_raw.append(Aln.reference_length)
								read_count += 1 # is this sufficient for finding number of reads? confused by read + mate_1 comment
				read_lengths = [x for x in read_lengths_raw if x is not None]
				avg_length = np.mean(read_lengths)
				# Finding number of bps in each chromosome using the cn2 bed file
				command = "awk '{ print $1,$3 - $2 }' "+cn2_bed+" | awk '$1 == "+str(chrom)+"' | awk '{s+=$2} END {print s}'"
				bp_count_raw = sp.check_output(command,shell=True)
				bp_count = int(str(bp_count_raw)[2:-3])
				chr_doc = (avg_length * read_count) / bp_count
				tsv_writer.writerow([chrom, avg_length, chr_doc])

