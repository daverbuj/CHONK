import pysam
import sys
import numpy as np
import subprocess as sp
import os
import csv

directory = sys.argv[1] # Here we use the directory with the .bam files
if directory[-1] != '/':
	directory = directory + '/' # Check that directory argument given ends with a forwardslash
cn2_bed = sys.argv[2]
sv_bed = sys.argv[3]

gen_doc = {} # Dictionary for the coverage of the .bam files in the working directory. 
			 # Key = genome name (HG00174), Value = list of tuples that include coverage and avg_length per chromosome by index-value+1 ([4.2,4.1,3.9,4.2,...]). (HG00174 chr 3 doc == gen_doc[HG00174][2][0]]

# Finding the DOC for each chromosome
	# L (Average read length) == avg_length
	# N (Number of reads) == read_count
	# G (Number of basepairs) == bp_count
	# Chromosome doc = L X N / G
for filename in os.listdir(directory):
	if filename.endswith(".bam"): 
		bamfile = filename
		genome = bamfile[:7]
		bamfile_abs_path = directory + bamfile
		bam = pysam.AlignmentFile(bamfile_abs_path,'r')
		for chrom in range(1,23): # loops through each chromosome
			read_lengths_raw = [] # Compiles read_length of all reads in the chromosome
			read_count = 0 # Counts up each read in the chromosome
			with open(cn2_bed) as tsv:
    			for cn2_loc in csv.reader(tsv, dialect="excel-tab"): # loops through each cn2 location 
    				cn2_chrom = cn2_loc[0]
    				cn2_start = cn2_loc[1]
    				cn2_end = cn2_loc[2]
    				if str(chrom) == cn2_chrom:
						for Aln in bam.fetch(reference = str(chrom), start = cn2_start, end = cn2_end): # loops through the alignments in the cn2 locations
							read_lengths_raw.append(Aln.reference_length)
							read_count += 1 # is this sufficient for finding number of reads? confused by read + mate_1 comment
			read_lengths = [x for x in read_lengths_raw if x is not None]
			avg_length = np.mean(read_lengths)
			# Finding number of bps in each chromosome using the cn2 bed file
			command = "awk '{ print $1,$3 - $2 }' "+cn2_bed+" | awk '$1 == "+str(chrom)+"' | awk '{s+=$2} END {print s}'"
			bp_count_raw = sp.check_output(command,shell=True)
			bp_count = int(str(bp_count_raw)[2:-3])
			chr_doc = (avg_length * read_count) / bp_count
			if genome not in gen_doc:
				gen_doc[genome] = [(chr_doc, avg_length)]
			else:
				gen_doc[genome].append((chr_doc, avg_length))


# Finding the DOC for each SV
	# L (Average read length) == sv_avg_read_length
	# N (Number of reads) == sv_read_count
	# G (Number of basepairs) == sv_bp_count
	# Chromosome doc = L X N / G
with open(sv_bed) as tsv:
    for sv in csv.reader(tsv, dialect="excel-tab"):
        sv_chrom = sv[0]
        sv_start =sv[1]
        sv_end = sv[2]
        sv_type = sv[3]
        sv_genome = sv[4]
        if sv[5] == '0|0':
        	sv_genotype = '0'
        elif sv[5] == '1|0' or sv[5] == '0|1':
        	sv_genotype = '1'
        elif sv[5] == '1|1':
        	sv_genotype = '2'
        sv_genotype = sv[5]
        sv_read_count = 0
        sv_bp_count = sv_end - sv_start
		sv_avg_read_length = gen_doc[sv_genome][int(sv_chrom)-1][1]
		for Aln in bam.fetch(reference = sv_chrom, start = sv_start, end = sv_end):
			sv_read_count += 1
		sv_doc = (sv_avg_read_length * sv_read_count) / sv_bp_count
		doc_fc = sv_doc / gen_doc[sv_genome][int(sv_chrom)-1][0]
		print(sv_type, sv_genotype, sv_doc, doc_fc)
