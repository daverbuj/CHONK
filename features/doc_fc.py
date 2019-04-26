import pysam
import sys
import numpy as np
import subprocess as sp
import os
import csv

directory = sys.argv[1] # Here we use the directory with the .bam files that have been intersected with the cn2 .bed file
if directory[-1] != '/':
	directory = directory + '/' # Check that directory argument given ends with a forwardslash
ctrl_bed = sys.argv[2]
sv_bed = sys.argv[3]

gen_doc = {} # Dictionary for the coverage of the .bam files in the working directory. 
			 # Key = genome name (HG00174), Value = list of tuples that include coverage and avg_length per chromosome by index-value+1 ([4.2,4.1,3.9,4.2,...]). (HG00174 chr 3 doc == gen_doc[HG00174][2])

# Finding the DOC for each chromosome
	# L (Average read length) == avg_length
	# N (Number of reads) == read_count
	# G (Number of basepairs) == bp_count
	# Chromosome doc = L X N / G
for filename in os.listdir(directory):
	if filename.endswith("cn2.bam"): 
		bamfile = filename
		genome = bamfile[:7]
		bamfile_abs_path = directory + bamfile
		bam = pysam.AlignmentFile(bamfile_abs_path,'r')
		for chrom in range(1,23):
			read_lengths_raw = [] # Compiles read_length of all reads in the chromosome
			read_count = 0 # Counts up each read in the chromosome
			for Aln in bam.fetch(reference = str(chrom)):
				read_lengths_raw.append(Aln.reference_length)
				read_count += 1 # is this sufficient for finding number of reads? confused by read + mate_1 comment
			read_lengths = [x for x in read_lengths_raw if x is not None]
			avg_length = np.mean(read_lengths)
			# Finding number of bps in each chromosome using the cn2 bd file
			command = "awk '{ print $1,$3 - $2 }' "+ctrl_bed+" | awk '$1 == "+str(chrom)+"' | awk '{s+=$2} END {print s}'"
			bp_count_raw = sp.check_output(command,shell=True)
			bp_count = int(str(bp_count_raw)[2:-3])
			if genome not in gen_doc:
				gen_doc[genome] = [((avg_length * read_count) / bp_count, avg_length)]
			else:
				gen_doc[genome].append(((avg_length * read_count) / bp_count, avg_length))
print(gen_doc)


# Finding the DOC for each SV
	# L (Average read length) == 
	# N (Number of reads) == 
	# G (Number of basepairs) == 
	# Chromosome doc = L X N / G
with open(sv_bed) as tsv:
    for sv in csv.reader(tsv, dialect="excel-tab"):
        sv_chrom = sv[0]
        sv_start =sv[1]
        sv_end = sv[2]
        sv_type = sv[3]
        sv_geno = sv[4]
		for Aln in bam.fetch(reference=sv_chrom, start=sv_start, end=sv_end):
			chrom = Aln.reference_name
			lpos = Aln.reference_start+1 #Aln.reference_start is in 0-base
			rpos = Aln.reference_end
			mapq = Aln.mapping_quality
			qpos = Aln.query_alignment_start+1 # 1-base position


			#print(chrom, qpos, lpos, rpos)
			#print(sv_start-lpos,sv_end-rpos)
