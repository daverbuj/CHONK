import pysam
import sys
import numpy as np
import subprocess as sp

bam_fh = sys.argv[1]
bed = sys.argv[2]

chrom_doc = {}
#file = open(file_out,'w')

bam = pysam.AlignmentFile(bam_fh,'r')

# Finding average read length and read count of chromosomes in .bam file
for chrom in range(1,23):
	read_lengths_raw = []
	read_count = 0
	for Aln in bam.fetch(reference = str(chrom)):
		read_lengths_raw.append(Aln.reference_length)
		read_count += 1 # is this sufficient for finding number of reads? confused by read + mate_1 comment
	read_lengths = [x for x in read_lengths_raw if x is not None]
	avg_length = np.mean(read_lengths)
	#print(chrom,avg_length,counter)
	# Finding basepairs in each chromosome
	command = "awk '{ print $1,$3 - $2 }' "+bed+" | awk '$1 == "+str(chrom)+"' | awk '{s+=$2} END {print s}'"
	bp_count_raw = subprocess.check_output(command,shell=True)
	bp_count = int(str(bp_count_raw)[2:-3])
	print(chrom, avg_length, read_count, bp_count)
	#chrom_doc[chrom] = [avg_length,counter,bp_count]





	#file.write(str(chrom)+'\t'+str(avg_length)+'\n') # I'm writing out to a file, but would it be better to store as a dict?
#file.close()
