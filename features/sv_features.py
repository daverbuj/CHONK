import pysam
import sys
import os
import csv

meta_directory = sys.argv[1] # Here we use the directory with the genome metadata tsv files
if meta_directory[-1] != '/':
	meta_directory = meta_directory + '/' # Check that directory argument given ends with a forwardslash
bam_directory = sys.argv[2]
if bam_directory[-1] != '/':
	bam_directory = bam_directory + '/'
sv_bed = sys.argv[2]


# Finding the DOC for each SV
	# L (Average read length) == sv_avg_read_length
	# N (Number of reads) == sv_read_count
	# G (Number of basepairs) == sv_bp_count
	# SV doc = L X N / G
# SV BED (sv_bp_count)
with open(sv_bed) as tsv:
	for sv in csv.reader(tsv, dialect='excel-tab'):
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
		else:
			sv_genotype = '?'
		sv_read_count = 0
		sv_bp_count = int(sv_end) - int(sv_start)
		# BAM (sv_read_count)
		for filename in os.listdir(bam_directory):
			if filename.startswith(sv_genome) and filename.endswith('.bam'):
				bamfile = filename
				bamfile_abs_path = bam_directory + bamfile
				bam = pysam.AlignmentFile(bamfile_abs_path,'r')
				for Aln in bam.fetch(reference = sv_chrom, start = int(sv_start), end = int(sv_end)):
					sv_read_count += 1
		# META (sv_avg_read_length, gen_chrom_doc)
		for meta_file in os.listdir(meta_directory):
			if meta_file.startswith(sv_genome) and meta_file.endswith('.meta.tsv'):
				meta_file_abspath = meta_directory + meta_file
				with open(meta_file_abspath) as tsv:
					for gen_chrom_meta in csv.reader(tsv, dialect='excel-tab'):
						if sv_chrom == gen_chrom_meta[0]:
							sv_avg_read_length = gen_chrom_meta[1]
							gen_chrom_doc = gen_chrom_meta[2]
		sv_doc = (sv_avg_read_length * sv_read_count) / sv_bp_count
		doc_fc = sv_doc / gen_chrom_doc
		print(sv_type, sv_genotype, doc_fc, sv_doc, sv_start, sv_end)
