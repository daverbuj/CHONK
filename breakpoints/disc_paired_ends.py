import sys
import pysam
import os
import csv

bam_directory = sys.argv[1]
if bam_directory[-1] != '/':
	bam_directory = bam_directory + '/'
meta_directory = sys.argv[2]
if meta_directory[-1] != '/':
	meta_directory = meta_directory + '/'
chrom = sys.argv[3]

# Metafile to dict
meta_dict = {}
for filename_rm in os.listdir(meta_directory):
	if filename_rm.endswith('.meta.tsv'): 
		metafile = filename_rm
		meta_genome = metafile[:7]
		metafile_abspath = meta_directory + metafile
		with open(metafile_abspath) as tsv:
			for meta in csv.reader(tsv, dialect="excel-tab"):
				meta_chrom = meta[0]
				meta_med_mpd = meta[3]
				meta_mad = meta[4]
				if chrom == meta_chrom:
					meta_dict[meta_genome] = [meta_med_mpd, meta_mad]


for filename_rb in os.listdir(bam_directory):
	if filename_rb.endswith('.bam'):
		bamfile = filename_rb
		genome = bamfile[:7]
		bamfile_abspath = bam_directory + bamfile
		bam = pysam.AlignmentFile(bamfile_abspath,'r')
		for Aln in bam.fetch(reference = chrom):
			if (Aln.is_paired # Read is Paired
				and Aln.next_reference_name == Aln.reference_name # Same chr
				and  Aln.is_reverse != Aln.mate_is_reverse # Opposite strands
				and not Aln.has_tag('SA') # Is not a split read
				and Aln.is_read1 # Each mate/pair once
				):
				dist = Aln.template_length
				med_mpd = meta_dict[genome][0]
				mad = meta_dict[genome][1]
				if ((not Aln.is_reverse and dist > 0) or (Aln.is_reverse and dist < 0)) and dist > (float(med_mpd) + (999 * float(mad))):
					sv = 'Del'
					breakpoint_start = Aln.reference_end
					breakpoint_end = Aln.next_reference_start
					print(chrom,'\t',breakpoint_start,'\t',breakpoint_end,'\t',sv,'\t',genome,'\t',Aln.query_name, sep='')
				if ((not Aln.is_reverse and dist < 0) or (Aln.is_reverse and dist > 0)) and dist < (float(med_mpd) - (999 * float(mad))):
					sv = 'Dup'
					breakpoint_start = Aln.next_reference_start
					breakpoint_end = Aln.reference_end
					print(chrom,'\t',breakpoint_start,'\t',breakpoint_end,'\t',sv,'\t',genome,'\t',Aln.query_name, sep='')
