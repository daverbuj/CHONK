import pysam, csv, sys, os, numpy as np, subprocess as sp

bam = '/projects/ps-gleesonlab5/reference_panels/1kgp/platinumgenomes/HG02568.wgs.ILLUMINA.bwa.GWD.high_cov_pcr_free.20140203.bam'
chrom = '21'
raw_mask = '/home/daverbuj/genomes/hg19.sv.mask.merged.parexcluded.bed'
#mask = '/home/daverbuj/scripts/HG02568.chr21.filtered.masked.bed'
fasta = '/home/dantakli/ref/human_g1k_v37.fasta'
genome = 'HG02568'
# USER WILL NEED TO ADD THEIR FASTA FILE AS AN ARGUMENT~~~
# AND GENOME NAME ~~~


AlnFile = pysam.AlignmentFile(bam,'r')

# Preparing mask file for extraction
# * Create chrom region bed file
chrom_length = AlnFile.lengths[int(chrom)-1]
cwd = os.getcwd()
chrom_bed = "{}/{}.{}.bed".format(cwd, genome, chrom)
with open(chrom_bed, 'wt') as bed_file:
	tsv_writer = csv.writer(bed_file, delimiter='\t')
	tsv_writer.writerow([chrom, 0, chrom_length])
	bed_file.close()
# * Inverse of mask to chrom of interest
masked_bed = chrom_bed[:-4] + '.masked.bed'
mask_abs_path = os.path.abspath(raw_mask)
mask_cmd = "awk subtract -a {} -b {} > {}".format(chrom_bed, mask_abs_path, masked_bed)
sp.check_output(mask_cmd, shell = True)
# * Filter out regions smaller than 10kb
filtered_bed = masked_bed[:-4] + '.filtered.bed'
filter_cmd = "awk 'BEGIN {OFS='\t'}; $3 - $2 > 10000 {print $1, $2, $3}' {} > {}".format(masked_bed, filtered_bed)
sp.check_output(filter_cmd, shell = True)
# * Make windows of set sizes (100, 1000)
fasta_abs_path = os.path.abspath(fasta)
for window_size in [100,1000]:
	windowed_bed = filtered_bed[:-4] + '.window' + str(window_size) + '.bed'
	win_cmd = "bedtools makewindows -w {} -b {} > {}".format(str(window_size), filtered_bed, windowed_bed)
	sp.check_output(win_cmd, shell = True)
# Get GC content for each window
	gc_bed = windowed_bed[:-4] + '.gc.bed'
	gc_cmd = "bedtools nuc -fi {} -bed {} | cut -f 1-3,5 > {}".format(fasta_abs_path, windowed_bed, gc_bed)
	sp.check_output(gc_cmd, shell = True)

#mask_gc_dict = {}
#mask_gc_dict[window_size] = gc_mask_window

AlnFile = pysam.AlignmentFile(bam,'r')

# Getting Chromosome DOC, average read length, median MPD, and MAD
#read_lengths_raw = []
#mpd = []
#read_count = 0
#bp_count = 0
#is_broken = False
#with open(mask) as tsv:
#	for mask_loc in csv.reader(tsv, dialect="excel-tab"):
#		if is_broken == True: break
#		mask_chrom = mask_loc[0]
#		mask_start = mask_loc[1]
#		mask_end = mask_loc[2]
#		if chrom == mask_chrom:
#			for Aln in AlnFile.fetch(reference = chrom, start = int(mask_start), end = int(mask_end)):
#				# if midpoint within the region(need to add later)
#				if read_count > 5e6:
#					bp_count += Aln.reference_end - int(mask_start)
#					is_broken = True
#					break
#				read_lengths_raw.append(Aln.reference_length)
#				read_count += 1
#				dist = abs(Aln.template_length)
#				mpd.append(dist)
#			if is_broken == False: bp_count += (int(mask_end) - int(mask_start))
#	read_lengths = [x for x in read_lengths_raw if type(x) == int]
#	avg_length = np.mean(read_lengths)
#	chr_doc = (avg_length * read_count) / bp_count
#	med_mpd = np.median(mpd)
#	abs_diff = []
#	for x in mpd:
#		abs_diff.append(abs(x - med_mpd))
#	mad = np.median(abs_diff)
#
#	metadata = [chrom, avg_length, chr_doc, med_mpd, mad]


# Getting coverage for the gc bins
read_lengths_raw = []
mpd = []
read_count = 0
bp_count = 0
is_broken = False
gc_bin_dict = {}
for bin_size in mask_gc_dict:
	gc_bin_dict[bin_size] = {}
	with open(mask_gc_dict[bin_size]) as tsv:
		for window in csv.reader(tsv, dialect="excel-tab"):
			if is_broken == True: break
			if chrom == window[0]:
				win_chrom = window[0]
				win_start = window[1]
				win_end = window[2]
				win_gc = window[3][:4]
				for Aln in AlnFile.fetch(reference = chrom, start = int(win_start), end = int(win_end)):
					if Aln.is_unmapped: continue
					midpoint = ((Aln.reference_end - Aln.reference_start) / 2) + Aln.reference_start
					if int(win_start) <= midpoint <= int(win_end):
						read_count += 1
						read_lengths_raw.append(Aln.reference_length)
						if read_count > 2e5:
							bp_count += Aln.reference_end - int(win_start)
							is_broken = True
							break
				if is_broken == False: bp_count += (int(win_end) - int(win_start))
				read_lengths = [x for x in read_lengths_raw if type(x) == int]
				avg_length = np.mean(read_lengths)
				win_doc = (avg_length * read_count) / bp_count
				if not win_gc in gc_bin_dict[bin_size]:
					gc_bin_dict[bin_size][win_gc] = [win_doc]
				else:
					gc_bin_dict[bin_size][win_gc].append(win_doc)
	for gc_perc in gc_bin_dict[bin_size]:
		gc_bin_dict[bin_size][gc_perc] = [np.mean(gc_bin_dict[bin_size][gc_perc]), np.median(gc_bin_dict[bin_size][gc_perc]), len(gc_bin_dict[bin_size][gc_perc])]
print(gc_bin_dict)


# Deleting BED files that were made
rm_cmd = "rm {} {} {}".format(chrom_bed, masked_bed, filtered_bed)
sp.check_output(rm_cmd, shell = True)
for window_size in [100, 1000]:
		windowed_bed = filtered_bed[:-4] + '.window' + str(window_size) + '.bed'
		gc_bed = windowed_bed[:-4] + '.gc.bed'
		rm_cmd = "rm {} {}".format(windowed_bed, gc_bed)
		sp.check_output(rm_cmd, shell = True)
