#!/usr/bin/env python3
from chonk.Bam import Bam
import chonk.Backend as Backend
import pysam, sys, csv, os, random, json, numpy as np, subprocess as sp
random.seed(42)

# inputs are bam file, mask file, fasta file, genome name, and chromosome
def metadataExtraction(Args):

	# check if alignment file exists
	AlnFile = Bam(Args.i)
	# mask file
	raw_mask = Args.e
	# fasta file
	fasta = Args.f
	# genome name
	genome = Args.g
	# ensure the chromosome prefix matches the one in the aln file
	chrom = Backend.check_chrom(Args.r,AlnFile.chrom_flag)

	# Create temporary files used to extract metadata
	# * Create genome chromosomal region bed file: chrom_bed = 'cwd/genome.chr.CHONKdel.bed'
	chrom_length = AlnFile.bam.lengths[int(chrom)-1]
	cwd = os.path.dirname(Args.o)
	chrom_bed = "{}/{}.{}.CHONKdel.bed".format(cwd, genome, chrom)
	with open(chrom_bed, 'wt') as cbed_file:
		cbed_writer = csv.writer(cbed_file, delimiter='\t')
		cbed_writer.writerow([chrom, 0, chrom_length])
		cbed_file.close()
	# * Inverse of mask to chrom of interest: masked_bed = 'cwd/genome.chr.CHONKdel.masked.bed'
	masked_bed = chrom_bed[:-4] + '.masked.bed'
	mask_abs_path = os.path.abspath(raw_mask)
	mask_cmd = "bedtools subtract -a {} -b {} > {}".format(chrom_bed, mask_abs_path, masked_bed)
	sp.check_output(mask_cmd, shell = True)
	# * Filter out regions smaller than 10kb: filtered_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.bed'
	filtered_bed = masked_bed[:-4] + '.filtered.bed'
	filter_cmd = """awk 'BEGIN {OFS="\t"}; $3 - $2 > 10000 {print $1, $2, $3}' """ + masked_bed + ' > ' + filtered_bed
	sp.check_output(filter_cmd, shell = True)
	# * Make windows of set sizes (25, 1000): windowed_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.window().bed'
	fasta_abs_path = os.path.abspath(fasta)
	for window_size in [25,1000]: #files made in this loop are used only for gc bins
		windowed_bed = filtered_bed[:-4] + '.window' + str(window_size) + '.bed'
		win_cmd = "bedtools makewindows -w {} -b {} > {}".format(str(window_size), filtered_bed, windowed_bed)
		sp.check_output(win_cmd, shell = True)
	# * Get GC content for each window: raw_gc_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.window().gc.bed'
		raw_gc_bed = windowed_bed[:-4] + '.gc.bed'
		gc_cmd = "bedtools nuc -fi {} -bed {} | cut -f 1-3,5 > {}".format(fasta_abs_path, windowed_bed, raw_gc_bed)
		sp.check_output(gc_cmd, shell = True)
	# * Create a bed file with grouped gc bins: grouped_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.window().gc.grouped.bed'
		grouped_bed = raw_gc_bed[:-4] + '.grouped.bed'
		with open(grouped_bed, 'wt') as gbed_file:
			gbed_writer = csv.writer(gbed_file, delimiter='\t')
			with open(raw_gc_bed) as tsv:
				for window in csv.reader(tsv, dialect="excel-tab"):
					win_chrom = window[0]
					win_start = window[1]
					win_end = window[2]
					try: win_gc = round(float(window[3]) * 100)
					except: continue
					gc_bin = 0
					if win_gc % 4 == 0 and win_gc != 0:
						gc_bin = int(win_gc / 4)
					elif win_gc % 4 == 1:
						gc_bin = int((win_gc + 3) / 4)
					elif win_gc % 4 == 2:
						gc_bin = int((win_gc + 2) / 4)
					elif win_gc % 4 == 3:
						gc_bin = int((win_gc + 1) / 4) #gc bin encompasses the percentages (gc_bin*4-3)--(gc_bin*4). example: gc_bin 4 == 13%-16%
					gbed_writer.writerow([win_chrom, win_start, win_end, gc_bin])
				gbed_file.close()


	# CHROMOSOME METADATA:
	# Getting Chromosome DOC, average read length, median MPD, and MAD
	mpd = []
	read_length_sum = 0
	read_count = 0
	bp_count = 0
	is_broken = False
	with open(filtered_bed) as tsv:
		for mask_loc in csv.reader(tsv, dialect="excel-tab"):
			if is_broken == True: break
			mask_chrom = mask_loc[0]
			mask_start = mask_loc[1]
			mask_end = mask_loc[2]
			if chrom == mask_chrom:
				for Aln in AlnFile.bam.fetch(reference = chrom, start = int(mask_start), end = int(mask_end)):
					if Aln.is_unmapped or type(Aln.reference_length) != int: continue
					midpoint = ((Aln.reference_end - Aln.reference_start) / 2) + Aln.reference_start
					if int(win_start) <= midpoint <= int(win_end):
						if read_count > 5e6:
							bp_count += Aln.reference_end - int(mask_start)
							is_broken = True
							break
						read_length_sum += Aln.reference_length
						read_count += 1
						dist = abs(Aln.template_length)
						mpd.append(dist)
				if is_broken == False: bp_count += (int(mask_end) - int(mask_start))
		avg_length = read_length_sum / read_count
		chr_doc = (avg_length * read_count) / bp_count
		med_mpd = np.median(mpd)
		abs_diff = []
		for x in mpd:
			abs_diff.append(abs(x - med_mpd))
		mad = np.median(abs_diff)
		chr_metadata = (genome, chrom, avg_length, chr_doc, med_mpd, mad)


	# GC METADATA
	# Shuffle the grouped_bedfile and write file that selects only the first 4,000 windows at each gcbin for 25bp-sized windows and 1000bp-sized windows
	gc_metadata = {}
	for win_size in [25, 1000]:
		num_wins = 10000
		grouped_bedfile = "{}/{}.{}.CHONKdel.masked.filtered.window{}.gc.grouped.bed".format(cwd, genome, chrom, win_size) 
		# * Shuffle the grouped_bedfile: shuffled_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.window().gc.grouped.shuffled.bed'
		line_count = len(open(grouped_bedfile).readlines())
		shuffled_bed = grouped_bedfile[:-4] + '.shuffled.bed'
		with open(grouped_bedfile) as infile, open(shuffled_bed, 'w') as outfile:
		    outfile.writelines(random.sample(infile.readlines(), line_count))
		infile.close()
		outfile.close()
		# * Write a bed file that only contains the first 10,000 instances of each gcbin: selected_bed = 'cwd/genome.chr.CHONKdel.masked.filtered.window().gc.grouped.shuffled.selected().bed'
		selected_bed = shuffled_bed[:-4] + '.selected{}.bed'.format(num_wins)
		wincount_dict = {} #used to count down remaining number of gc windows needed at each bin
		aln_count = {} #used to accumulate number of reads in each window for each gcbin
		for binid in range(26):
			wincount_dict[str(binid)] = num_wins
			aln_count[str(binid)] = []
		with open(selected_bed, 'wt') as soutfile:
			writer = csv.writer(soutfile, delimiter='\t')	
			with open(shuffled_bed) as tsv:
				for window in csv.reader(tsv, dialect="excel-tab"):
					swin_chrom, swin_start, swin_end, swin_gcbin = window[0], window[1], window[2], window[3]
					if wincount_dict[swin_gcbin] != 0:
						writer.writerow([swin_chrom, swin_start, swin_end, swin_gcbin])
						wincount_dict[swin_gcbin] -= 1
		soutfile.close()
	# Count the number of alignments in each gc bin, add to dict, and compute statistics on read counts
		with open(selected_bed) as tsv:
			for window in csv.reader(tsv, dialect="excel-tab"):
				count = 0
				if chrom == window[0]:
					win_start = window[1]
					win_end = window[2]
					win_gc = window[3]
					for Aln in AlnFile.bam.fetch(reference = chrom, start = int(win_start), end = int(win_end)):
						if Aln.is_unmapped: continue
						midpoint = ((Aln.reference_end - Aln.reference_start) / 2) + Aln.reference_start
						if int(win_start) <= midpoint <= int(win_end):
							count += 1
				aln_count[win_gc].append(count)
		gc_metadata[win_size] = {}
		for binid in aln_count:
			mean_count = np.mean(aln_count[binid])
			std_count = np.std(aln_count[binid])
			gc_metadata[win_size][binid] = (mean_count, std_count)

	# Adding chromosomal DOC, MPD, MAD data and gc bin data to a tuple to output one object as json file 
	metadata = (chr_metadata, gc_metadata)

	# Export chrom_metadata and gc_metadata as JSON file
	meta_json = '{}.metadata.chonk.json'.format(Args.o)
	with open(meta_json, 'w') as json_outfile:
		json.dump(metadata, json_outfile)
		json_outfile.close()

	# Deleting files that were made, except for metadata json output file
	chonkdel_cwd = "{}/{}.{}.CHONKdel*".format(cwd, genome, chrom)
	rm_cmd = "rm {}".format(chonkdel_cwd)
	sp.check_output(rm_cmd, shell = True)