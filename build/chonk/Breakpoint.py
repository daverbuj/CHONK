#!/usr/bin/env python3
from chonk.Bam import Bam
import chonk.Backend as Backend
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import pysam, json

def breakpoints(Args):
	"""
	main function for SV breakpoint detection
	"""
	
	# metadata file
	meta = Args.m
	# genome
	genome = Args.g
	# check if alignment file exists
	with open(meta) as json_file:
		data = json.load(json_file)
		gen_metadata, chrom_metadata = data[0], data[1]
		meta_bam = gen_metadata[2]
	AlnFile = Bam(meta_bam)
	# ensure the chromosome prefix(es) matches the one(s) in the aln file
	chromos = ()
	for chromosome in Args.r:
		chromos += (Backend.check_chrom(chromosome,AlnFile.chrom_flag),)
	# creating chromosome label to add to file names
	if len(chromos) == 1:
		chrom_label = chromos[0]
	else:
		chrom_label = 'all'
	# output file
	output_file = '{}/{}.{}.breakpoints.chonk.bed'.format(Args.o, genome, chrom_label)
	out = open(output_file,'w')
	out.write('#chrom\tstart\tend\tsvtype\n')

	for chrom in chromos:
		processed_reads=[]
		for Aln in AlnFile.bam.fetch(region=chrom):
			# skip reads with mapping quality of 0
			if Aln.mapq == 0: continue 
			
			# split-read
			s_alns = Splitread.get_split(Aln)
			
			breaks = []
			if len(s_alns)>0:
				# create a unique read name
				read_name=Aln.query_name+'1'
				if Aln.is_read2: read_name=Aln.query_name+'2' 

				# skip if we looked at the read before
				if read_name in processed_reads: 
					del processed_reads[processed_reads.index(read_name)]
					continue

				breaks = Splitread.splitread(Aln,s_alns,chrom)
				for sv in breaks:
					out.write('\t'.join(map(str,sv))+'\n')

				processed_reads.append(read_name)

			# discordant paired end
			if len(breaks) == 0:
				sv = DiscordantPE.discordantpe(Aln,meta,chrom)
				if sv:
					out.write('\t'.join(map(str,sv))+'\n')
			

	out.close()
