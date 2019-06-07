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
	
	# check if alignment file exists
	AlnFile = Bam(Args.i)
	# metadata file
	meta = Args.m
	# chromosome from meta file
	with open(meta) as json_file:
		data = json.load(json_file)
		metadata = data[0]
		chrom = metadata[1]
	# output file
	out = open(Args.o,'w')
	out.write('#chrom\tstart\tend\tsvtype\n')

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
