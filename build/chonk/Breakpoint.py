#!/usr/bin/env python3
from chonk.Bam import Bam
from chonk.Metadata import Metadata
import chonk.Backend as Backend
import chonk.Splitread as Splitread
import chonk.DiscordantPE as DiscordantPE
import pysam, json, os, sys, csv


def breakpoints(Args):
	"""
	main function for SV breakpoint detection
	"""
	
	# metadata file
	Meta = Metadata()
	Meta.read_json(Args.m)

	# check if alignment file exists
	AlnFile = Bam(Meta.bam_path)

	# check that contigs selected match those in meta file
	# if no contigs selected, use the contigs available in the meta file
	contigs = selected_contigs(Args.r, Meta.user_contigs, AlnFile)

	# output bed file
	output_file = Args.o + Meta.sample + '_' + '_'.join(contigs) + '_breakpoints.chonk.bed'
	if contigs == Meta.contigs:
		output_file = Args.o + Meta.sample + '_breakpoints.chonk.bed'
	bed_fh = open(output_file,'wt')
	bed_writer = csv.writer(bed_fh,delimiter='\t')
	bed_writer.writerow( ['#chrom', 'start', 'end', 'svtype'] )

	## call sv breakpoints using split reads and discordant paired ends
	for chrom in contigs:
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
					bed_writer.writerow(sv)

				processed_reads.append(read_name)

			# discordant paired end
			if len(breaks) == 0:
				sv = DiscordantPE.discordantpe(Aln,Meta,chrom)
				if sv:
					bed_writer.writerow(sv)

	bed_fh.close()

	Backend.reporter('breakpoints complete.\n    bed output ---> {}'.format(os.path.abspath(Args.o)))

##########################

# determines contigs to use based on user input and what contigs are available to analyze from the previous subcommand
def selected_contigs(selected_contigs, available_contigs, AlnFile):
	contigs = ()

	# ensure user-defined contigs are found in file
	if selected_contigs != None:
		for c in selected_contigs:
			# ensure contig name matches those of bam file
			c = Backend.check_chrom(c, AlnFile.chrom_flag)
			if c not in available_contigs:
				# the contig selected is not available in the metadata file. Run chonk metadata again with selected contigs
				sys.stderr.write('FATAL ERROR: contig {} not extracted to metadata file\n'.format(c))
				sys.exit(1)
			else:
				contigs += (c,)
	else: 
		# if the user did not specify contigs, use the contigs available in the meta file
		contigs = available_contigs
		
	if len(contigs)==0:
		sys.stderr.write('FATAL ERROR: all user defined contigs were not found in the bam header\n')
		sys.exit(1)

	return contigs
