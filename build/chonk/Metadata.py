#!/usr/bin/env python3
from chonk.Bam import Bam
from chonk.Backend import alignment_midpoint, check_chrom, get_gc_bin, reporter,tuple2string, Welford

import csv, json, os, pysam, random, sys
import subprocess as sp
from operator import itemgetter


#############
class Metadata(object):
	def __init__(self):
		
		# metadata variables
		self.window_lengths = (25, 1000)
		self.n_windows = int(1e4)
		self.n_reads = int(5e6)
		self.tlen_cap = int(1e4)

		# genomic regions to extract metadata
		self.binned_loci = {}  # [ contig ] -> [ (start, end, gc_bin, window_length) ]
		self.loci = {}  # [ contig ] -> [ (start, end) ]

		self.bam_path=None
		self.sample=None
		self.outdir=None
		self.tmpdir=None
		self.fasta=None
		self.contigs=None
		self.user_contigs=None
		self.json_file=None

		# contig depth of coverage
		self.doc = {} # [ contig ] = DOC

		# contig mean read length
		self.read_len = {} # [ contig ] = mean read length

		# contig mean template length
		self.tlen = {} # [ contig ] = mean template length

		# contig template length standard deviation
		self.tlen_std = {} # [ contig ] = template length std

		# contig gc binned mean read counts
		self.gc_rc = {} # [ (contig, window length, GC bin ) ] = mean read counts

		# contig gc binned read count standard deviation
		self.gc_std = {} # [ (contig, window length, GC bin ) ] = read count std

	def __setattr__(self,name,value):
		super().__setattr__(name,value)
	
	def init_meta(self,Bam,Args):
		# this function is used when the Metadata instance is 
		# initalized for the first time

		# BAM file path
		self.bam_path = Bam.bam_path
		# sample name
		self.sample=Bam.sample
		
		# output directory
		self.outdir = Args.o
	
		# FASTA file
		self.fasta=os.path.abspath(Args.f)

		# contigs as a tuple 
		self.contigs = tuple(Bam.bam.references)
		self.user_contigs = self.init_contigs(Args.r,Bam.chrom_flag)

		# temporary directory
		tmpdir = Args.o + self.sample +'_tmp' + '/'
		if self.contigs != self.user_contigs:
			tmpdir = Args.o + self.sample + '_' + '_'.join(self.user_contigs) + '_tmp' + '/'	
		if not os.path.exists(tmpdir):
			os.mkdir(tmpdir)
		self.tmpdir=tmpdir

		# init loci dict values as empty lists
		for x in self.user_contigs: 
			self.binned_loci[x]=[]
			self.loci[x]=[]
		
		#output json file
		self.json_file = self.outdir + self.sample + '_' + '_'.join(self.user_contigs) + '_metadata.chonk.json' 
		if self.user_contigs == self.contigs:
			self.json_file = self.outdir + self.sample + '_metadata.chonk.json'

	def init_contigs(self,contigs,chrom_flag):
		# returns a list of contigs to iterate over

		user_contigs = []
		
		if contigs != None:
			for c in contigs:
				# ensures user defined contig matches contig name in Bam file
				c = check_chrom(c,chrom_flag)
				if c not in self.contigs: 
					sys.stderr.write('WARNING: {} not found in the bam header\n'.format(c))
				else:
					user_contigs.append(c)
		else: 
			# if the user did not specify contigs, use the entire genome
			user_contigs = self.contigs
	
		if len(user_contigs)==0:
			sys.stderr.write('FATAL ERROR: all user defined contigs were not found in the bam header\n')
			sys.exit(1)

		return tuple(user_contigs)

	def output_json(self):
		# save metadata object as a json file
		#purge loci and binned_loci
		self.loci = {}
		self.binned_loci = {}

		out_file_handle = open(self.json_file,'w')
		json_str = json.dump(self.__dict__, out_file_handle, indent=4)
		out_file_handle.close()


		# remove tmpdir
		#sp.call("rm -rf "+self.tmpdir, shell = True)

	def read_json(self,json_file):
		# load metadata object from a json file
		with open(json_file) as json_fh:
			data = json.load(json_fh)
			for attr in data:
				self.__setattr__(attr,data[attr])

#############
def metadata(Args):
	"""
	main function for genome metadata extraction
	"""

	# init the Bam object
	AlnFile = Bam(Args.i)

	# init the Meta object
	Meta = Metadata()
	Meta.init_meta(AlnFile,Args)

	# define and filter the mask according to 
	mask = init_mask(AlnFile,Args.x, Meta)

	# two gc binned bed files
	gc_binned_beds = bin_gc_content( window_mask(mask,Args.f,Meta.tmpdir,Meta.window_lengths) )

	# load the metadata loci
	load_loci(Meta,gc_binned_beds)

	compute_contig_metadata(Meta,AlnFile)
	compute_gc_binned_metadata(Meta,AlnFile)
	
	# save to a json file
	Meta.output_json()

	reporter('metadata extraction complete.\n    JSON output ---> {}'.format(os.path.abspath(Meta.json_file)))

	#New_Meta = Metadata()
	#New_Meta.read_json(Meta.json_file)

########
def init_mask(Bam,exclude,Meta):
	# returns a bed file to extract metadata from
	
	mask=Meta.tmpdir + "masked.bed"
	# make bed file of just contigs and lengths
	contig_bed = Meta.tmpdir + "contigs.bed"
	contig_bed_fh = open(contig_bed,'wt')
	contig_bed_writer = csv.writer(contig_bed_fh, delimiter = '\t', lineterminator='\n')
	for i,e in enumerate(Bam.bam.references):
		# only output contigs the user defined
		if e in Meta.user_contigs:
			contig_bed_writer.writerow( [e,0,Bam.bam.lengths[i]] )
	contig_bed_fh.close()

	# if the user did NOT define an exclude file, use the entire genome
	if exclude==None: mask=contig_bed
	else: 
		
		# if the user defined an exclude file, subtract the mask from the genome
		subtract_cmd = "bedtools subtract -a {} -b {} > {}".format(contig_bed, os.path.abspath(exclude), mask)
		sp.call(subtract_cmd, shell = True)
		filtered_mask = mask.replace('.bed','.gt10kb.bed')
		# * Filter out regions smaller than 10kb
		filter_cmd = """awk 'BEGIN {OFS="\t"}; $3 - $2 >= 10000 {print $1, $2, $3}' """ + mask + ' > ' + filtered_mask
		sp.call(filter_cmd, shell = True)
		mask=filtered_mask

	with open(mask,'r') as file_handle:
		for line in file_handle:
			chrom, start, end = line.rstrip().split('\t')
			Meta.loci[chrom].append( (int(start), int(end)) )
	
	return mask

def window_mask(mask,fasta,tmpdir,window_lengths):
	# returns two bed files windowed into 1000bp and 25bp intervals

	windowed = {} # return this dict containing the paths to the gc_bed files
	# [ window length ] = bed file
	window_bed_pre = tmpdir + 'masked.windowed.'
	
	for wlen in window_lengths:
		window_bed = window_bed_pre + str(wlen) + 'bp.bed'

		# make windows using a 20bp step
		win_cmd = "bedtools makewindows -w {} -s 20 -b {} > {}".format(wlen, mask, window_bed)
		sp.call(win_cmd, shell = True)

		# * Get GC content for each window
		gc_bed = window_bed.replace('.bed','.gc.bed')
		gc_cmd = "bedtools nuc -fi {} -bed {} | cut -f 1-3,5 > {}".format(os.path.abspath(fasta), window_bed, gc_bed)
		sp.call(gc_cmd, shell = True)

		windowed[wlen] = gc_bed

	return windowed

def bin_gc_content(windowed):
	# returns two bed files binned according to GC content

	binned = {} # return this dict containing the paths to the shuffled gc_bed files
	# [ window length ] -> bed file path

	for wlen in windowed:
		gc_bed = windowed[wlen]
		# * Create a bed file with grouped gc bins
		binned_bed = gc_bed.replace('.bed','.binned.bed')
		gbed_file = open(binned_bed, 'wt')
		gbed_writer = csv.writer(gbed_file, delimiter='\t')
		tsv =  open(gc_bed,'r')

		for row in csv.reader(tsv, dialect="excel-tab"):
			gc_bin = get_gc_bin(row[3])
			if gc_bin == -9: continue #skip header 
			gbed_writer.writerow([row[0],row[1], row[2], gc_bin])
		
		tsv.close()
		gbed_file.close()

		# * Shuffle the binned bed file
		line_count = len(open(binned_bed).readlines())
		shuffled_bed = binned_bed.replace('.bed','.shuffled.bed')
		with open(binned_bed) as infile, open(shuffled_bed, 'w') as outfile:
			outfile.writelines(random.sample(infile.readlines(), line_count))
		infile.close()
		outfile.close()

		binned[wlen] = shuffled_bed

	return binned

def load_loci(Meta,binned):
	# take the shuffled binned regions and load the first 10k loci for each GC bin
	
	counts={} # keep track of the finished gc bins
	# [ (chrom, gc_bin, window_length) ] = counts

	for wlen in binned:
		bed_file = binned[wlen]
		with open(bed_file) as file_handle:
			for line in file_handle:
				chrom, start, end, gc_bin = line.rstrip().split('\t')
				start, end, gc_bin = int(start), int(end), int(gc_bin)
				
				# skip loci if we have enough windows
				if counts.get( (chrom, gc_bin, wlen) ) == Meta.n_windows: continue
				
				if counts.get( (chrom, gc_bin, wlen)) == None: counts[(chrom,gc_bin,wlen)]=1
				else: counts[(chrom,gc_bin,wlen)]+=1

				Meta.binned_loci[chrom].append( (start,end,gc_bin,wlen) )

	# for each contig, sort the positions in genome-order
	for contig in Meta.loci:
		Meta.binned_loci[contig].sort(key=itemgetter(0))
		Meta.loci[contig].sort(key=itemgetter(0))

def compute_contig_metadata(Meta,Bam):

	tlen_welford = {}

	# for each contig
	for contig in Meta.user_contigs:
		## reset these variables for each contig
		if tlen_welford.get(contig)==None: tlen_welford[contig] = Welford()
		read_length_sum,read_count, bp_span = 0,0,0
		terminator = False
		##
		for (start,end) in Meta.loci[contig]:
			if terminator==True: break
			########
			for Aln in Bam.bam.fetch(reference=contig,start=start,end=end):

				if Aln.is_unmapped or type(Aln.reference_length) != int: continue
				
				midpoint=alignment_midpoint(Aln)
				
				if start <= midpoint <= end:
					read_count+=1
					read_length_sum+= Aln.reference_length
					
					if abs(Aln.template_length)<= Meta.tlen_cap:
						tlen_welford[contig].update( abs(Aln.template_length) )
					
				if read_count > Meta.n_reads:
					bp_span+= Aln.reference_end - start
					terminator=True
					break

			if terminator==False:
				bp_span += (end-start)

			########		
		if read_count!=0:
			Meta.read_len[contig] = read_length_sum / float(read_count)
		else: 
			Meta.read_len[contig]='nan'
		if bp_span!=0:	
			Meta.doc[contig] = (Meta.read_len[contig] * read_count) / float(bp_span)
		else:
			Meta.doc[contig]='nan'
		Meta.tlen[contig] = tlen_welford[contig].mean
		Meta.tlen_std[contig] = tlen_welford[contig].std

def compute_gc_binned_metadata(Meta,Bam):
	
	gc_welford = {}
	
	# for each contig
	for contig in Meta.user_contigs:
		# for each window 
		for (start,end,gc_bin,wlen) in Meta.binned_loci[contig]:
			
			gc_count=0
			key = tuple2string( ( contig, wlen, gc_bin ) )
			if gc_welford.get(key)==None: gc_welford[key] = Welford()

			# for each alignment
			for Aln in Bam.bam.fetch(reference=contig, start=start, end=end):
				if Aln.reference_start == None or Aln.reference_end==None: continue
				midpoint = alignment_midpoint(Aln)
				# if the mapped midpoint of the alignment is within the window
				if start <= midpoint <= end: 
					gc_count+=1
			gc_welford[key].update(gc_count)

		for key in gc_welford:
			Meta.gc_rc[key] = gc_welford[key].mean
			Meta.gc_std[key] = gc_welford[key].std