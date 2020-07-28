import argparse, os, csv, sys, datetime, subprocess as sp
from chonk.Bam import Bam
from chonk.Metadata import Metadata
from coverage import cov_features
from supp_frag import sf_features
from kmer_junction import kmer_features
from context import cont_features
from overlap import overlap_features

############################################## INPUT FILE SHOULD LOOK LIKE THIS
# CHROM 	START 	END 	SVTYPE 	CIPOS 	CIEND 	IID 	GT

############################################## ARGUMENTS

parser = argparse.ArgumentParser(description='Extract kmer junction features.')
parser.add_argument('-bp'      ,type=str,   required=True, default=None, help="Breakpoint file")
parser.add_argument('-k'       ,type=int,   required=True, default=None, help="Kmers size for kmer junction")
parser.add_argument('-rlen'    ,type=float, required=True, default=None, help="Read length modifier for kmer junction")
parser.add_argument('-pk'      ,type=int,   required=True, default=None, help="Kmers size for pseudoalignment")
parser.add_argument('-o'       ,type=str,   required=True, default=None, help="Output file w/ directory")
parser.add_argument('-metadir' ,type=str,   required=True, default=None, help="metadata directory")
parser.add_argument('-ref' 	   ,type=int,   required=True, default=None, help="reference genome")
parser.add_argument('-v'       , '--verbose'  , help="increase output verbosity"                     , action='store_true')
parser.add_argument('-t'       , '--timestamp', help="add timestamp to debug messages"               , action='store_true')
args = parser.parse_args()

# Inputs
meta_dir = args.metadir
if not meta_dir.endswith('/'):
	meta_dir+='/'
bp_file = args.bp
ofile = args.o
if not ofile.endswith('.txt'):
	sys.stderr.write("FATAL ERROR: -o argument must end in '.txt'\n")
	sys.exit(1)
k = args.k
rlen_modifier = args.rlen
pk = args.pk
if args.ref == 19:
	fasta_file = '/home/daverbuj/resources/fastas/human_g1k_v37_decoy.fasta'
	merged_rm = '/home/daverbuj/resources/repeat_tracks/hg19_repeatmasker_merged.bed'
	merged_sd = '/home/daverbuj/resources/repeat_tracks/hg19_segdup_merged.bed'
	merged_str = '/home/daverbuj/resources/repeat_tracks/hg19_str_merged.bed'
elif args.ref == 38:
	fasta_file = '/home/daverbuj/resources/fastas/GRCh38_full_analysis_set_plus_decoy_hla.fa'
	merged_rm = '/home/daverbuj/resources/repeat_tracks/hg38_repeatmasker_merged.bed'
	merged_sd = '/home/daverbuj/resources/repeat_tracks/hg38_segdup_merged.bed'
	merged_str = '/home/daverbuj/resources/repeat_tracks/hg38_str_merged.bed'

##################################################

def report(message):
	time = '{:%H:%M:%S}'.format(datetime.datetime.now())
	stopwatch = time 
	if args.verbose:
		if args.timestamp:
			print('{} DAN_MSG: '.format(time), message)
		else:
			print('DAN_MSG: ', message)

def create_meta_dict(meta_dir):
	meta_dict = {}
	for filename in os.listdir(meta_dir):
		if filename.endswith("chonk.json"):
			meta_file = meta_dir + filename
			Meta = Metadata()
			Meta.read_json(meta_file)
			# sample needs to match the value in the 7th column (genome)
			sample = filename[:filename.find('_')] # USED WHEN FID.IID IN BP FILE # PTCD and HGSV and GIAB
			#sample = Meta.bam_path[Meta.bam_path.find('trios/')+6 : Meta.bam_path.find('.polaris.bam')] # USED WHEN LONG AWKWARD SAMPLE NAME IN BP FILE # USED FOR POLARIS
			meta_dict[sample] = Meta
	return meta_dict

def new_sv_bed(bp_file, meta_dict, rlen_modifier, ofile, fasta_file):

	## create new bed file which includes flanking region in start and end, modified flank size derived from read length, and fasta sequence
	# open original bed file
	sv_bed_fh = open(bp_file)

	# create new bed file with flanking region in start and end, and modified read length
	odir = os.path.dirname(bp_file)
	file_name = ofile[ofile.rfind('/')+1:]
	uniq_id = file_name[:file_name.find('.')] # used when splitting input file into multiple files. use a "." to differentiate between files in the output(ex. output = sample1.output.txt)
	new_bed_ofile = odir + '/{}.flanked_bpfile.TEMP.bed'.format(uniq_id)
	new_bed_fh = open(new_bed_ofile, 'w', newline='')
	new_bed_writer = csv.writer(new_bed_fh, delimiter='\t')
	for sv in csv.reader(sv_bed_fh, dialect='excel-tab'):
		chrom, start, end, svtype, ci_start, ci_end, genome, genotype = sv
		start, end = int(start), int(end)
		#sample = genome.replace('_','.') # POLARIS
		sample = genome.split('.')[1] # PTCD and HGSV
		#sample = genome $ used for somatic when == chrxx
		mod_rlen = int(float(meta_dict[sample].read_len[chrom]) * float(rlen_modifier))
		new_start = start - 1000
		new_end = end + 1000

		# write to file
		new_bed_writer.writerow( [chrom, new_start, new_end, svtype, ci_start, ci_end, sample, genotype, mod_rlen] )
	sv_bed_fh.close()
	new_bed_fh.close()

	# append the fasta sequence to the end of each row using 'bedtools getfasta'
	sv_fasta_bed = odir + '/{}.fasta_bpfile.TEMP.bed'.format(uniq_id)
	fasta_cmd = 'bedtools getfasta -bedOut -fi {} -bed {} > {}'.format(fasta_file, new_bed_ofile, sv_fasta_bed)
	sp.call(fasta_cmd, shell = True)

	# delete flanked bed file and unsorted fasta file
	del_cmd = 'rm {}'.format(new_bed_ofile)
	sp.call(del_cmd, shell = True)

	return sv_fasta_bed

def check_ci(ci):
	if ci == '.' or ci == '0,0': return '-20,20'
	start_value = int(ci[:ci.find(',')])
	end_value = int(ci[ci.find(',')+1:])
	if start_value > 0:
		start_value *= -1
	if end_value < 0:
		end_value *= -1
	new_ci = '{},{}'.format(start_value, end_value)
	return new_ci

def parse_seq(sequence, start, end, seq_zero, ci_start, ci_end, mod_rlen):
	# Get SV body sequence 
	sv_start = start - seq_zero
	sv_end = end - seq_zero +1
	sv_seq = sequence[sv_start:sv_end]

	# Get flanking sequence
	lf_start = sv_start - mod_rlen
	rf_end = sv_end + mod_rlen
	rlen_seq = sequence[lf_start:rf_end]

	# left/right flank seq
	lf_seq = sequence[lf_start:sv_start]
	rf_seq = sequence[sv_end:rf_end]

	# Get left overlapping sequence
	lo_start = int(sv_start) + int(float(ci_start.split(',')[0]))
	lo_end = int(sv_start) + int(float(ci_start.split(',')[1]))
	lo_seq = sequence[lo_start:lo_end]

	# Get left/right 1000 seq
	l1k = sequence[:1000]
	r1k = sequence[-1000:]

	# Get right overlapping sequence
	ro_start = int(sv_end) + int(float(ci_end.split(',')[0]))
	ro_end = int(sv_end) + int(float(ci_end.split(',')[1]))
	ro_seq = sequence[ro_start:ro_end]

	return sv_seq, rlen_seq, lo_seq, l1k, ro_seq, r1k, lf_seq, rf_seq

##########################################################################################################
# MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN #
##########################################################################################################

report('Starting extraction script')

# Create 3 outfiles for Dels, Dups, and Invs
# DEL outfile
del_ofile = ofile[:ofile.rfind('.')] + '.DEL' + ofile[ofile.rfind('.'):]
del_ofile_fh = open(del_ofile, 'w', newline='')
del_ofile_writer = csv.writer(del_ofile_fh, delimiter='\t')

# DUP outfile
dup_ofile = ofile[:ofile.rfind('.')] + '.DUP' + ofile[ofile.rfind('.'):]
dup_ofile_fh = open(dup_ofile, 'w', newline='')
dup_ofile_writer = csv.writer(dup_ofile_fh, delimiter='\t')

# INV outfile
inv_ofile = ofile[:ofile.rfind('.')] + '.INV' + ofile[ofile.rfind('.'):]
inv_ofile_fh = open(inv_ofile, 'w', newline='')
inv_ofile_writer = csv.writer(inv_ofile_fh, delimiter='\t')

# Create outfile headers
# Dels and Dups
report("Writing headers to the output files")
deldup_header = ['chrom', 'start', 'end', 'svtype', 'iid', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
				 'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std', 'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std', 'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
				 'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
				 'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
				 'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
				 'll', 'lr', 'la', 'rl', 'rr', 'ra', 'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
				 'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc', 'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp', 'log_sv_len', 'bp_start_ci', 'bp_end_ci']
# INVs
inv_header = ['chrom', 'start', 'end', 'svtype', 'id', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
				 'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std', 'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std', 'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
				 'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
				 'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
				 'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
				 'll', 'lr', 'lla', 'lra', 'rl', 'rr', 'rla', 'rra', 'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
				 'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc', 'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp', 'log_sv_len', 'bp_start_ci', 'bp_end_ci']

del_ofile_writer.writerow( deldup_header )
dup_ofile_writer.writerow( deldup_header )
inv_ofile_writer.writerow( inv_header )

# create metadata dictionary for all samples in the metadata directory
report("Creating metadata dictionary")
meta_dict = create_meta_dict(meta_dir)

# Create a new bed file which now includes the left and right flanks for start and end, the new modified rlen, and the fasta sequence

sv_fasta_bed = new_sv_bed(bp_file, meta_dict, rlen_modifier, ofile, fasta_file)
sv_fasta_bed_fh = open(sv_fasta_bed)

csv.field_size_limit(int(sys.maxsize/100))
report('Starting first extraction')
count = 0
for sv in csv.reader(sv_fasta_bed_fh, dialect='excel-tab'):
	if int(flanked_end) - int(flanked_start) > 5e6: continue
	count+=1
	report('~~~~~~~~~~~~~~~ SV {}'.format(count))

	# SV info from row
	try:
		chrom, flanked_start, flanked_end, svtype, cistart, ciend, sample, genotype, mod_rlen, sequence = sv
	except ValueError:
		print(sv)
		continue
	start, end = int(flanked_start)+1000, int(flanked_end)-1000
	ci_start = check_ci(cistart)
	ci_end = check_ci(ciend)
	rlen = int(int(mod_rlen)/float(rlen_modifier))
	sv_seq, rlen_seq, lo_seq, l1k, ro_seq, r1k, lf_seq, rf_seq = parse_seq(sequence, start, end, int(flanked_start), ci_start, ci_end, int(mod_rlen))
	# get metadata file for SV
	Meta = meta_dict[sample]
	# check if alignment file exists
	AlnFile = Bam(Meta.bam_path)

	# extract coverage features
	report('Cov Start')
	coverage_features = cov_features(AlnFile, Meta, chrom, start, end, sv_seq, lf_seq, rf_seq)

	# extract supporting fragment features
	report('SF Start')	
	suppfrag_features = sf_features(AlnFile, Meta, chrom, start, end, ci_start, ci_end, svtype)

	# extract kmer junction features
	report('KJ Start')	
	km_features = kmer_features(chrom, start, end, AlnFile, rlen_seq, svtype, int(mod_rlen), int(k), int(pk), ci_start, ci_end, lo_seq, ro_seq)

	# extract context-based features
	report('Cont Start')
	context_features = cont_features(start, end, sequence, flanked_start, rlen, ci_start, ci_end, lo_seq, ro_seq, sv_seq, lf_seq, rf_seq, rlen_seq)

	## start to output features to respective outfile
	# determine genotype label
	if '1' not in genotype: gt = 0
	elif '0' not in genotype: gt = 2
	elif '1' in genotype and '0' in genotype: gt = 1
	else: gt = '.' # place holder for genotype if we dont have genotypes at the moment

	base_output = (chrom, start, end, svtype, sample, gt, ci_start, ci_end, k, rlen_modifier)
	output_data = base_output + coverage_features + suppfrag_features + km_features + context_features

	if svtype == 'DEL':
		del_ofile_writer.writerow( output_data )
	elif svtype == 'DUP':
		dup_ofile_writer.writerow( output_data )
	elif svtype == 'INV':
		inv_ofile_writer.writerow( output_data )

del_ofile_fh.close()
dup_ofile_fh.close()
inv_ofile_fh.close()
sv_fasta_bed_fh.close()

# delete fasta bed file
del_cmd = 'rm {}'.format(sv_fasta_bed)
sp.call(del_cmd, shell = True)

overlap_features(del_ofile, dup_ofile, inv_ofile, merged_rm, merged_sd, merged_str)
report('DONE')