import csv, os, subprocess as sp
from collections import defaultdict

def split_bed(svtype_file):

	# read file
	svtype_file_fh = open(svtype_file)

	# write to file
	regioned_file = svtype_file[:svtype_file.find('.txt')] + '.TMPregioned.txt'
	regioned_fh = open(regioned_file, 'w', newline='')
	regioned_writer = csv.writer(regioned_fh, delimiter='\t')

	for sv in csv.reader(svtype_file_fh, dialect='excel-tab'):
		if sv[0] == 'chrom': continue
		chrom, start, end, svtype, ci_start, ci_end = sv[0], int(sv[1]), int(sv[2]), sv[3], sv[6], sv[7]

		cistart_start = int(ci_start[:ci_start.find(',')])
		cistart_end = int(ci_start[ci_start.find(',')+1:])
		ciend_start = int(ci_end[:ci_end.find(',')])
		ciend_end = int(ci_end[ci_end.find(',')+1:])

		sv_start, sv_end = start, end
		lf_start, lf_end = start-1000, start 
		rf_start, rf_end = end, end+1000
		lo_start, lo_end = (start + cistart_start), (start + cistart_end)
		ro_start, ro_end = (end + ciend_start), (end + ciend_end)

		# sv body
		regioned_writer.writerow( [chrom, sv_start, sv_end, svtype, 'sv'] )
		# left flank
		regioned_writer.writerow( [chrom, lf_start, lf_end, svtype, 'lf'] )
		# right flank
		regioned_writer.writerow( [chrom, rf_start, rf_end, svtype, 'rf'] )
		# left overlap
		regioned_writer.writerow( [chrom, lo_start, lo_end, svtype, 'lo'] )
		# right overlap
		regioned_writer.writerow( [chrom, ro_start, ro_end, svtype, 'ro'] )

	regioned_fh.close()

	return regioned_file

def add_features(feat_file, svtype, overlap_dict):
	feat_ofile = feat_file[:feat_file.find('.{}.'.format(svtype))+4] + '.final.txt'
	feat_ofile_fh = open(feat_ofile, 'w', newline='')
	feat_writer = csv.writer(feat_ofile_fh, delimiter='\t')

	svtype_file_fh = open(feat_file)

	if svtype == 'DEL' or svtype == 'DUP':
		header = ['chrom', 'start', 'end', 'svtype', 'iid', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
			 'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std', 'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std', 'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
			 'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
			 'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
			 'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
			 'll', 'lr', 'la', 'rl', 'rr', 'ra', 'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
			 'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc', 'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp', 'log_sv_len', 'bp_start_ci', 'bp_end_ci',
			 'sv_rm', 'sv_sd', 'sv_str', 'lf_rm', 'lf_sd', 'lf_str', 'rf_rm', 'rf_sd', 'rf_str', 'lo_rm', 'lo_sd', 'lo_str', 'ro_rm', 'ro_sd', 'ro_str']
	elif svtype == 'INV':
		header = ['chrom', 'start', 'end', 'svtype', 'id', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
			 'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std', 'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std', 'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
			 'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
			 'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
			 'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
			 'll', 'lr', 'lla', 'lra', 'rl', 'rr', 'rla', 'rra', 'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
			 'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc', 'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp', 'log_sv_len', 'bp_start_ci', 'bp_end_ci',
			 'sv_rm', 'sv_sd', 'sv_str', 'lf_rm', 'lf_sd', 'lf_str', 'rf_rm', 'rf_sd', 'rf_str', 'lo_rm', 'lo_sd', 'lo_str', 'ro_rm', 'ro_sd', 'ro_str']

	feat_writer.writerow( header )

	feat_fh = open(feat_file)
	for sv in csv.reader(svtype_file_fh, dialect='excel-tab'):
		if sv[0] == 'chrom': continue
		chrom, start, end, svtype, ci_start, ci_end = sv[0], sv[1], sv[2], sv[3], sv[6], sv[7]

		lf_start, lf_end = str( int(start)-1000 ), start
		rf_start, rf_end = end, str( int(end)+1000 )
		lo_start, lo_end = str(int(start) + int(float(ci_start.split(',')[0]))), str(int(start) + int(float(ci_start.split(',')[1])))
		ro_start, ro_end = str(int(end) + int(float(ci_end.split(',')[0]))), str(int(end) + int(float(ci_end.split(',')[1])))

		sv_len, flank_len = int(end)-int(start), 1000
		left_overlap_len = int(float(ci_start.split(',')[1])) + -int(float(ci_start.split(',')[0]))
		right_overlap_len = int(float(ci_end.split(',')[1])) + -int(float(ci_end.split(',')[0]))

		sv_rm =  overlap_dict['rm'][(chrom, start, end, svtype, 'sv')] / sv_len
		sv_sd =  overlap_dict['sd'][(chrom, start, end, svtype, 'sv')] / sv_len
		sv_str = overlap_dict['str'][(chrom, start, end, svtype, 'sv')] / sv_len
		lf_rm =  overlap_dict['rm'][(chrom, lf_start, lf_end, svtype, 'lf')] / flank_len
		lf_sd =  overlap_dict['sd'][(chrom, lf_start, lf_end, svtype, 'lf')] / flank_len
		lf_str = overlap_dict['str'][(chrom, lf_start, lf_end, svtype, 'lf')] / flank_len
		rf_rm =  overlap_dict['rm'][(chrom, rf_start, rf_end, svtype, 'rf')] / flank_len
		rf_sd =  overlap_dict['sd'][(chrom, rf_start, rf_end, svtype, 'rf')] / flank_len
		rf_str = overlap_dict['str'][(chrom, rf_start, rf_end, svtype, 'rf')] / flank_len
		lo_rm =  overlap_dict['rm'][(chrom, lo_start, lo_end, svtype, 'lo')] / left_overlap_len
		lo_sd =  overlap_dict['sd'][(chrom, lo_start, lo_end, svtype, 'lo')] / left_overlap_len
		lo_str = overlap_dict['str'][(chrom, lo_start, lo_end, svtype, 'lo')] / left_overlap_len
		ro_rm =  overlap_dict['rm'][(chrom, ro_start, ro_end, svtype, 'ro')] / right_overlap_len
		ro_sd =  overlap_dict['sd'][(chrom, ro_start, ro_end, svtype, 'ro')] / right_overlap_len
		ro_str = overlap_dict['str'][(chrom, ro_start, ro_end, svtype, 'ro')] / right_overlap_len

		sv.extend( [sv_rm, sv_sd, sv_str, lf_rm, lf_sd, lf_str, rf_rm, rf_sd, rf_str, lo_rm, lo_sd, lo_str, ro_rm, ro_sd, ro_str] )

		feat_writer.writerow(sv)

	feat_ofile_fh.close()

	#remove command for incomplete feat_file
	rm_feat_cmd = "rm {}".format(feat_file)
	sp.call(rm_feat_cmd, shell=True)

##############################################################################
# MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN #
##############################################################################

def overlap_features(del_file, dup_file, inv_file, rm_merged, sd_merged, str_merged):

	# split the bed files into regions (sv, lf, rf, lo, ro)
	del_regioned = split_bed(del_file)
	dup_regioned = split_bed(dup_file)
	inv_regioned = split_bed(inv_file)
	file_name = del_regioned[del_regioned.rfind('/')+1:]
	uniq_id = file_name[:file_name.find('.')]

	for file in del_regioned, dup_regioned, inv_regioned:

		# sort the bedfile for quicker intersect
		sorted_file = file[:file.find('.txt')] + '.sorted.txt'
		sort_cmd = "sort {} | uniq | sortBed > {}".format(file, sorted_file)
		sp.call(sort_cmd, shell=True)

		# remove previous unsorted file	
		rm_cmd = 'rm {}'.format(file)
		sp.call(rm_cmd, shell=True)

		for track, tname in zip((rm_merged, sd_merged, str_merged), ('.rm.txt', '.sd.txt', '.str.txt')):
			ofile = sorted_file[:sorted_file.find('.txt')] + tname

			# use the intersectBed command to get the intersection of the rm, sd, str tracks on the sorted regioned files
			intersect_cmd = 'intersectBed -a {} -b {} -wao -sorted > {}'.format(sorted_file, track, ofile)
			sp.call(intersect_cmd, shell=True)

	# create a dictionary with the number of overlaps (rm, sd, str) for all svs
	overlap_dict = {'rm':defaultdict(int), 'sd':defaultdict(int), 'str':defaultdict(int)}
	directory = os.path.dirname(del_file)
	for filename in os.listdir(directory):
		if (filename.endswith('.rm.txt') or filename.endswith('.sd.txt') or filename.endswith('.str.txt')) and uniq_id in filename:
			filename_abspath = directory +'/'+ filename
			fh = open(filename_abspath)
			trname = filename[filename.find('.TMPregioned.sorted.')+20:filename.find('.txt')]
			for sv in csv.reader(fh, dialect='excel-tab'):
				chrom, start, end, svtype, region, n_overlap = sv[0], sv[1], sv[2], sv[3], sv[4], int(sv[8])
				sv_key = (chrom, start, end, svtype, region)
				overlap_dict[trname][sv_key]+=n_overlap

	rm_regioned_cmd = "rm {}*.TMPregioned.*".format(directory+'/')
	sp.call(rm_regioned_cmd, shell=True)

	# add the overlap featues to the existing del, dup, and inv features
	add_features(del_file, 'DEL', overlap_dict)
	add_features(dup_file, 'DUP', overlap_dict)
	add_features(inv_file, 'INV', overlap_dict)