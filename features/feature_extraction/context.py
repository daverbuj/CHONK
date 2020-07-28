from math import log10
import numpy as np
import sys
import zlib

###############################################################

def get_gc_perc(sequence):
	gc_count = 0
	region_size = 0

	# count up the number of g's and c's in the sequence
	for base in sequence:
		region_size += 1
		if base.upper() == 'G' or base.upper() == 'C':
			gc_count += 1
	try:
		gc_perc = gc_count / region_size
	except ZeroDivisionError:
		print(sequence)
		gc_perc = np.nan

	return gc_perc

def sv_gc_content(sv_seq, lf_seq, rf_seq, lo_seq, ro_seq):
	sv_gc = get_gc_perc(sv_seq)
	lf_gc = get_gc_perc(lf_seq)
	rf_gc = get_gc_perc(rf_seq)
	lo_gc = get_gc_perc(lo_seq)
	ro_gc = get_gc_perc(ro_seq)

	return sv_gc, lf_gc, rf_gc, lo_gc, ro_gc

####################################################################

def complexity(seq):
    seq = seq.upper()
    comp_sz = sys.getsizeof(zlib.compress(seq.encode('utf-8')))-sys.getsizeof(zlib.compress(''.encode('utf-8')))
    uncomp_sz = sys.getsizeof(seq.encode('utf-8'))-sys.getsizeof(''.encode('utf-8'))
    if uncomp_sz > 0: return comp_sz/uncomp_sz
    else: return 'NaN'

def sv_complexity(sv_seq, lf_seq, rf_seq, lo_seq, ro_seq):
	sv_comp = complexity(sv_seq)
	lf_comp = complexity(lf_seq)
	rf_comp = complexity(rf_seq)
	lo_comp = complexity(lo_seq)
	ro_comp = complexity(ro_seq)

	return sv_comp, lf_comp, rf_comp, lo_comp, ro_comp

###################################################################

def sv_length(start, end):

	log_sv_len = log10( end-start )
	return log_sv_len

########################################################################

def bp_ci(start, end, ci_start, ci_end, rlen):

	# parse start and end positions from ci_start
	cistart_start = int(ci_start[:ci_start.find(',')])
	cistart_end = int(ci_start[ci_start.find(',')+1:])

	# parse start and end positions from ci_end
	ciend_start = int(ci_end[:ci_end.find(',')])
	ciend_end = int(ci_end[ci_end.find(',')+1:])

	# calculate breakpoint start and end confidence intervals
	bp_start_ci = (cistart_end - cistart_start) / rlen
	bp_end_ci = (ciend_end - ciend_start) / rlen

	return bp_start_ci, bp_end_ci, cistart_start, cistart_end, ciend_start, ciend_end

#####################################################################################
# MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN #
#####################################################################################

def cont_features(start, end, sequence, seq_start, rlen, ci_start, ci_end, lo_seq, ro_seq, sv_seq, lf_seq, rf_seq, rlen_seq):

	# extract gc content features
	sv_gc, lf_gc, rf_gc, lo_gc, ro_gc = sv_gc_content(sv_seq, lf_seq, rf_seq, lo_seq, ro_seq)

	# extract complexity features
	sv_comp, lf_comp, rf_comp, lo_comp, ro_comp = sv_complexity(sv_seq, lf_seq, rf_seq, lo_seq, ro_seq)

	# extract sv length feature
	log_sv_len = sv_length(start, end)

	# extract breakpoint confidence interval features
	bp_start_ci, bp_end_ci, cistart_start, cistart_end, ciend_start, ciend_end = bp_ci(start, end, ci_start, ci_end, rlen)

	context_features = (sv_gc, lf_gc, rf_gc, lo_gc, ro_gc, sv_comp, lf_comp, rf_comp, lo_comp, ro_comp, log_sv_len, bp_start_ci, bp_end_ci)

	return context_features