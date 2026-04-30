#!/usr/bin/env python3
"""Repeat-element and segmental-duplication overlap feature extraction."""
import csv
import os
import subprocess as sp
from collections import defaultdict

FLANK_SIZE = 1000
_REGION_TYPES = ('sv', 'lf', 'rf', 'lo', 'ro')


def _sv_regions(sv_row):
    """Yield (chrom, start, end, svtype, region_label) for all regions of one SV."""
    chrom = sv_row[0]
    start, end = int(sv_row[1]), int(sv_row[2])
    svtype = sv_row[3]
    ci_start, ci_end = sv_row[6], sv_row[7]

    def _parse_ci(ci_str, pos):
        lo, hi = ci_str.split(',')
        return int(pos) + int(float(lo)), int(pos) + int(float(hi))

    lo_s, lo_e = _parse_ci(ci_start, start)
    ro_s, ro_e = _parse_ci(ci_end, end)

    yield (chrom, start, end, svtype, 'sv')
    yield (chrom, start - FLANK_SIZE, start, svtype, 'lf')
    yield (chrom, end, end + FLANK_SIZE, svtype, 'rf')
    yield (chrom, lo_s, lo_e, svtype, 'lo')
    yield (chrom, ro_s, ro_e, svtype, 'ro')


def _write_regioned_bed(svtype_file):
    """Expand each SV row into five region rows and write to a temporary BED file."""
    regioned = svtype_file.replace('.txt', '.TMPregioned.txt')
    with open(svtype_file) as fh, open(regioned, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        for row in csv.reader(fh, dialect='excel-tab'):
            if row[0] == 'chrom':
                continue
            for region in _sv_regions(row):
                writer.writerow(region)
    return regioned


def _intersect_tracks(sorted_bed, tracks, outdir, uniq_id):
    """Run bedtools intersect against each repeat track and return result paths."""
    results = {}
    for track, tname in zip(tracks, ('.rm.txt', '.sd.txt', '.str.txt')):
        ofile = sorted_bed.replace('.txt', tname)
        sp.call(f'intersectBed -a {sorted_bed} -b {track} -wao -sorted > {ofile}', shell=True)
        results[tname.strip('.')] = ofile
    return results


def _build_overlap_dict(directory, uniq_id):
    """Parse intersect output files into a nested default-dict of overlap lengths."""
    overlap = {'rm': defaultdict(int), 'sd': defaultdict(int), 'str': defaultdict(int)}
    for fname in os.listdir(directory):
        if uniq_id not in fname:
            continue
        for key in overlap:
            if fname.endswith(f'.{key}.txt'):
                with open(os.path.join(directory, fname)) as fh:
                    for row in csv.reader(fh, dialect='excel-tab'):
                        sv_key = (row[0], row[1], row[2], row[3], row[4])
                        overlap[key][sv_key] += int(row[8])
    return overlap


def _region_length(sv_row, region_label):
    """Return the expected length for a given region label."""
    start, end = int(sv_row[1]), int(sv_row[2])
    ci_start, ci_end = sv_row[6], sv_row[7]

    def _ci_len(ci_str):
        lo, hi = ci_str.split(',')
        return int(float(hi)) - int(float(lo))

    lengths = {
        'sv': end - start,
        'lf': FLANK_SIZE,
        'rf': FLANK_SIZE,
        'lo': _ci_len(ci_start),
        'ro': _ci_len(ci_end),
    }
    return lengths[region_label]


def _add_overlap_features(feat_file, svtype, overlap_dict, header):
    """Append overlap columns to *feat_file*, writing a new *.final.txt file."""
    final = feat_file.replace(f'.{svtype}.', '.') + '.final.txt'
    with open(feat_file) as fh_in, open(final, 'w', newline='') as fh_out:
        writer = csv.writer(fh_out, delimiter='\t')
        writer.writerow(header + [
            'sv_rm', 'sv_sd', 'sv_str',
            'lf_rm', 'lf_sd', 'lf_str',
            'rf_rm', 'rf_sd', 'rf_str',
            'lo_rm', 'lo_sd', 'lo_str',
            'ro_rm', 'ro_sd', 'ro_str',
        ])
        for row in csv.reader(fh_in, dialect='excel-tab'):
            if row[0] == 'chrom':
                continue
            extras = []
            for region in _REGION_TYPES:
                for track in ('rm', 'sd', 'str'):
                    key = (row[0], row[1], row[2], row[3], region)
                    length = _region_length(row, region)
                    n_overlap = overlap_dict[track].get(key, 0)
                    extras.append(n_overlap / length if length else 0.0)
            writer.writerow(row + extras)

    sp.call(f'rm {feat_file}', shell=True)
    return final


# --------------------------------------------------------------------------- #
#  Public API                                                                   #
# --------------------------------------------------------------------------- #

def extract_overlap_features(del_file, dup_file, inv_file,
                              rm_bed, sd_bed, str_bed,
                              del_header, dup_header, inv_header):
    """Append repeat-element and segmental-duplication overlap fractions to feature files.

    Uses ``bedtools intersect`` to compute overlap between each SV region
    (body, flanks, confidence intervals) and three genomic annotation tracks.

    Args:
        del_file, dup_file, inv_file: paths to per-svtype feature TSV files.
        rm_bed, sd_bed, str_bed: paths to merged repeat-masker, segdup, STR BEDs.
        del_header, dup_header, inv_header: column headers for each output file.

    Returns:
        Tuple of (del_final, dup_final, inv_final) output file paths.
    """
    directory = os.path.dirname(del_file)
    uniq_id = os.path.basename(del_file).split('.')[0]

    final_files = []
    for feat_file, svtype, header in (
        (del_file, 'DEL', del_header),
        (dup_file, 'DUP', dup_header),
        (inv_file, 'INV', inv_header),
    ):
        regioned = _write_regioned_bed(feat_file)
        sorted_bed = regioned.replace('.txt', '.sorted.txt')
        sp.call(f'sort {regioned} | uniq | sortBed > {sorted_bed}', shell=True)
        sp.call(f'rm {regioned}', shell=True)

        for track, tname in ((rm_bed, 'rm'), (sd_bed, 'sd'), (str_bed, 'str')):
            ofile = sorted_bed.replace('.txt', f'.{tname}.txt')
            sp.call(
                f'intersectBed -a {sorted_bed} -b {track} -wao -sorted > {ofile}',
                shell=True,
            )

    overlap = _build_overlap_dict(directory, uniq_id)
    sp.call(f'rm {directory}/*.TMPregioned.*', shell=True)

    for feat_file, svtype, header in (
        (del_file, 'DEL', del_header),
        (dup_file, 'DUP', dup_header),
        (inv_file, 'INV', inv_header),
    ):
        final_files.append(_add_overlap_features(feat_file, svtype, overlap, header))

    return tuple(final_files)
