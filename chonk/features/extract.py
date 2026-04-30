#!/usr/bin/env python3
"""Feature extraction pipeline orchestrating all feature modules."""
import csv
import os
import sys
import subprocess as sp

from chonk.bam import BamFile
from chonk.metadata import Metadata
from chonk.breakpoints import select_contigs
from chonk.features.coverage import extract_coverage_features, bp_contigs
from chonk.features.fragments import extract_fragment_features
from chonk.features.kmers import extract_kmer_features
from chonk.features.context import extract_context_features
from chonk.features.overlap import extract_overlap_features
from chonk.utils import reporter


# Column headers for output files
_DELDUP_HEADER = [
    'chrom', 'start', 'end', 'svtype', 'iid', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
    'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std',
    'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std',
    'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
    'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
    'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
    'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
    'll', 'lr', 'la', 'rl', 'rr', 'ra',
    'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
    'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc',
    'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp',
    'log_sv_len', 'bp_start_ci', 'bp_end_ci',
]

_INV_HEADER = [
    'chrom', 'start', 'end', 'svtype', 'id', 'gt', 'ci_start', 'ci_end', 'k', 'rlen_mod',
    'sv_doc_fc', 'sv_gc_mean', 'sv_gc_std',
    'lf_doc_fc', 'lf_gc_mean', 'lf_gc_std',
    'rf_doc_fc', 'rf_gc_mean', 'rf_gc_std',
    'sf_ratio', 'split_ratio', 'disc_ratio', 'clip_ratio',
    'sf_mapq_mean', 'sf_mapq_median', 'nonsf_mapq_mean', 'nonsf_mapq_median',
    'sf_baseq_mean', 'sf_baseq_median', 'nonsf_baseq_mean', 'nonsf_baseq_median',
    'll', 'lr', 'lla', 'lra', 'rl', 'rr', 'rla', 'rra',
    'start_ratio_left', 'start_ratio_right', 'end_ratio_left', 'end_ratio_right',
    'sv_gc', 'lf_gc', 'rf_gc', 'lo_gc', 'ro_gc',
    'sv_comp', 'lf_comp', 'rf_comp', 'lo_comp', 'ro_comp',
    'log_sv_len', 'bp_start_ci', 'bp_end_ci',
]


# --------------------------------------------------------------------------- #
#  Internal helpers                                                             #
# --------------------------------------------------------------------------- #

def _load_metadata_dict(meta_dir):
    """Build a sample-name → Metadata mapping from all JSON files in *meta_dir*."""
    meta_dict = {}
    for fname in os.listdir(meta_dir):
        if not fname.endswith('chonk.json'):
            continue
        meta = Metadata()
        meta.load(os.path.join(meta_dir, fname))
        sample = fname[:fname.find('_')]
        meta_dict[sample] = meta
    return meta_dict


def _normalise_ci(ci_str):
    """Ensure CI string has the form 'negative,positive'."""
    if ci_str in ('.', '0,0'):
        return '-20,20'
    lo, hi = int(ci_str.split(',')[0]), int(ci_str.split(',')[1])
    return f'{min(lo, -abs(lo))},{max(hi, abs(hi))}'


def _encode_genotype(gt_str):
    """Map a VCF genotype string to an integer label (0, 1, 2)."""
    if '1' not in gt_str:
        return 0
    if '0' not in gt_str:
        return 2
    return 1


def _build_flanked_bed(bp_file, meta_dict, rlen_modifier, fasta, outdir, uniq_id):
    """Create a temporary BED file with flanked SV coordinates and attached FASTA sequence.

    Flanks of 1000 bp are added to each side so that sequence context is available
    for k-mer and GC-content features.

    Returns:
        Path to the temporary FASTA-annotated BED file.
    """
    flanked_bed = os.path.join(outdir, f'{uniq_id}.flanked_bpfile.TEMP.bed')
    with open(bp_file) as fh, open(flanked_bed, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        for row in csv.reader(fh, dialect='excel-tab'):
            if row[0].startswith('#'):
                continue
            chrom, start, end, svtype, ci_start, ci_end, sample, genotype = row[:8]
            start, end = int(start), int(end)
            mod_rlen = int(float(meta_dict[sample].read_len[chrom]) * float(rlen_modifier))
            writer.writerow([chrom, start - 1000, end + 1000,
                             svtype, ci_start, ci_end, sample, genotype, mod_rlen])

    fasta_bed = os.path.join(outdir, f'{uniq_id}.fasta_bpfile.TEMP.bed')
    sp.call(f'bedtools getfasta -bedOut -fi {fasta} -bed {flanked_bed} > {fasta_bed}', shell=True)
    sp.call(f'rm {flanked_bed}', shell=True)
    return fasta_bed


def _parse_sequences(sequence, start, end, seq_zero, ci_start, ci_end, mod_rlen):
    """Slice a long sequence string into SV-region and flanking sub-sequences.

    Args:
        sequence: full sequence string (spanning start-1000 to end+1000).
        start, end: SV coordinates (absolute).
        seq_zero: start coordinate of *sequence* (used to convert to local indices).
        ci_start, ci_end: parsed confidence-interval strings.
        mod_rlen: modified read length for kmer windows.

    Returns:
        (sv_seq, rlen_seq, lo_seq, l1k, ro_seq, r1k, lf_seq, rf_seq)
    """
    sv_s = start - seq_zero
    sv_e = end - seq_zero + 1
    sv_seq = sequence[sv_s:sv_e]

    lf_s = sv_s - mod_rlen
    rf_e = sv_e + mod_rlen
    rlen_seq = sequence[lf_s:rf_e]
    lf_seq = sequence[lf_s:sv_s]
    rf_seq = sequence[sv_e:rf_e]

    cs_lo = int(float(ci_start.split(',')[0]))
    cs_hi = int(float(ci_start.split(',')[1]))
    ce_lo = int(float(ci_end.split(',')[0]))
    ce_hi = int(float(ci_end.split(',')[1]))

    lo_seq = sequence[sv_s + cs_lo: sv_s + cs_hi]
    ro_seq = sequence[sv_e + ce_lo: sv_e + ce_hi]

    l1k = sequence[:1000]
    r1k = sequence[-1000:]

    return sv_seq, rlen_seq, lo_seq, l1k, ro_seq, r1k, lf_seq, rf_seq


# --------------------------------------------------------------------------- #
#  Pipeline entry point                                                         #
# --------------------------------------------------------------------------- #

def run_features(args):
    """Extract all features for every SV in the breakpoint file.

    Writes three output files (one per SV type: DEL, DUP, INV).
    If repeat-track BED files are provided, appends overlap features.

    Args:
        args: namespace with attributes:
            bp          - path to breakpoint BED file
            metadir     - directory containing per-sample metadata JSON files
            fasta       - path to reference FASTA
            k           - k-mer size for KJ features
            rlen        - read-length modifier
            pk          - k-mer size for KP features
            o           - output file path (base name, must end in .txt)
            rm          - (optional) repeat-masker BED
            sd          - (optional) segdup BED
            str_bed     - (optional) STR BED
    """
    if not args.o.endswith('.txt'):
        sys.stderr.write("FATAL ERROR: -o must end in '.txt'\n")
        sys.exit(1)

    outdir = os.path.dirname(os.path.abspath(args.o))
    fname = os.path.basename(args.o)
    uniq_id = fname[:fname.find('.')]

    reporter('Loading metadata')
    meta_dict = _load_metadata_dict(args.metadir)

    reporter('Building flanked BED with FASTA sequences')
    fasta_bed = _build_flanked_bed(
        args.bp, meta_dict, args.rlen, args.fasta, outdir, uniq_id
    )

    # Open output files
    base = args.o[:args.o.rfind('.')]
    del_path = base + '.DEL.txt'
    dup_path = base + '.DUP.txt'
    inv_path = base + '.INV.txt'

    del_fh = open(del_path, 'w', newline='')
    dup_fh = open(dup_path, 'w', newline='')
    inv_fh = open(inv_path, 'w', newline='')
    del_writer = csv.writer(del_fh, delimiter='\t')
    dup_writer = csv.writer(dup_fh, delimiter='\t')
    inv_writer = csv.writer(inv_fh, delimiter='\t')
    del_writer.writerow(_DELDUP_HEADER)
    dup_writer.writerow(_DELDUP_HEADER)
    inv_writer.writerow(_INV_HEADER)

    csv.field_size_limit(int(sys.maxsize / 100))

    reporter('Extracting features')
    count = 0
    with open(fasta_bed) as fh:
        for row in csv.reader(fh, dialect='excel-tab'):
            try:
                chrom, f_start, f_end, svtype, cistart, ciend, sample, genotype, mod_rlen, sequence = row
            except ValueError:
                sys.stderr.write(f'WARNING: skipping malformed row: {row}\n')
                continue

            f_start, f_end = int(f_start), int(f_end)
            if f_end - f_start > 5e6:
                continue

            start, end = f_start + 1000, f_end - 1000
            ci_start = _normalise_ci(cistart)
            ci_end = _normalise_ci(ciend)
            mod_rlen = int(mod_rlen)
            rlen = int(mod_rlen / float(args.rlen))
            gt = _encode_genotype(genotype)

            sv_seq, rlen_seq, lo_seq, l1k, ro_seq, r1k, lf_seq, rf_seq = _parse_sequences(
                sequence, start, end, f_start, ci_start, ci_end, mod_rlen
            )

            meta = meta_dict[sample]
            bam = BamFile(meta.bam_path)

            cov_feats = extract_coverage_features(
                bam, meta, chrom, start, end, sv_seq, lf_seq, rf_seq
            )
            sf_feats = extract_fragment_features(
                bam, meta, chrom, start, end, ci_start, ci_end, svtype
            )
            km_feats = extract_kmer_features(
                chrom, start, end, bam, rlen_seq, svtype,
                mod_rlen, args.k, args.pk, ci_start, ci_end, lo_seq, ro_seq
            )
            ctx_feats = extract_context_features(
                start, end, ci_start, ci_end, rlen,
                sv_seq, lf_seq, rf_seq, lo_seq, ro_seq
            )

            base_row = (chrom, start, end, svtype, sample, gt,
                        ci_start, ci_end, args.k, args.rlen)
            out_row = base_row + cov_feats + sf_feats + km_feats + ctx_feats

            if svtype == 'DEL':
                del_writer.writerow(out_row)
            elif svtype == 'DUP':
                dup_writer.writerow(out_row)
            elif svtype == 'INV':
                inv_writer.writerow(out_row)

            count += 1

    del_fh.close()
    dup_fh.close()
    inv_fh.close()
    sp.call(f'rm {fasta_bed}', shell=True)

    reporter(f'Feature extraction complete — {count} SVs processed')

    # Optional overlap features
    if getattr(args, 'rm', None) and getattr(args, 'sd', None) and getattr(args, 'str_bed', None):
        reporter('Computing overlap features')
        extract_overlap_features(
            del_path, dup_path, inv_path,
            args.rm, args.sd, args.str_bed,
            _DELDUP_HEADER, _DELDUP_HEADER, _INV_HEADER,
        )
        reporter('Overlap features complete')
