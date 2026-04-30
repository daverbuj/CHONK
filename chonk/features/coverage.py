#!/usr/bin/env python3
"""Depth-of-coverage feature extraction."""
import csv

from chonk.utils import (
    alignment_midpoint, create_windows, get_gc_bin,
    get_gc_perc_seq, tuple2key, Welford,
)

FLANK_SIZE = 1000


def doc_fold_change(bam, meta, contig, start, end):
    """Compute depth-of-coverage fold-change for a region vs. chromosome baseline.

    Args:
        bam: BamFile object.
        meta: Metadata object for the sample.
        contig: chromosome name.
        start, end: 0-based half-open coordinates.

    Returns:
        Float fold-change (sv_doc / chrom_doc), or 0.0 if region is empty.
    """
    meta_doc = meta.doc[contig]
    region_size = end - start
    if region_size == 0 or meta_doc == 0:
        return 0.0

    rl_sum = 0
    for aln in bam.fetch(reference=contig, start=start, end=end):
        if aln.is_duplicate and aln.is_unmapped:
            continue
        try:
            mid = alignment_midpoint(aln)
        except Exception:
            continue
        if start <= mid <= end:
            rl_sum += aln.reference_length

    sv_doc = rl_sum / region_size
    return sv_doc / meta_doc


def gcbin_fold_change(bam, meta, contig, start, end, sequence):
    """Compute per-window GC-normalised read-count fold-change statistics.

    Divides the region into windows, assigns each a GC bin, counts reads,
    then normalises against the genome-wide null for that bin.

    Args:
        bam: BamFile object.
        meta: Metadata object.
        contig, start, end: region coordinates.
        sequence: DNA sequence string for the region (used for GC binning).

    Returns:
        (mean_fc, std_fc) across all windows.
    """
    region_size = end - start
    wsize = meta.window_lengths[0]
    if region_size >= 2 * meta.window_lengths[1]:
        wsize = meta.window_lengths[1]

    accumulator = Welford()
    for local_start, local_end in create_windows(wsize, 0, len(sequence)):
        win_seq = sequence[local_start:local_end]
        gc_bin = get_gc_bin(get_gc_perc_seq(win_seq))
        key = tuple2key((contig, wsize, gc_bin))

        win_start = start + local_start
        win_end = start + local_end
        count = 0
        for aln in bam.fetch(reference=contig, start=win_start, end=win_end):
            if aln.is_duplicate and aln.is_unmapped:
                continue
            try:
                mid = alignment_midpoint(aln)
            except Exception:
                continue
            if start <= mid <= end:
                count += 1

        try:
            null = float(meta.gc_rc[key])
            gc_fc = count / null if null != 0 else 0.0
        except KeyError:
            gc_fc = 0.0

        accumulator.update(gc_fc)

    return accumulator.mean, accumulator.std


def extract_coverage_features(bam, meta, chrom, start, end, sv_seq, lf_seq, rf_seq):
    """Extract nine coverage features for an SV body and its flanking regions.

    Args:
        bam: BamFile object.
        meta: Metadata object.
        chrom: chromosome name.
        start, end: SV coordinates.
        sv_seq, lf_seq, rf_seq: DNA sequences for SV body, left flank, right flank.

    Returns:
        9-tuple:
            (sv_doc_fc, sv_gc_mean, sv_gc_std,
             lf_doc_fc, lf_gc_mean, lf_gc_std,
             rf_doc_fc, rf_gc_mean, rf_gc_std)
    """
    sv_doc_fc = doc_fold_change(bam, meta, chrom, start, end)
    sv_gc_mean, sv_gc_std = gcbin_fold_change(bam, meta, chrom, start, end, sv_seq)

    lf_doc_fc = doc_fold_change(bam, meta, chrom, start - FLANK_SIZE, start)
    lf_gc_mean, lf_gc_std = gcbin_fold_change(
        bam, meta, chrom, start - FLANK_SIZE, start, lf_seq
    )

    rf_doc_fc = doc_fold_change(bam, meta, chrom, end, end + FLANK_SIZE)
    rf_gc_mean, rf_gc_std = gcbin_fold_change(
        bam, meta, chrom, end, end + FLANK_SIZE, rf_seq
    )

    return (sv_doc_fc, sv_gc_mean, sv_gc_std,
            lf_doc_fc, lf_gc_mean, lf_gc_std,
            rf_doc_fc, rf_gc_mean, rf_gc_std)


def bp_contigs(bp_file):
    """Return a tuple of distinct contigs present in a breakpoint BED file."""
    contigs = []
    with open(bp_file) as fh:
        for row in csv.reader(fh, dialect='excel-tab'):
            if row[0].startswith('#'):
                continue
            if row[0] not in contigs:
                contigs.append(row[0])
    return tuple(contigs)
