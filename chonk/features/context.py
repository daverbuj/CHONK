#!/usr/bin/env python3
"""Sequence context feature extraction (GC content, complexity, SV length, CI)."""
import sys
import zlib
from math import log10

import numpy as np


def _gc_fraction(sequence):
    """Return GC fraction for a sequence string, or NaN if empty."""
    if not sequence:
        return np.nan
    gc = sum(1 for b in sequence if b.upper() in ('G', 'C'))
    return gc / len(sequence)


def _complexity(sequence):
    """Return sequence complexity as compressed_size / uncompressed_size."""
    seq = sequence.upper().encode('utf-8')
    empty = ''.encode('utf-8')
    comp_sz = sys.getsizeof(zlib.compress(seq)) - sys.getsizeof(zlib.compress(empty))
    uncomp_sz = sys.getsizeof(seq) - sys.getsizeof(empty)
    return comp_sz / uncomp_sz if uncomp_sz > 0 else float('nan')


def _parse_ci(ci_str):
    """Parse a confidence-interval string (e.g. '-20,20') into (lo, hi) ints."""
    lo, hi = ci_str.split(',')
    return int(float(lo)), int(float(hi))


def extract_context_features(start, end, ci_start, ci_end, rlen,
                              sv_seq, lf_seq, rf_seq, lo_seq, ro_seq):
    """Extract sequence-context features for an SV.

    Args:
        start, end: SV coordinates.
        ci_start, ci_end: confidence-interval strings.
        rlen: read length (for normalising CI width).
        sv_seq: SV body sequence.
        lf_seq, rf_seq: left/right flank sequences.
        lo_seq, ro_seq: left/right overlap sequences.

    Returns:
        13-tuple:
            (sv_gc, lf_gc, rf_gc, lo_gc, ro_gc,
             sv_comp, lf_comp, rf_comp, lo_comp, ro_comp,
             log_sv_len, bp_start_ci, bp_end_ci)
    """
    # GC content
    sv_gc = _gc_fraction(sv_seq)
    lf_gc = _gc_fraction(lf_seq)
    rf_gc = _gc_fraction(rf_seq)
    lo_gc = _gc_fraction(lo_seq)
    ro_gc = _gc_fraction(ro_seq)

    # Sequence complexity
    sv_comp = _complexity(sv_seq)
    lf_comp = _complexity(lf_seq)
    rf_comp = _complexity(rf_seq)
    lo_comp = _complexity(lo_seq)
    ro_comp = _complexity(ro_seq)

    # SV length (log10)
    log_sv_len = log10(end - start) if end > start else float('nan')

    # Breakpoint confidence interval widths, normalised by read length
    cs_lo, cs_hi = _parse_ci(ci_start)
    ce_lo, ce_hi = _parse_ci(ci_end)
    bp_start_ci = (cs_hi - cs_lo) / rlen if rlen else float('nan')
    bp_end_ci = (ce_hi - ce_lo) / rlen if rlen else float('nan')

    return (sv_gc, lf_gc, rf_gc, lo_gc, ro_gc,
            sv_comp, lf_comp, rf_comp, lo_comp, ro_comp,
            log_sv_len, bp_start_ci, bp_end_ci)
