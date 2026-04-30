#!/usr/bin/env python3
"""K-mer junction feature extraction."""
from itertools import islice


# --------------------------------------------------------------------------- #
#  Sequence utilities                                                           #
# --------------------------------------------------------------------------- #

def _kmers(sequence, k):
    """Yield all k-mers of length *k* from *sequence*."""
    it = iter(sequence)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)


def _reverse_complement(seq):
    """Return the reverse complement of a DNA sequence string."""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(comp.get(b.upper(), 'N') for b in reversed(seq))


# --------------------------------------------------------------------------- #
#  Reference / alt k-mer construction                                          #
# --------------------------------------------------------------------------- #

def _build_ref_alt_kmers(sequence, mod_rlen, svtype, k):
    """Build sets of reference and alt junction k-mers.

    Args:
        sequence: full read-length-extended sequence spanning the SV.
        mod_rlen: modified read length used to slice flanking regions.
        svtype: 'DEL', 'DUP', or 'INV'.
        k: k-mer size.

    Returns:
        For DEL/DUP: (l_ref_kmers, r_ref_kmers, alt_kmers)
        For INV:     (l_ref_kmers, r_ref_kmers, l_alt_kmers, r_alt_kmers)
    """
    l_ref = set(_kmers(sequence[:mod_rlen * 2], k))
    r_ref = set(_kmers(sequence[-mod_rlen * 2:], k))

    if svtype == 'DEL':
        alt = set(_kmers(sequence[:mod_rlen] + sequence[-mod_rlen:], k))
        return l_ref, r_ref, alt

    if svtype == 'DUP':
        alt = set(_kmers(sequence[-mod_rlen * 2:-mod_rlen] + sequence[mod_rlen:mod_rlen * 2], k))
        return l_ref, r_ref, alt

    if svtype == 'INV':
        l_alt = set(_kmers(
            sequence[:mod_rlen] + _reverse_complement(sequence[-mod_rlen * 2:-mod_rlen]), k
        ))
        r_alt = set(_kmers(
            _reverse_complement(sequence[mod_rlen:mod_rlen * 2]) + sequence[-mod_rlen:], k
        ))
        return l_ref, r_ref, l_alt, r_alt

    return ()


# --------------------------------------------------------------------------- #
#  Read k-mer extraction                                                       #
# --------------------------------------------------------------------------- #

def _read_kmers(bam, chrom, start, end, k):
    """Collect k-mers from all reads mapped to a region."""
    kmers = set()
    for aln in bam.fetch(chrom, start, end):
        if aln.is_duplicate and aln.is_unmapped:
            continue
        seq = aln.query_sequence
        if seq is None or 'N' in seq or aln.cigarstring is None:
            continue
        kmers.update(_kmers(seq, k))
    return kmers


def _clip_kmers(bam, chrom, start, end, k, left_flank, svtype):
    """Collect k-mers from soft-clipped portions of reads.

    For DEL: left flank uses right clips; right flank uses left clips.
    For DUP: left flank uses left clips; right flank uses right clips.
    For INV: both clip orientations are included.
    """
    kmers = set()
    for aln in bam.fetch(chrom, start, end):
        if aln.cigarstring is None or 'S' not in aln.cigarstring:
            continue
        ct = aln.cigartuples

        if ct[0][0] == 4:  # left soft-clip
            if ((svtype == 'DEL' and not left_flank)
                    or (svtype == 'DUP' and left_flank)
                    or svtype == 'INV'):
                kmers.update(_kmers(aln.query_sequence[:aln.query_alignment_start], k))

        if ct[-1][0] == 4:  # right soft-clip
            if ((svtype == 'DEL' and left_flank)
                    or (svtype == 'DUP' and not left_flank)
                    or svtype == 'INV'):
                kmers.update(_kmers(aln.query_sequence[aln.query_alignment_end:], k))

    return kmers


# --------------------------------------------------------------------------- #
#  K-mer comparison                                                             #
# --------------------------------------------------------------------------- #

def _compare_kmers(l_read_kmers, r_read_kmers, ref_alt):
    """Compute overlap fractions between read k-mers and ref/alt k-mer sets.

    Returns:
        For DEL/DUP: 6-tuple (ll, lr, la, rl, rr, ra)
        For INV:     8-tuple (ll, lr, lla, lra, rl, rr, rla, rra)
    """
    lref, rref = ref_alt[0], ref_alt[1]

    def frac(query, ref):
        return len(query & ref) / len(query) if query else 0.0

    ll = frac(l_read_kmers, lref)
    lr = frac(l_read_kmers, rref)
    rl = frac(r_read_kmers, lref)
    rr = frac(r_read_kmers, rref)

    if len(ref_alt) == 3:
        alt = ref_alt[2]
        la = frac(l_read_kmers, alt)
        ra = frac(r_read_kmers, alt)
        return ll, lr, la, rl, rr, ra

    if len(ref_alt) == 4:
        lalt, ralt = ref_alt[2], ref_alt[3]
        lla = frac(l_read_kmers, lalt)
        lra = frac(l_read_kmers, ralt)
        rla = frac(r_read_kmers, lalt)
        rra = frac(r_read_kmers, ralt)
        return ll, lr, lla, lra, rl, rr, rla, rra

    return ()


def _kp_ratios(start_kmers, end_kmers, start_ref, end_ref):
    """Compute pseudo-alignment ratios between soft-clip k-mers and reference k-mers."""
    def _count(query, ref):
        return len(query & ref)

    suk_l = _count(start_kmers, start_ref)
    suk_r = _count(start_kmers, end_ref)
    euk_l = _count(end_kmers, start_ref)
    euk_r = _count(end_kmers, end_ref)

    def ratio(num, denom):
        return num / float(denom) if denom else 0.0

    return (
        ratio(suk_l, len(start_kmers)),
        ratio(suk_r, len(start_kmers)),
        ratio(euk_l, len(end_kmers)),
        ratio(euk_r, len(end_kmers)),
    )


# --------------------------------------------------------------------------- #
#  Public API                                                                   #
# --------------------------------------------------------------------------- #

def extract_kmer_features(chrom, start, end, bam, sequence, svtype,
                          mod_rlen, k, kp_k, ci_start, ci_end, lo_seq, ro_seq):
    """Extract k-mer junction (KJ) and k-mer pseudo-alignment (KP) features.

    Args:
        chrom, start, end: SV coordinates.
        bam: BamFile object.
        sequence: extended reference sequence spanning the SV plus flanks.
        svtype: 'DEL', 'DUP', or 'INV'.
        mod_rlen: modified read length (controls window size).
        k: k-mer size for KJ features.
        kp_k: k-mer size for KP features.
        ci_start, ci_end: confidence interval strings (e.g. '-20,20').
        lo_seq, ro_seq: left/right overlap sequences for KP reference k-mers.

    Returns:
        Tuple of KJ features concatenated with KP features.
        DEL/DUP: 10-tuple; INV: 12-tuple.
    """
    # --- KJ features ---
    ref_alt = _build_ref_alt_kmers(sequence, mod_rlen, svtype, k)
    l_read_kmers = _read_kmers(bam, str(chrom), start - mod_rlen, start, k)
    r_read_kmers = _read_kmers(bam, str(chrom), end, end + mod_rlen, k)
    kj = _compare_kmers(l_read_kmers, r_read_kmers, ref_alt)

    # --- KP features ---
    def _parse_ci(ci_str):
        lo, hi = ci_str.split(',')
        return int(lo), int(hi)

    cs_lo, cs_hi = _parse_ci(ci_start)
    ce_lo, ce_hi = _parse_ci(ci_end)

    start_kmers = _clip_kmers(bam, chrom, start + cs_lo, start + cs_hi, kp_k, True, svtype)
    end_kmers = _clip_kmers(bam, chrom, end + ce_lo, end + ce_hi, kp_k, False, svtype)

    start_ref = set(_kmers(lo_seq, kp_k)) if svtype != 'INV' else set(_kmers(_reverse_complement(lo_seq), kp_k))
    end_ref = set(_kmers(ro_seq, kp_k)) if svtype != 'INV' else set(_kmers(_reverse_complement(ro_seq), kp_k))

    kp = _kp_ratios(start_kmers, end_kmers, start_ref, end_ref)

    return kj + kp
