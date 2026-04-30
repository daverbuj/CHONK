#!/usr/bin/env python3
"""Supporting fragment (split-read, discordant PE, clip) feature extraction."""
from chonk.alignment import Alignment
from chonk.utils import get_secondary_alns, safe_mean, safe_median


# --------------------------------------------------------------------------- #
#  Fragment class                                                               #
# --------------------------------------------------------------------------- #

class Fragment:
    """Tracks SV-supporting evidence and quality scores for one sequencing fragment.

    Evidence types are prioritised: SR > DPE > CLIP > non-supporting.
    """

    _PRIORITY = {'SR': 3, 'DPE': 2, 'CLIP': 1, None: 0}

    def __init__(self):
        self.supp = False
        self.split = False
        self.disc = False
        self.clip = False
        self.mapq = []
        self.baseq = []
        self.n_forward = 0
        self.n_reverse = 0
        self.n_aln = 0
        self.qname = None

    def update(self, evidence):
        """Update this fragment with new SV evidence.

        Args:
            evidence: (caller, (l_strand, l_mapq), (r_strand, r_mapq), basequals, qname)
                      caller is one of 'SR', 'DPE', 'CLIP', or None.
        """
        caller, left, right, basequals, qname = evidence
        l_strand, l_mapq = left
        r_strand, r_mapq = right

        current = self._PRIORITY[
            'SR' if self.split else ('DPE' if self.disc else ('CLIP' if self.clip else None))
        ]
        incoming = self._PRIORITY[caller]

        if caller is None:
            # Non-supporting read: accumulate base qualities and strand counts
            self.baseq.extend(basequals or [])
            self.qname = qname
            self._count_strand(l_strand)
            self.n_aln += 1
            return

        if incoming < current:
            return

        # Append additional SR evidence when the fragment already has a split read
        if caller == 'SR' and self.split:
            self.mapq.extend([l_mapq, r_mapq])
            self.baseq.extend(basequals or [])
            self.n_aln += 2
            self._count_strand(l_strand)
            self._count_strand(r_strand)
            return

        # Write / overwrite with the new (higher-priority) evidence
        self.supp = True
        self.split = caller == 'SR'
        self.disc = caller == 'DPE'
        self.clip = caller == 'CLIP'
        self.n_forward = self.n_reverse = 0
        self.baseq = list(basequals or [])
        self.qname = qname

        if caller == 'CLIP':
            self.mapq = [l_mapq]
            self.n_aln = 1
            self._count_strand(l_strand)
        else:
            self.mapq = [l_mapq, r_mapq]
            self.n_aln = 2
            self._count_strand(l_strand)
            self._count_strand(r_strand)

    def update_baseq(self, qualities):
        """Append additional base-quality scores to this fragment."""
        self.baseq.extend(qualities or [])

    def _count_strand(self, strand):
        if strand == '+':
            self.n_forward += 1
        elif strand == '-':
            self.n_reverse += 1


# --------------------------------------------------------------------------- #
#  Read classification helpers                                                  #
# --------------------------------------------------------------------------- #

def _classify_split_read(aln, secondary_alns, chrom, owin_start, owin_end, true_svtype):
    """Return SR evidence tuple if *aln* supports *true_svtype*, else empty tuple."""
    left = Alignment().from_primary(aln)

    for sa in secondary_alns:
        if not sa:
            continue
        fields = sa.split(',')
        if fields[0] != chrom:
            continue
        right = Alignment().from_sa_tag(fields)

        in_window = right.rpos > owin_start and right.lpos < owin_end
        overlap = set(range(left.lpos, left.rpos + 1)) & set(range(right.lpos, right.rpos + 1))
        if not in_window or overlap:
            continue

        if left.strand == right.strand:
            svtype = 'DUP' if left.qpos > right.qpos else 'DEL'
        else:
            svtype = 'INV'

        if svtype == true_svtype:
            return (
                'SR',
                (left.strand, aln.mapping_quality),
                (right.strand, aln.mapping_quality),
                aln.query_alignment_qualities,
                aln.query_name,
            )
    return ()


def _classify_discordant_pe(aln, meta, chrom, owin_start, owin_end, true_svtype):
    """Return DPE evidence tuple if *aln* supports *true_svtype*, else empty tuple."""
    tlen_mean = float(meta.tlen[chrom])
    tlen_std = float(meta.tlen_std[chrom])

    if aln.mate_is_unmapped:
        return ()
    try:
        _ = aln.reference_end < aln.next_reference_start
    except Exception:
        return ()

    if not (aln.is_paired
            and aln.next_reference_name == aln.reference_name
            and aln.next_reference_start < owin_end
            and (aln.next_reference_start + aln.query_length) > owin_start
            and aln.reference_end < aln.next_reference_start):
        return ()

    if abs(aln.template_length) <= tlen_mean + 3.5 * tlen_std:
        return ()

    if aln.is_reverse != aln.mate_is_reverse:
        svtype = 'DUP' if aln.is_reverse else 'DEL'
        l_strand = '-' if aln.is_reverse else '+'
        r_strand = '+' if aln.is_reverse else '-'
    else:
        svtype = 'INV'
        l_strand = r_strand = '-' if aln.is_reverse else '+'

    if svtype != true_svtype:
        return ()

    l_mapq = aln.mapping_quality
    try:
        r_mapq = aln.get_tag('MQ')
    except Exception:
        r_mapq = 0

    return ('DPE', (l_strand, l_mapq), (r_strand, r_mapq),
            aln.query_alignment_qualities, aln.query_name)


def _is_clipped(aln):
    """Return (left_clip, right_clip) booleans for soft/hard clips > 20 bp."""
    cigar = aln.cigartuples
    if cigar is None:
        return False, False
    try:
        left_clip = cigar[0][0] in (4, 5) and cigar[0][1] > 20
        right_clip = cigar[-1][0] in (4, 5) and cigar[-1][1] > 20
    except (IndexError, TypeError):
        return False, False
    return left_clip, right_clip


def _classify_clipped(aln, true_svtype, left_window=True):
    """Return CLIP evidence tuple if *aln* is informatively clipped, else empty tuple."""
    if aln.has_tag('SA'):
        return ()
    left_clip, right_clip = _is_clipped(aln)
    if not left_clip and not right_clip:
        return ()

    if true_svtype == 'DEL':
        informative = (left_window and right_clip) or (not left_window and left_clip)
    elif true_svtype == 'DUP':
        informative = (left_window and left_clip) or (not left_window and right_clip)
    elif true_svtype == 'INV':
        informative = left_clip or right_clip
    else:
        informative = False

    if not informative:
        return ()

    strand = '-' if aln.is_reverse else '+'
    return ('CLIP', (strand, aln.mapping_quality), (None, None),
            aln.query_alignment_qualities, aln.query_name)


# --------------------------------------------------------------------------- #
#  Window scanning                                                              #
# --------------------------------------------------------------------------- #

def _scan_left_window(bam, meta, contig, lwin_s, lwin_e, rwin_s, rwin_e, svtype, fragments):
    """Scan the left breakpoint window, classifying each read."""
    for aln in bam.fetch(reference=contig, start=int(lwin_s), end=int(lwin_e)):
        if aln.is_duplicate and aln.is_unmapped:
            continue
        if aln.query_name not in fragments:
            fragments[aln.query_name] = Fragment()

        secondary = get_secondary_alns(aln)
        if secondary:
            ev = _classify_split_read(aln, secondary, contig, rwin_s, rwin_e, svtype)
            if ev:
                fragments[aln.query_name].update(ev)
                continue

        ev = _classify_discordant_pe(aln, meta, contig, rwin_s, rwin_e, svtype)
        if ev:
            fragments[aln.query_name].update(ev)
            continue

        ev = _classify_clipped(aln, svtype, left_window=True)
        if ev:
            fragments[aln.query_name].update(ev)
            continue

        strand = '-' if aln.is_reverse else '+'
        fragments[aln.query_name].update(
            (None, (strand, None), (None, None), aln.query_alignment_qualities, aln.query_name)
        )


def _scan_right_window(bam, contig, rwin_s, rwin_e, svtype, fragments):
    """Scan the right breakpoint window, updating base qualities and adding new reads."""
    for aln in bam.fetch(reference=contig, start=int(rwin_s), end=int(rwin_e)):
        if aln.is_duplicate and aln.is_unmapped:
            continue

        if aln.query_name in fragments:
            # Already processed in left window: just accumulate base qualities
            fragments[aln.query_name].update_baseq(aln.query_alignment_qualities)
        else:
            fragments[aln.query_name] = Fragment()
            ev = _classify_clipped(aln, svtype, left_window=True)
            if ev:
                fragments[aln.query_name].update(ev)
            else:
                strand = '-' if aln.is_reverse else '+'
                fragments[aln.query_name].update(
                    (None, (strand, None), (None, None),
                     aln.query_alignment_qualities, aln.query_name)
                )


# --------------------------------------------------------------------------- #
#  Feature aggregation                                                          #
# --------------------------------------------------------------------------- #

def _aggregate_fragments(fragments):
    """Summarise a dict of Fragment objects into scalar features.

    Returns:
        12-tuple:
            (sf_ratio, split_ratio, disc_ratio, clip_ratio,
             sf_mapq_mean, sf_mapq_median, nonsf_mapq_mean, nonsf_mapq_median,
             sf_baseq_mean, sf_baseq_median, nonsf_baseq_mean, nonsf_baseq_median)
    """
    n_total = len(fragments)
    n_sf = n_split = n_disc = n_clip = 0
    sf_mapq = []
    nonsf_mapq = []
    sf_baseq = []
    nonsf_baseq = []

    for frag in fragments.values():
        if frag.supp:
            n_sf += 1
            sf_mapq.extend(frag.mapq)
            sf_baseq.extend(frag.baseq)
        else:
            nonsf_mapq.extend(frag.mapq)
            nonsf_baseq.extend(frag.baseq)
        if frag.split:
            n_split += 1
        if frag.disc:
            n_disc += 1
        if frag.clip:
            n_clip += 1

    def ratio(n):
        return n / n_total if n_total else 0.0

    return (
        ratio(n_sf), ratio(n_split), ratio(n_disc), ratio(n_clip),
        safe_mean(sf_mapq), safe_median(sf_mapq),
        safe_mean(nonsf_mapq), safe_median(nonsf_mapq),
        safe_mean(sf_baseq), safe_median(sf_baseq),
        safe_mean(nonsf_baseq), safe_median(nonsf_baseq),
    )


def _build_ci_windows(start, end, ci_start, ci_end, tlen_mean, read_len):
    """Compute left/right confidence-interval windows around SV breakpoints."""
    def parse(ci_str, pos):
        lo, hi = ci_str.split(',')
        return int(pos) + int(float(lo)), int(pos) + int(float(hi))

    lwin_s, lwin_e = parse(ci_start, start)
    rwin_s, rwin_e = parse(ci_end, end)

    if lwin_e > rwin_s:
        half = int(read_len / 2)
        lwin_s, lwin_e = int(start) - half, int(start) + half
        rwin_s, rwin_e = int(end) - half, int(end) + half

    if lwin_e > rwin_s:
        lwin_s, lwin_e = int(start) - 20, int(start) + 20
        rwin_s, rwin_e = int(end) - 20, int(end) + 20

    return lwin_s, lwin_e, rwin_s, rwin_e


# --------------------------------------------------------------------------- #
#  Public API                                                                   #
# --------------------------------------------------------------------------- #

def extract_fragment_features(bam, meta, contig, start, end, ci_start, ci_end, svtype):
    """Extract supporting-fragment features for a structural variant.

    Args:
        bam: BamFile object.
        meta: Metadata object.
        contig, start, end: SV coordinates.
        ci_start, ci_end: confidence-interval strings (e.g. '-20,20').
        svtype: 'DEL', 'DUP', or 'INV'.

    Returns:
        12-tuple of supporting-fragment features.
    """
    lwin_s, lwin_e, rwin_s, rwin_e = _build_ci_windows(
        start, end, ci_start, ci_end,
        meta.tlen[contig], meta.read_len[contig],
    )

    fragments = {}
    _scan_left_window(bam, meta, contig, lwin_s, lwin_e, rwin_s, rwin_e, svtype, fragments)
    _scan_right_window(bam, contig, rwin_s, rwin_e, svtype, fragments)

    return _aggregate_fragments(fragments)
