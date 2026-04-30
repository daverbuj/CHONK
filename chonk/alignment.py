#!/usr/bin/env python3
"""Data structures for parsing split-read alignments."""
import re

_CODE2CIGAR = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']
_CIGAR2CODE = {ord(c): i for i, c in enumerate(_CODE2CIGAR)}
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=XB])')


def _parse_cigar(cigar_str):
    """Parse a CIGAR string into a list of (op_code, length) tuples."""
    return [(_CIGAR2CODE[ord(op)], int(n)) for n, op in _CIGAR_RE.findall(cigar_str)]


def _query_start(cigar_tuples):
    """Return the 1-based first mapped position on the query sequence."""
    q_idx = 0
    for op, length in cigar_tuples:
        if op in (0, 1, 4, 7, 8):
            q_idx += length
        if op in (0, 7) and q_idx == length:
            return q_idx - length + 1
    return q_idx


def _ref_end(cigar_tuples, lpos):
    """Return the rightmost reference position given a left anchor."""
    span = sum(length for op, length in cigar_tuples if op in (0, 2, 3, 7, 8))
    return lpos + span - 1


class Alignment:
    """A single aligned segment, used in split-read SV analysis."""

    __slots__ = ('chrom', 'lpos', 'rpos', 'strand', 'mapq', 'qpos')

    def __init__(self):
        self.chrom = None
        self.lpos = None   # 0-based left reference position
        self.rpos = None   # right reference position
        self.strand = '+'
        self.mapq = None
        self.qpos = None   # 1-based first mapped query position

    def from_primary(self, aln):
        """Populate from a pysam AlignedSegment (primary alignment)."""
        self.chrom = aln.reference_name
        self.lpos = aln.reference_start
        self.rpos = aln.reference_end
        self.strand = '-' if aln.is_reverse else '+'
        self.mapq = aln.mapping_quality
        self.qpos = aln.query_alignment_start + 1
        return self

    def from_sa_tag(self, sa_fields):
        """Populate from a parsed SA tag field list [chrom, pos, strand, cigar, mapq, ...]."""
        self.chrom = sa_fields[0]
        self.lpos = int(sa_fields[1]) - 1
        self.strand = sa_fields[2]
        self.mapq = int(sa_fields[4])
        cigar = _parse_cigar(sa_fields[3])
        self.rpos = _ref_end(cigar, self.lpos)
        self.qpos = _query_start(cigar)
        return self
