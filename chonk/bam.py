#!/usr/bin/env python3
"""BAM file access with header validation."""
import os
import sys
import pysam

from chonk.utils import check_path


class BamFile:
    """Wraps a pysam AlignmentFile and exposes validated sample metadata."""

    def __init__(self, path):
        check_path(path)
        self.path = os.path.abspath(path)
        self._aln = pysam.AlignmentFile(path, 'r')
        self.chrom_prefix = 'chr' in self._aln.references[0]
        self.sample = self._extract_sample()

    def _extract_sample(self):
        rg = self._aln.header.get('RG')
        if not rg:
            sys.stderr.write('FATAL ERROR: BAM is missing @RG header entry\n')
            sys.exit(1)
        samples = tuple(sorted({r['SM'] for r in rg}))
        if not samples:
            sys.stderr.write('FATAL ERROR: @RG header is missing SM field\n')
            sys.exit(1)
        if len(samples) > 1:
            sys.stderr.write(f'WARNING: multiple SM names in @RG; using {samples[0]}\n')
        return samples[0]

    @property
    def references(self):
        return self._aln.references

    @property
    def lengths(self):
        return self._aln.lengths

    def fetch(self, **kwargs):
        return self._aln.fetch(**kwargs)
