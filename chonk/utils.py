#!/usr/bin/env python3
"""Shared utility functions and statistical classes."""
import math
import os
import sys
import subprocess as sp


# --------------------------------------------------------------------------- #
#  Alignment helpers                                                            #
# --------------------------------------------------------------------------- #

def alignment_midpoint(aln):
    """Return the mapped midpoint of an alignment on the reference."""
    return int((aln.reference_end - aln.reference_start) / 2 + aln.reference_start)


def get_secondary_alns(aln):
    """Return a list of SA-tag strings, or an empty list if none."""
    return aln.get_tag('SA').split(';') if aln.has_tag('SA') else []


# --------------------------------------------------------------------------- #
#  Path / contig helpers                                                        #
# --------------------------------------------------------------------------- #

def check_path(path):
    """Exit with a fatal error message if *path* does not exist."""
    if not os.path.isfile(path):
        sys.stderr.write(f'FATAL ERROR: {path} not found\n')
        sys.exit(1)


def check_chrom(contig, chrom_prefix):
    """Normalise *contig* to match the BAM header naming convention."""
    if 'chr' in contig and not chrom_prefix:
        return contig.replace('chr', '')
    if 'chr' not in contig and chrom_prefix:
        return 'chr' + contig
    return contig


def tokenize_contigs(contig_str):
    """Split a comma-separated contig string into a tuple."""
    return tuple(contig_str.split(','))


def tuple2key(t):
    """Serialise a tuple to a comma-separated string for use as a dict key."""
    return ','.join(map(str, t))


def reporter(message):
    """Print a bordered status message to stderr."""
    pad = ' ' * 4
    bar = '-' * 80
    sys.stderr.write(f'\n{bar}\n{pad}{message}{pad}\n{bar}\n\n')


# --------------------------------------------------------------------------- #
#  GC content                                                                   #
# --------------------------------------------------------------------------- #

def get_gc_perc_fasta(fasta, contig, start, end):
    """Compute GC fraction for a genomic region using ``samtools faidx``."""
    region = f'{contig}:{start}-{end}'
    raw = sp.check_output(f'samtools faidx {fasta} {region}', shell=True)
    seq = raw.decode('utf-8').replace(f'>{region}', '').strip()
    if not seq:
        return 0.0
    return sum(1 for b in seq if b.upper() in ('G', 'C')) / len(seq)


def get_gc_perc_seq(sequence):
    """Compute GC fraction directly from a sequence string."""
    if not sequence:
        return 0.0
    return sum(1 for b in sequence if b.upper() in ('G', 'C')) / len(sequence)


def get_gc_bin(gc_content):
    """Map a GC fraction (0.0–1.0) to an integer bin index (0–25).

    Bins span 4 percentage-point intervals: bin 1 = 1–4 %, bin 2 = 5–8 %, …
    Returns -9 for non-numeric input (e.g. header lines).
    """
    try:
        pct = round(float(gc_content) * 100)
    except (ValueError, TypeError):
        return -9
    if pct % 4 == 0 and pct != 0:
        return int(pct / 4)
    return int((pct + (4 - pct % 4)) / 4)


# --------------------------------------------------------------------------- #
#  Window generation                                                            #
# --------------------------------------------------------------------------- #

def create_windows(window_size, start, end):
    """Return a list of (start, end) tuples tiling the interval by *window_size*."""
    windows = []
    pos = start
    for _ in range(int((end - start) / window_size)):
        windows.append((pos, pos + window_size - 1))
        pos += window_size
    return windows


# --------------------------------------------------------------------------- #
#  Statistics                                                                   #
# --------------------------------------------------------------------------- #

class Welford:
    """Online mean and standard deviation using Welford's algorithm.

    References:
        http://www.johndcook.com/standard_deviation.html
        https://gist.github.com/alexalemi/2151722
    """

    def __init__(self, iterable=None):
        self.k = 0
        self._mean = 0.0
        self._s = 0.0
        if iterable is not None:
            for x in iterable:
                self.update(x)

    def update(self, x):
        if x is None:
            return
        self.k += 1
        delta = x - self._mean
        self._mean += delta / self.k
        self._s += delta * (x - self._mean)

    def __call__(self, x):
        if hasattr(x, '__iter__'):
            for item in x:
                self.update(item)
        else:
            self.update(x)

    @property
    def mean(self):
        return self._mean

    @property
    def std(self):
        return math.sqrt(self._s / (self.k - 1)) if self.k > 1 else 0.0


def safe_mean(values):
    """Mean of *values*, ignoring ``None`` entries. Returns 0.0 for empty input."""
    filtered = [v for v in values if v is not None]
    return sum(filtered) / len(filtered) if filtered else 0.0


def safe_median(values):
    """Median of *values*, ignoring ``None`` entries. Returns 0.0 for empty input."""
    filtered = sorted(v for v in values if v is not None)
    n = len(filtered)
    if not n:
        return 0.0
    mid = (n - 1) // 2
    return filtered[mid] if n % 2 else (filtered[mid] + filtered[mid + 1]) / 2.0
