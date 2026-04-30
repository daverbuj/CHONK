#!/usr/bin/env python3
"""SV breakpoint detection via split-reads and discordant paired-end reads."""
import csv
import operator
import os
import sys

from chonk.alignment import Alignment
from chonk.bam import BamFile
from chonk.metadata import Metadata
from chonk.utils import check_chrom, get_secondary_alns, reporter


# --------------------------------------------------------------------------- #
#  Split-read calling                                                           #
# --------------------------------------------------------------------------- #

def call_split_read(aln, secondary_alns, chrom):
    """Infer SV breakpoints from a split-read alignment.

    Args:
        aln: pysam AlignedSegment (primary alignment).
        secondary_alns: list of SA-tag strings from the same read.
        chrom: contig being processed.

    Returns:
        List of (chrom, start, end, svtype) tuples.
    """
    alns = [Alignment().from_primary(aln)]
    for sa in secondary_alns:
        if not sa:
            continue
        fields = sa.split(',')
        if fields[0] != chrom:
            continue
        alns.append(Alignment().from_sa_tag(fields))

    if len(alns) < 2:
        return []

    alns.sort(key=operator.attrgetter('qpos'))
    breaks = []

    for i in range(len(alns) - 1):
        left, right = alns[i], alns[i + 1]
        if left.chrom != right.chrom:
            continue
        # Skip if reference positions overlap
        if set(range(left.lpos, left.rpos + 1)) & set(range(right.lpos, right.rpos + 1)):
            continue

        if left.strand == right.strand:
            if left.lpos <= right.rpos:
                start, end, svtype = left.rpos, right.lpos, 'DEL'
            else:
                start, end, svtype = right.lpos, left.rpos, 'DUP'
        else:
            svtype = 'INV'
            start, end = left.rpos, right.rpos
            if aln.mate_is_reverse:
                start, end = left.lpos, right.lpos
            if start > end:
                start, end = end, start

        breaks.append((chrom, start, end, svtype))

    return breaks


# --------------------------------------------------------------------------- #
#  Discordant paired-end calling                                               #
# --------------------------------------------------------------------------- #

def call_discordant_pe(aln, meta, chrom):
    """Detect an SV from a discordant paired-end read.

    Args:
        aln: pysam AlignedSegment.
        meta: Metadata object for the sample.
        chrom: contig being processed.

    Returns:
        (chrom, start, end, svtype) tuple, or ``None``.
    """
    tlen_mean = float(meta.tlen[chrom])
    tlen_std = float(meta.tlen_std[chrom])
    threshold = tlen_mean + 3.5 * tlen_std

    # Opposite-strand pairs → DEL or DUP
    if (aln.is_paired
            and aln.next_reference_name == aln.reference_name
            and aln.is_reverse != aln.mate_is_reverse
            and not aln.is_reverse):
        if abs(aln.template_length) > threshold:
            if aln.reference_start < aln.next_reference_start:
                return (chrom, aln.reference_end, aln.next_reference_start, 'DEL')
            return (chrom, aln.next_reference_start, aln.reference_end, 'DUP')

    # Same-strand pairs → INV
    if (aln.is_paired
            and aln.next_reference_name == aln.reference_name
            and aln.is_reverse == aln.mate_is_reverse
            and aln.is_read1):
        if abs(aln.template_length) > threshold:
            if not aln.is_reverse and not aln.mate_is_reverse:
                start, end = sorted([
                    aln.reference_end,
                    aln.next_reference_start + aln.query_length,
                ])
                return (chrom, start, end, 'INV')
            if aln.is_reverse and not aln.mate_is_reverse:
                start, end = sorted([aln.reference_start, aln.next_reference_start])
                return (chrom, start, end, 'INV')

    return None


# --------------------------------------------------------------------------- #
#  Pipeline entry point                                                         #
# --------------------------------------------------------------------------- #

def run_breakpoints(args):
    """Detect SV breakpoints and write results to a BED file."""
    meta = Metadata()
    meta.load(args.m)

    bam = BamFile(meta.bam_path)
    contigs = select_contigs(args.r, meta.user_contigs, bam)

    if contigs == meta.contigs:
        out_path = os.path.join(args.o, f'{meta.sample}_breakpoints.chonk.bed')
    else:
        out_path = os.path.join(args.o, f'{meta.sample}_{"_".join(contigs)}_breakpoints.chonk.bed')

    processed = []
    with open(out_path, 'wt') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['#chrom', 'start', 'end', 'svtype'])

        for chrom in contigs:
            for aln in bam.fetch(reference=chrom):
                if aln.mapq == 0:
                    continue

                secondary = get_secondary_alns(aln)
                breaks = []

                if secondary:
                    read_id = aln.query_name + ('2' if aln.is_read2 else '1')
                    if read_id in processed:
                        processed.remove(read_id)
                        continue
                    breaks = call_split_read(aln, secondary, chrom)
                    for sv in breaks:
                        writer.writerow(sv)
                    processed.append(read_id)

                if not breaks:
                    sv = call_discordant_pe(aln, meta, chrom)
                    if sv:
                        writer.writerow(sv)

    reporter(f'breakpoints complete  →  {os.path.abspath(out_path)}')


def select_contigs(user_contigs, available, bam):
    """Return the tuple of contigs to process, validated against the metadata.

    Args:
        user_contigs: tuple of user-specified contigs, or ``None`` for all.
        available: tuple of contigs present in the metadata file.
        bam: BamFile used to normalise chromosome naming.

    Returns:
        Tuple of validated contig names.
    """
    if user_contigs is None:
        return available
    selected = []
    for c in user_contigs:
        c = check_chrom(c, bam.chrom_prefix)
        if c not in available:
            sys.stderr.write(
                f'FATAL ERROR: {c} not found in metadata; re-run chonk metadata\n'
            )
            sys.exit(1)
        selected.append(c)
    if not selected:
        sys.stderr.write('FATAL ERROR: no valid contigs specified\n')
        sys.exit(1)
    return tuple(selected)
