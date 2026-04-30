#!/usr/bin/env python3
"""Per-sample genome metadata extraction for SV calling."""
import csv
import json
import os
import random
import sys
from operator import itemgetter
import subprocess as sp

from chonk.bam import BamFile
from chonk.utils import (
    alignment_midpoint, check_chrom, get_gc_bin,
    reporter, tuple2key, Welford,
)


class Metadata:
    """Per-sample sequencing metadata required by downstream SV detection.

    Stores depth-of-coverage, read-length, template-length, and GC-normalised
    read-count statistics computed once per sample and persisted as JSON.
    """

    def __init__(self):
        # Tuneable parameters (serialised to JSON so downstream steps stay in sync)
        self.window_lengths = (25, 1000)
        self.n_windows = int(1e4)
        self.n_reads = int(5e6)
        self.tlen_cap = int(1e4)

        # File / sample identity
        self.bam_path = None
        self.sample = None
        self.outdir = None
        self.tmpdir = None
        self.fasta = None
        self.json_file = None

        # Contig lists
        self.contigs = None       # all contigs from BAM header
        self.user_contigs = None  # contigs actually processed

        # Per-contig sequencing statistics
        self.doc = {}       # depth of coverage
        self.read_len = {}  # mean read length
        self.tlen = {}      # mean insert (template) length
        self.tlen_std = {}  # insert-length standard deviation

        # GC-binned read-count null model (key = tuple2key(contig, window_len, gc_bin))
        self.gc_rc = {}
        self.gc_std = {}

    # ----------------------------------------------------------------------- I/O

    def save(self):
        """Serialise the metadata object to its JSON file."""
        with open(self.json_file, 'w') as fh:
            json.dump(self.__dict__, fh, indent=4)

    def load(self, json_file):
        """Populate this object from a previously saved JSON file."""
        with open(json_file) as fh:
            data = json.load(fh)
        for key, value in data.items():
            setattr(self, key, value)

    # ----------------------------------------------------------------- Setup

    def initialize(self, bam: BamFile, args):
        """Configure paths and contig lists from a BamFile and parsed CLI args."""
        self.bam_path = bam.path
        self.sample = bam.sample
        self.outdir = args.o
        self.fasta = os.path.abspath(args.f)
        self.contigs = tuple(bam.references)
        self.user_contigs = self._resolve_contigs(args.r, bam.chrom_prefix)

        # Temporary working directory
        tag = '_'.join(self.user_contigs) if self.user_contigs != self.contigs else ''
        tmpdir_name = f'{self.sample}_{tag}_tmp' if tag else f'{self.sample}_tmp'
        tmpdir = os.path.join(args.o, tmpdir_name)
        os.makedirs(tmpdir, exist_ok=True)
        self.tmpdir = tmpdir + os.sep

        # Output JSON path
        name = (f'{self.sample}_{tag}_metadata.chonk.json' if tag
                else f'{self.sample}_metadata.chonk.json')
        self.json_file = os.path.join(self.outdir, name)

    def _resolve_contigs(self, user_contigs, chrom_prefix):
        if user_contigs is None:
            return self.contigs
        resolved = []
        for c in user_contigs:
            c = check_chrom(c, chrom_prefix)
            if c not in self.contigs:
                sys.stderr.write(f'WARNING: {c} not found in BAM header\n')
            else:
                resolved.append(c)
        if not resolved:
            sys.stderr.write('FATAL ERROR: no valid contigs specified\n')
            sys.exit(1)
        return tuple(resolved)


# --------------------------------------------------------------------------- #
#  Pipeline entry point                                                         #
# --------------------------------------------------------------------------- #

def run_metadata(args):
    """Extract genome metadata from a BAM file and save to JSON."""
    random.seed(args.s)

    bam = BamFile(args.i)
    meta = Metadata()
    meta.initialize(bam, args)

    loci = {c: [] for c in meta.user_contigs}
    binned_loci = {c: [] for c in meta.user_contigs}

    mask = _build_mask(bam, args.x, meta, loci)
    gc_beds = _bin_gc_content(_window_mask(mask, args.f, meta.tmpdir, meta.window_lengths))
    _load_loci(meta, binned_loci, gc_beds)
    _compute_contig_metadata(meta, bam, loci)
    _compute_gc_metadata(meta, bam, binned_loci)

    meta.save()
    reporter(f'metadata complete  →  {os.path.abspath(meta.json_file)}')


# --------------------------------------------------------------------------- #
#  Internal pipeline helpers                                                    #
# --------------------------------------------------------------------------- #

def _build_mask(bam, exclude, meta, loci):
    """Create a BED mask of accessible regions and populate *loci*."""
    contig_bed = os.path.join(meta.tmpdir, 'contigs.bed')
    with open(contig_bed, 'wt') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        for name, length in zip(bam.references, bam.lengths):
            if name in meta.user_contigs:
                writer.writerow([name, 0, length])

    if exclude is None:
        mask = contig_bed
    else:
        raw = os.path.join(meta.tmpdir, 'masked.bed')
        sp.call(f'bedtools subtract -a {contig_bed} -b {os.path.abspath(exclude)} > {raw}', shell=True)
        mask = os.path.join(meta.tmpdir, 'masked.gt10kb.bed')
        sp.call(
            f"""awk 'BEGIN {{OFS="\t"}}; $3-$2>=10000 {{print $1,$2,$3}}' {raw} > {mask}""",
            shell=True,
        )

    with open(mask) as fh:
        for line in fh:
            parts = line.rstrip().split('\t')
            loci[parts[0]].append((int(parts[1]), int(parts[2])))

    return mask


def _window_mask(mask, fasta, tmpdir, window_lengths):
    """Tile the mask into windows of each size and annotate with GC content."""
    windowed = {}
    for wlen in window_lengths:
        win_bed = os.path.join(tmpdir, f'masked.windowed.{wlen}bp.bed')
        gc_bed = win_bed.replace('.bed', '.gc.bed')
        sp.call(f'bedtools makewindows -w {wlen} -s 20 -b {mask} > {win_bed}', shell=True)
        sp.call(
            f'bedtools nuc -fi {os.path.abspath(fasta)} -bed {win_bed} | cut -f 1-3,5 > {gc_bed}',
            shell=True,
        )
        windowed[wlen] = gc_bed
    return windowed


def _bin_gc_content(windowed):
    """Assign each window a GC bin and shuffle for random sampling."""
    binned = {}
    for wlen, gc_bed in windowed.items():
        binned_bed = gc_bed.replace('.bed', '.binned.bed')
        with open(gc_bed) as tsv, open(binned_bed, 'wt') as out:
            writer = csv.writer(out, delimiter='\t')
            for row in csv.reader(tsv, dialect='excel-tab'):
                gc_bin = get_gc_bin(row[3])
                if gc_bin != -9:
                    writer.writerow([row[0], row[1], row[2], gc_bin])

        lines = open(binned_bed).readlines()
        shuffled_bed = binned_bed.replace('.bed', '.shuffled.bed')
        with open(shuffled_bed, 'w') as out:
            out.writelines(random.sample(lines, len(lines)))
        binned[wlen] = shuffled_bed
    return binned


def _load_loci(meta, binned_loci, binned):
    """Populate *binned_loci* with up to n_windows entries per GC bin."""
    counts = {}
    for wlen, bed_file in binned.items():
        with open(bed_file) as fh:
            for line in fh:
                chrom, start, end, gc_bin = line.rstrip().split('\t')
                start, end, gc_bin = int(start), int(end), int(gc_bin)
                key = (chrom, gc_bin, wlen)
                if counts.get(key, 0) >= meta.n_windows:
                    continue
                counts[key] = counts.get(key, 0) + 1
                binned_loci[chrom].append((start, end, gc_bin, wlen))

    for contig in meta.user_contigs:
        binned_loci[contig].sort(key=itemgetter(0))


def _compute_contig_metadata(meta, bam, loci):
    """Compute per-contig DOC, mean read length, and insert-length statistics."""
    for contig in meta.user_contigs:
        welford = Welford()
        rl_sum = rc = bp_span = 0
        done = False

        for start, end in loci[contig]:
            if done:
                break
            for aln in bam.fetch(reference=contig, start=start, end=end):
                if aln.is_unmapped or not isinstance(aln.reference_length, int):
                    continue
                mid = alignment_midpoint(aln)
                if not (start <= mid <= end):
                    continue
                rc += 1
                rl_sum += aln.reference_length
                if abs(aln.template_length) <= meta.tlen_cap:
                    welford.update(abs(aln.template_length))
                if rc > meta.n_reads:
                    bp_span += aln.reference_end - start
                    done = True
                    break
            if not done:
                bp_span += end - start

        meta.read_len[contig] = rl_sum / rc if rc else 'nan'
        meta.doc[contig] = (meta.read_len[contig] * rc) / bp_span if bp_span else 'nan'
        meta.tlen[contig] = welford.mean
        meta.tlen_std[contig] = welford.std


def _compute_gc_metadata(meta, bam, binned_loci):
    """Compute mean and std of read counts per GC bin across the genome."""
    welfords = {}
    for contig in meta.user_contigs:
        for start, end, gc_bin, wlen in binned_loci[contig]:
            key = tuple2key((contig, wlen, gc_bin))
            if key not in welfords:
                welfords[key] = Welford()
            count = 0
            for aln in bam.fetch(reference=contig, start=start, end=end):
                if aln.reference_start is None or aln.reference_end is None:
                    continue
                if start <= alignment_midpoint(aln) <= end:
                    count += 1
            welfords[key].update(count)

    for key, w in welfords.items():
        meta.gc_rc[key] = w.mean
        meta.gc_std[key] = w.std
