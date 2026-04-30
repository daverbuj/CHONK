#!/usr/bin/env python3
"""Command-line interface for CHONK."""
import argparse
import os
import sys
from argparse import RawTextHelpFormatter

from chonk import __version__
from chonk.utils import tokenize_contigs

_LOGO = r"""
  oooooooo8   ooooo ooooo    ooooooo    oooo   oooo   oooo   oooo
o888     88    888   888   o888   888o   8888o  88     888  o88
888            888ooo888   888     888   88 888o88     888888
888o     oo    888   888   888o   o888   88   8888     888  88o
 888oooo88    o888o o888o    88ooo88    o88o    88    o888o o888o

       detect germline and somatic structural variants
---------------------------------------------------------------
Version {version}
Authors: Danny Antaki  ·  Dan Averbuj
"""


def _outdir(path, fallback):
    """Return the output directory, defaulting to the directory of *fallback*."""
    if path is None:
        path = os.path.dirname(os.path.abspath(fallback))
    os.makedirs(path, exist_ok=True)
    return path if path.endswith('/') else path + '/'


class CLI:
    """Top-level argument dispatcher for the chonk command."""

    _SUBCOMMANDS = ('metadata', 'breakpoints', 'features')

    def __init__(self):
        usage = _LOGO.format(version=__version__) + (
            '\n  chonk <subcommand> [options]\n\n'
            '  subcommands:\n'
            '    metadata     compute per-sample sequencing metadata\n'
            '    breakpoints  detect SV breakpoints\n'
            '    features     extract SV features\n\n'
            '  -h  show this message and exit\n'
        )
        parser = argparse.ArgumentParser(
            formatter_class=RawTextHelpFormatter,
            usage=usage,
            add_help=False,
        )
        parser.add_argument('subcommand', choices=self._SUBCOMMANDS)
        args = parser.parse_args(sys.argv[1:2])
        self.args = getattr(self, args.subcommand)()

    # ---------------------------------------------------------------- metadata

    def metadata(self):
        usage = """
  chonk metadata -i BAM -f FASTA [-r CONTIGS] [-x EXCLUDE] [-s SEED] [-o OUTDIR]

  Compute per-sample metadata (depth of coverage, insert-length statistics,
  GC-binned read counts) required by the breakpoints and features subcommands.

Required:
  -i   BAM alignment file (must be indexed)
  -f   reference FASTA file (must be indexed with samtools faidx)

Optional:
  -r   comma-separated list of contigs to process  [all]
  -x   BED file of regions to exclude              [none]
  -s   random seed for region sampling             [42]
  -o   output directory                            [BAM directory]
"""
        p = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                    usage=usage, add_help=False)
        p.add_argument('-i', required=True)
        p.add_argument('-f', required=True)
        p.add_argument('-r', default=None)
        p.add_argument('-x', default=None)
        p.add_argument('-s', type=int, default=42)
        p.add_argument('-o', default=None)
        p.add_argument('-h', action='store_true', default=False)
        args = p.parse_args(sys.argv[2:])
        if args.h:
            sys.stderr.write(usage + '\n')
            sys.exit(0)
        args.o = _outdir(args.o, args.i)
        if args.r is not None:
            args.r = tokenize_contigs(args.r)
        return args

    # ------------------------------------------------------------- breakpoints

    def breakpoints(self):
        usage = """
  chonk breakpoints -m METADATA [-r CONTIGS] [-o OUTDIR]

  Call SV breakpoints using split-reads and discordant paired-end reads.

Required:
  -m   metadata JSON produced by 'chonk metadata'

Optional:
  -r   comma-separated list of contigs to process  [all from metadata]
  -o   output directory                            [metadata directory]
"""
        p = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                    usage=usage, add_help=False)
        p.add_argument('-m', required=True)
        p.add_argument('-r', default=None)
        p.add_argument('-o', default=None)
        p.add_argument('-h', action='store_true', default=False)
        args = p.parse_args(sys.argv[2:])
        if args.h:
            sys.stderr.write(usage + '\n')
            sys.exit(0)
        args.o = _outdir(args.o, args.m)
        if args.r is not None:
            args.r = tokenize_contigs(args.r)
        return args

    # --------------------------------------------------------------- features

    def features(self):
        usage = """
  chonk features -bp BED -metadir DIR -fasta FASTA -k K -rlen RLEN -pk PK -o OUT
                 [-rm RM_BED] [-sd SD_BED] [-str STR_BED]

  Extract SV features (coverage, supporting fragments, k-mers, sequence context,
  and optionally genomic-context overlap) from a breakpoint BED file.

Required:
  -bp       breakpoint BED file (tab-separated: chrom start end svtype cipos ciend iid gt)
  -metadir  directory containing per-sample metadata JSON files
  -fasta    reference FASTA file
  -k        k-mer size for junction k-mer features
  -rlen     read-length modifier (float, e.g. 1.0)
  -pk       k-mer size for pseudo-alignment features
  -o        output file prefix (must end in .txt)

Optional overlap tracks (all three required to enable overlap features):
  -rm       repeat-masker merged BED
  -sd       segmental-duplication merged BED
  -str      short tandem repeat merged BED
"""
        p = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                    usage=usage, add_help=False)
        p.add_argument('-bp', required=True)
        p.add_argument('-metadir', required=True)
        p.add_argument('-fasta', required=True)
        p.add_argument('-k', type=int, required=True)
        p.add_argument('-rlen', type=float, required=True)
        p.add_argument('-pk', type=int, required=True)
        p.add_argument('-o', required=True)
        p.add_argument('-rm', default=None)
        p.add_argument('-sd', default=None)
        p.add_argument('-str', dest='str_bed', default=None)
        p.add_argument('-h', action='store_true', default=False)
        args = p.parse_args(sys.argv[2:])
        if args.h:
            sys.stderr.write(usage + '\n')
            sys.exit(0)
        if not args.metadir.endswith('/'):
            args.metadir += '/'
        return args


# --------------------------------------------------------------------------- #
#  Entry point                                                                  #
# --------------------------------------------------------------------------- #

def main():
    """Main entry point dispatched by the ``chonk`` console script."""
    import random
    cli = CLI()

    if sys.argv[1] == 'metadata':
        random.seed(cli.args.s)
        from chonk.metadata import run_metadata
        run_metadata(cli.args)

    elif sys.argv[1] == 'breakpoints':
        from chonk.breakpoints import run_breakpoints
        run_breakpoints(cli.args)

    elif sys.argv[1] == 'features':
        from chonk.features.extract import run_features
        run_features(cli.args)


if __name__ == '__main__':
    main()
