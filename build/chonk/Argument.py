#!/usr/bin/env python3
from chonk.Backend import tokenize_user_contigs
from argparse import RawTextHelpFormatter
import argparse,os,sys,json

__version__='0.0.2'

__usage__="""

  oooooooo8   ooooo ooooo    ooooooo    oooo   oooo   oooo   oooo 
o888     88    888   888   o888   888o   8888o  88     888  o88   
888            888ooo888   888     888   88 888o88     888888     
888o     oo    888   888   888o   o888   88   8888     888  88o   
 888oooo88    o888o o888o    88ooo88    o88o    88    o888o o888o 

          detect germline and somatic structural variants 
------------------------------------------------------------------
Version {}
Authors: Danny Antaki dantaki at ucsd dot edu
         Dan Averbuj 

  chonk
         subcommands:

         <breakpoints>  <metadata>  <features>   
         <genotype>     <somatic>


optional arguments:

  -h        show this message and exit

""".format(__version__)

class Argument(object):

	def __init__(self):
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=__usage__, add_help=False)
		parser.add_argument('subcommand', help='subcommand to run')
		args = parser.parse_args(sys.argv[1:2])
		if not hasattr(self, args.subcommand):
			sys.stderr.write('FATAL ERROR: unrecognized command\n'+__usage__+'\n')
			sys.exit(1)
		self.command = args.subcommand
		self.args = getattr(self, args.subcommand)()

	def metadata(self):
		"""
		subcommand to get the metadata for a bam file
		"""
		_help = """

    chonk metadata <-if> [-rxso]

about: Compute metadata for SV calling 

Required arguments:
    -i        bam file
    -f        fasta file

Optional arguments:
    -r        restrict to comma-separated list of contigs
    -x        bed file of regions to exclude
    -s        random seed [default: 42]
    -o        output directory   
		
		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# bam file
		parser.add_argument('-i',type=str, required=True,default=None)
		# fasta
		parser.add_argument('-f',type=str, required=True,default=None)
		# contigs
		parser.add_argument('-r',type=str,required=False,default=None)
		# mask .bed file
		parser.add_argument('-x',type=str, required=False,default=None)
		# random seed
		parser.add_argument('-s',type=int, required=False,default=42) 
		# output directory
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])

		## if help is defined, display usage and exit
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)
		
		## if the output directory is None, use the same directory as the BAM file
		if args.o == None:
			args.o = os.path.dirname(args.i)
		# add "/" to path if it does not end with "/""
		if not args.o.endswith('/'): args.o += '/'
		# make output directory if it does not exist
		if not os.path.exists(args.o):
			os.mkdir(args.o)
		
		## split contigs into a tuple
		if args.r !=None:
			args.r = tokenize_user_contigs(args.r)

		return args

	def breakpoints(self):
		"""
		subcommand breakpoints to discover SVs
		arguments for sv calling
		"""
		_help = """

	chonk breakpoints <-m> [-ro]

about: Call SV breakpoints 

Required arguments:
    -m        metadata file

Optional arguments:
    -r        restrict to comma-separated list of contigs
    -o        output directory

		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# metadata file
		parser.add_argument('-m',type=str, required=True,default=None)
		# contigs
		parser.add_argument('-r',type=str,required=False,default=None)
		# output directory
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])

		## if help is defined, display usage and exit
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)

		## if the output directory is None, use the same directory as the metadata file
		if args.o == None:
			args.o = os.path.dirname(args.m)
		# add "/" to path if it does not end with "/"
		if not args.o.endswith('/'): args.o += '/'
		# make output directory if it does not exist
		if not os.path.exists(args.o):
			os.mkdir(args.o)

		## split contigs into a tuple
		if args.r != None:
			args.r = tokenize_user_contigs(args.r)
		
		return args

	def features(self):
		"""
		subcommand to get the features for each sv call
		"""
		_help = """

	chonk features <-mb> [-ro]

about: Call SV breakpoints 

Required arguments:
    -m        metadata file
    -b        breakpoint file

Optional arguments:
    -r        restrict to comma-separated list of contigs
    -o        output directory
		
		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# metadata file
		parser.add_argument('-m',type=str, required=True,default=None)
		# breakpoint file
		parser.add_argument('-b',type=str, required=True,default=None)
		# contigs
		parser.add_argument('-r',type=str,required=False,default=None)
		# output directory
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])

		## if help is defined, display usage and exit
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)

		## if the output directory is None, use the same directory as the breakpoint file
		if args.o == None:
			args.o = os.path.dirname(args.b)
		# add "/" to path if it does not end with "/"
		if not args.o.endswith('/'): args.o += '/'
		# make output directory if it does not exist
		if not os.path.exists(args.o):
			os.mkdir(args.o)

		## split contigs into a tuple
		if args.r != None:
			args.r = tokenize_user_contigs(args.r)

		return args