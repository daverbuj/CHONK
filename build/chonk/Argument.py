#!/usr/bin/env python3
from argparse import RawTextHelpFormatter
import argparse,os,sys

__version__='0.0.1'

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

class Argument():

	def __init__(self):
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=__usage__, add_help=False)
		parser.add_argument('subcommand', help='subcommand to run')
		args = parser.parse_args(sys.argv[1:2])
		if not hasattr(self, args.subcommand):
			sys.stderr.write('FATAL ERROR: unrecognized command\n'+__usage__+'\n')
			sys.exit(1)
		self.command = args.subcommand
		self.args = getattr(self, args.subcommand)()

	def breakpoints(self):
		"""
		subcommand breakpoints to discover SVs
		arguments for sv calling
		"""
		_help = """

    chonk breakpoints -i <bam> -m <meta> -o <output prefix.>

		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# bam file (This is not necesarry since bam file could be included in metadata)
		parser.add_argument('-i',type=str, required=True,default=None)
		# meta file
		parser.add_argument('-m',type=str, required=True,default=None)
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)
		if args.o == None:
			args.o = '{}.{}.chonk.sv.bed'.format(
				args.i.replace('.bam','').replace('.sam','').replace('.cram',''),
				args.r
				)
		return args

	def metadata(self):
		"""
		subcommand to get the metadata for a bam file

		"""
		_help = """

    chonk metadata -i <bam> -e <mask> -f <fasta> -g <genome name> -r <chromosome> -o <output prefix> -x <exclude region>
		
		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# bam file
		parser.add_argument('-i',type=str, required=True,default=None)
		# mask .bed file
		parser.add_argument('-e',type=str, required=True,default=None)
		# fasta
		parser.add_argument('-f',type=str, required=True,default=None)
		# genome name
		parser.add_argument('-g',type=str, required=True,default=None)
		# chromosome
		parser.add_argument('-r',type=str,required=True,default=None)
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)
		if args.o == None:
			args.o = '{}.{}.metadata.chonk.json'.format(
				args.i.replace('.bam','').replace('.sam','').replace('.cram',''),
				args.r
				)
		return args