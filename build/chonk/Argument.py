#!/usr/bin/env python3
from argparse import RawTextHelpFormatter
import argparse,os,sys,json

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

	def metadata(self):
		"""
		subcommand to get the metadata for a bam file
		"""
		_help = """

    chonk metadata -i <bam> -e <mask> -f <fasta> -g <genome name> -r <chromosome> -o <output directory> -x <exclude region>
		
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
		parser.add_argument('-r',type=str,required=False,default=None)
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)
		if args.r == None:
			args.r = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')
			r_label = 'all'
		else:
			args.r = (args.r,)
			r_label = args.r
		if args.o == None:
			args.o = os.path.dirname(args.i)
		elif not os.path.exists(args.o):
			os.mkdir(args.o)
			if args.o[-1] == '/':
				args.o = args.o[:-1]
		elif args.o[-1] == '/' and os.path.exists(args.o):
			args.o = args.o[:-1]
			#args.o = '{}.{}.metadata.chonk.json'.format(
			#	args.i.replace('.bam','').replace('.sam','').replace('.cram',''),
			#	r_label
			#	)
		return args

	def breakpoints(self):
		"""
		subcommand breakpoints to discover SVs
		arguments for sv calling
		"""
		_help = """

    chonk breakpoints -m <meta> -g <genome name> -r <chromosome> -o <output directory>

		"""
		parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=_help, add_help=False)
		# meta file
		parser.add_argument('-m',type=str, required=True,default=None)
		parser.add_argument('-g',type=str, required=True,default=None)
		parser.add_argument('-r',type=str,required=False,default=None)
		parser.add_argument('-o',type=str,required=False,default=None)
		parser.add_argument('-h', '-help', required=False, action="store_true", default=False)
		args = parser.parse_args(sys.argv[2:])
		if args.h == True:
			sys.stderr.write(_help+'\n')
			sys.exit(1)
		with open(args.m) as json_file:
			data = json.load(json_file)
			metadata = data[0]
			meta_genome = metadata[0]
			meta_chrom = metadata[1]
			meta_bam = metadata[2]
		if args.g != meta_genome:
			sys.stderr.write('FATAL ERROR: metadata file does not contain genome '+args.g+' data'+'\n'+__usage__+'\n')
			sys.exit(1)
		if args.r == None:
			args.r = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')
			r_label = 'all'
		else:
			args.r = (args.r,)
			r_label = args.r
		for input_chrom in args.r:
			if not input_chrom in meta_chrom:
				sys.stderr.write('FATAL ERROR: chromosome ('+input_chrom+') not in metadata\n'+__usage__+'\n')
				sys.exit(1)
		if args.o == None:
			args.o = os.path.dirname(meta_bam)
		elif not os.path.exists(args.o):
			os.mkdir(args.o)
			if args.o[-1] == '/':
				args.o = args.o[:-1]
		elif args.o[-1] == '/' and os.path.exists(args.o):
			args.o = args.o[:-1]
		return args

	def features(self):
		"""
		subcommand to get the features for each sv call
		"""
		_help = """

    chonk features -i <bam> -e <mask> -f <fasta> -g <genome name> -r <chromosome> -o <output prefix>
		
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