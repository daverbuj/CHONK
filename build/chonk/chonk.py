#!/usr/bin/env python3
from chonk.Argument import Argument
from chonk.Metadata import metadata
from chonk.Breakpoint import breakpoints
from chonk.Features import features
import random

def main():
	Args = Argument()
	if Args.command == 'metadata':
		# init the random seed
		random.seed(Args.args.s)

		metadata(Args.args)


	if Args.command == 'breakpoints':
		breakpoints(Args.args)


	if Args.command == 'features':
		features(Args.args)