#!/usr/bin/env python3
from chonk.Argument import Argument
from chonk.Breakpoint import breakpoints
from chonk.Metadata import metadataExtraction
def main():
	Args = Argument()
	if Args.command == 'breakpoints':
		breakpoints(Args.args)
	if Args.command == 'metadata':
		metadataExtraction(Args.args)

