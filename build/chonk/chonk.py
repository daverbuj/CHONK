#!/usr/bin/env python3
from chonk.Argument import Argument
from chonk.Breakpoint import breakpoint
def main():
	Args = Argument()
	if Args.command == 'breakpoints':
		breakpoint(Args.args)