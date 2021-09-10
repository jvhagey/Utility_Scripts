#!/usr/bin/env python

## Jill Hagey, PhD
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2021

from argparse import ArgumentParser

def parse_cmdline():
	"""This will converte a wrapped fasta file to a one line fasta file."""
	parser = ArgumentParser(prog="MasksSites.py", description="""Given an input fasta file this will convert the fasta to a one-line fasta.""")
	parser.add_argument("-f", "--fasta", dest="fastafile", action="store", default=None, required=True, help="""Fasta file that will converted""")
	parser.add_argument("-o", "--outfile", dest="outfile",action="store", default=None, required=True, help="""File name for output to be written to.""")
	args = parser.parse_args()
	return args

def read_in(fastaFile):
	final = {}
	with open(fastaFile, "r") as file_in:
		for line in file_in:
			line = line.strip()
			if line.startswith(">"):
				id = line
				final[id] = " "
			else:
				final[id] += line
	return final

def write_out(final, outfile):
	with open(outfile, "a") as f:
		for key, val in final.items():
			f.write(key)
			f.write("\n")
			f.write(val)
			f.write("\n")

def main():
	args = parse_cmdline()
	final = read_in(args.fastafile)
	write_out(final, args.outfile)

if __name__ == '__main__':
	main()
