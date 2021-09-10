#!/usr/bin/env python

## Jill Hagey, PhD
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2021
## Original Author: Shatavia Morrison @ cdc
## https://github.com/SMorrison42/legionella_pneumophila_genomics/blob/master/scripts/maskSites.py
## Converted to python3 with some added bells and whistles...

from numpy import *
from Bio import SeqIO
import pandas as pd
from argparse import ArgumentParser

def parse_cmdline():
	"""Parameters in scripts to parse sites below 25x and consensus alignment."""
	parser = ArgumentParser(prog="MasksSites.py", description="""Given one line fasta file and a coverage file this will mask (replace with N) sites with a depth coverage below the choosen level """)
	parser.add_argument("-f", "--file", dest="fastaFile", action="store", default=None, required=True, help="""fasta file that will checked for sites to mask""")
	parser.add_argument("-o", "--outfile", dest="outputFile",action="store", default=None, required=True, help="""File name for output to be written to.""")
	parser.add_argument("-d", "--depth", dest="depth",action="store", default=None, required=True, help="""Choosen depth coverag for the cutoff for masking.""")
	parser.add_argument("-c", "--coverage", type=str, dest="coverage", action="store", required=True, default=None,
						help="""A file with per site depth coverage. Created by `bedtools genomecov -bga -split -ibam input.bam > coverage.txt`""")
	args = parser.parse_args()
	return args


def replaceNucs(seqString,coorList):
	"""The sites that need to be changed from nucleotide to N"""
	seq = list(seqString)
	for v in range(0,len(coorList)):
		rangeLen = coorList[v][1] - coorList[v][0]
		for m in range(0, rangeLen):
			pos = coorList[v][0] + m
			seq[pos] ='N'
	return seq

def read_fasta(fastaFile):
	"""Read fasta file with BioPython to read into temp directory."""
	tempHold=[]
	sequence = list(SeqIO.parse(fastaFile, "fasta"))
	for i in range(0,len(sequence)):
		tempHold.append(sequence[i].id)
		tempHold.append(sequence[i].seq)
	return tempHold

def extract_sites_from_coverage(cutoff, coverage, chromosome):
	""" Read in the sites from bam file """
	sitesChange=[]
	with open(coverage,'r') as e:
		for line in e:
			if chromosome in line:
				coorInfo=[]
				tempInfo =  line.rstrip("\n")
				info = tempInfo.split("\t")
				start = int(info[1])
				end = int(info[2])
				depth = int(info[3])
				if depth < int(cutoff):
					coorInfo.append(start)
					coorInfo.append(end)
					sitesChange.append(coorInfo) # contains list of start and stop position for locations with low coverage
			else:
				pass
	return sitesChange

def create_outfile(tempHold, sitesChange, outputFile, chromosome, count ,only_odd):
	if count == 0:
		updateSeq = replaceNucs(tempHold[1], sitesChange)
		finalSeq = ''.join(updateSeq)
		f = open(outputFile,"w")
		f.write(chromosome +"\n")
		f.write(finalSeq)
	elif count >= 1:
		updateSeq = replaceNucs(tempHold[only_odd], sitesChange)
		finalSeq = ''.join(updateSeq)
		with open(outputFile, "a") as f:
			f.write(chromosome + "\n")
			f.write(finalSeq)

def main():
	args = parse_cmdline()
	count = 0
	df = pd.read_csv(args.coverage, sep='\t', header=None)
	chromosomes = list(unique(df[0]))
	only_odd = list(range(1, len(chromosomes) * 2, 2)) #make a list of odd numbers to get sequences later
	tempHold = read_fasta(args.fastaFile)
	for chromosome in chromosomes:
		odd_count = only_odd[count]
		sitesChange = extract_sites_from_coverage(args.depth, args.coverage, chromosome)
		create_outfile(tempHold, sitesChange, args.outputFile, chromosome, count, odd_count)
		count = count + 1
		print("Finished with sequence {}.".format(chromosome))


if __name__ == '__main__':
	main()
