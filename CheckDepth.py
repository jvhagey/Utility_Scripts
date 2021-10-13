#!/usr/bin/env python

## Jill Hagey, PhD
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2021

import pandas as pd
import glob, os
from numpy import *
from argparse import ArgumentParser

# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CEND = '\033[0m'

def parse_cmdline():
	"""Parameters in scripts to parse sites below 25x and consensus alignment."""
	parser = ArgumentParser(prog="CheckDepth.py", description="""Checks the to see what percent of a sequence is above your depth cutoff""")
	parser.add_argument("-d", "--depth", dest="depth",action="store", default=None, required=True, help="""Choosen depth coverage for the cutoff for masking.""")
	parser.add_argument("-c", "--coverag-dir", type=str, dest="coverage_dir", action="store", required=True, default=None,
						help='A directory where files with per site depth coverage. Created by `bedtools genomecov -bga -split -ibam input.bam > coverage.txt`')
	args = parser.parse_args()
	return args

def calculate_good_depth(cutoff, coverage, chromosome):
	""" Read in the sites from bam file """
	with open(coverage,'r') as e:
		total = 0
		good_depth = 0
		bad_depth = 0
		for line in e:
			if chromosome in line:
				total = total + 1
				tempInfo =  line.rstrip("\n")
				info = tempInfo.split("\t")
				depth = int(info[3])
				if depth >= int(cutoff):
					good_depth = good_depth + 1
				else:
					bad_depth = bad_depth + 1
			else:
				pass
		percent_good = (good_depth / total) * 100
		percent_bad = (bad_depth / total) * 100
	return percent_good, percent_bad


def check_df(df,depth):
	df_bad = df[df["Percent_Good"] <= 70]
	if df_bad.empty:
		print(CYEL + " Yay, looks like all samples have at least 70% of each chromosome had a depth of {} at each nucleotide.".format(depth) + CEND)
	else:
		print(CRED + " A few samples have chromosomes with a less than 70% of their chromosome with a depth a depth of {} at each nucleotide.".format(depth) + CEND)
		print(df_bad)
		df_bad.to_csv('Depth_Bad_Stats.csv', sep='\t', index=True)


def main():
	args = parse_cmdline()
	#depth = 26
	#coverage = "17861_Calf-per-site-depth-coverage.txt"
	os.chdir(args.coverage_dir)
	dataframe = pd.DataFrame(columns=["Sample", "Chromosome", "Percent_Bad", "Percent_Good"])
	coverage_files = glob.glob("*-per-site-depth-coverage.txt")
	for coverage in coverage_files:
		file_name =  coverage.replace("-per-site-depth-coverage.txt", "")
		df = pd.read_csv(coverage, sep='\t', header=None)
		chromosomes = list(unique(df[0]))
		dic = {}
		for chromosome in chromosomes:
			percent_good, percent_bad = calculate_good_depth(args.depth, coverage, chromosome)
			dic1 = {chromosome: int(percent_bad)}
			dic.update(dic1)
			new_row = {"Sample": file_name,"Chromosome":chromosome, "Percent_Bad":percent_bad, "Percent_Good":percent_good}
			dataframe = dataframe.append(new_row, ignore_index = True)
		print("Finished with file {}.".format(coverage))
	check_df(dataframe, args.depth)
	dataframe.to_csv('Depth_stats.csv', sep='\t', index=True)


if __name__ == '__main__':
	main()