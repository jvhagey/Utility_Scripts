#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2020
## Given Project_dataset.txt and input directory it will produce a with host name, reads, average read lengths and number of reads in a fastq file.

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import glob,os
import re
import gzip


def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Get_Sequence_Files.py", description="""Given Project_dataset.txt and input directory it will produce a with host name, 
     reads, average read lengths and number of reads in a fastq file""")
	parser.add_argument("-f", "--file", dest="infile", action="store", default=None, required=True, help="""Full location of Project_dataset.txt file. Include the file name in the path. Needs to be tab delimited (required)""")
	#parser.add_argument("-o", "--outdir", dest="outdir",action="store", default=None, required=True, help="File name for output to be written to.")
	parser.add_argument("-i", "--indir", type=str, dest="indir", action="store", required=True, default=None,
						help='A file with seqID on each line (required)')
	args = parser.parse_args()
	return args


def get_host(seqID, infile):
	df = pd.read_csv(infile, sep='\t', header=0, dtype='unicode')  # Has sequence ID and path information
	host = df.loc[df['SampleName'] == str(seqID), 'Host'].iloc[0]
	return host


def get_seq_length(indir, infile):
	os.chdir(indir)  # change directory
	# files with extention that will be looped through code
	fastq_files = glob.glob("*.fastq.gz")
	df = pd.DataFrame(columns=['SeqID', 'Lane', 'Host', 'ave_seq_length', "num_seqs"])  # make a empty dataframe
	for file in fastq_files:
		if "R2" in file:
			pass
		else:
			if "L001" in file:
				Lane = "L001"
			if "L002" in file:
				Lane = "L002"
			else:
				Lane = "L003"
			seqID = re.search("^[^_]*", file).group(0)  # Capture the Sequence ID at the beginning
			num_seqs = 0
			Host = get_host(seqID, infile).capitalize()
			lengths = []
			with gzip.open(file, "r") as f:
				lines = f.readlines()
				for line in lines:
					line = line.decode("utf8").strip('\n')
					if line.startswith(('A', "G", "C", "T")):
						num_seqs = num_seqs + 1
						lengths.append(len(line))
					else:
						pass
			ave_seq_length = int(sum(lengths)/len(lengths))
			new_row = {'SeqID': seqID, 'Lane': Lane, 'Host': Host, 'ave_seq_length': ave_seq_length, "num_seqs": num_seqs}
			new_row["Theoretical_coverage"] = (new_row["num_seqs"]*new_row["ave_seq_length"])/(9.2*1000000)
			df = df.append(new_row, ignore_index=True)
			print(new_row)
	df.to_csv('read_lengths.csv', sep=',', index=False)


def main():
	args = parse_cmdline()
	get_seq_length(args.indir, args.infile)


if __name__ == '__main__':
	main()
