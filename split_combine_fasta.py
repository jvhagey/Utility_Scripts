#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2022

from numpy import *
import glob
import os
import re
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CEND = '\033[0m'
# allowing more output
pd.set_option("max_columns", 20)
pd.set_option("max_rows", None)

def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Split_combine_fasta.py", description="""This will filter the inital SNP file to return a new file with only the important sites.""")
	parser.add_argument("-d", "--dir", dest="dir", action="store", default=None, required=True, help="Path for a directory where all fasta files are found. Expected to have an extention of: 'Important_Genes.fasta'")
	parser.add_argument("-e", "--extension", dest="extension", action="store", default="Important_Genes.fasta", required=False, help="Path for a directory where all fasta files are found. Expected to have an extention of: 'Important_Genes.fasta'")
	args = parser.parse_args()
	return args


def get_fastas(str_extension): 
	file_name_extension = "*-" + str(str_extension) #add "-" to the file name
	fasta_files = glob.glob(file_name_extension) # gather files
	return file_name_extension, fasta_files

def get_seq_names(fasta_file):
	seq_ids = []
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		seq_ids.append(seq_record.id)
	print(CRED + " This script assumes the sequences are in the same order in all files." + CEND)
	return seq_ids

def generate_fastas(file_name_extension, fasta_files, seq_ids):
	file_name_extension = file_name_extension.replace("*","")
	for seq_id in seq_ids:
		seq_name = re.search('ID=(.*);', seq_id)[0]
		seq_name = seq_name.replace("ID=", "")
		seq_name = seq_name.replace(";", "")
		new_records = []
		for fasta in fasta_files:
			sample_name = fasta.replace(file_name_extension,"")
			for seq_record in SeqIO.parse(fasta, "fasta"):
				if seq_record.id == seq_id:
					rec = SeqRecord(seq_record.seq, description="", id=sample_name)
					new_records.append(rec)
		new_fasta_name = seq_name + ".fasta"
		SeqIO.write(new_records, new_fasta_name, "fasta")
		print(CYEL + "The file {} was written.".format(new_fasta_name) + CEND)

def main():
	args = parse_cmdline()
	os.chdir(args.dir)
	#str_extension ="Important_Genes.fasta"
	str_extension = str(args.extension)
	file_name_extension, fasta_files = get_fastas(str_extension)
	seq_ids=get_seq_names(fasta_files[0])
	generate_fastas(file_name_extension, fasta_files, seq_ids)

if __name__ == '__main__':
	main()
