#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## Given a .abi file from sanger sequencing, script will convert to a .fasta file


#https://biopython.org/wiki/Converting_sequence_files
#https://biopython.org/wiki/SeqIO

#importing packages
import glob,os
from Bio import SeqIO
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="abi_to_fasta.py", description="Converts .ab1 files into .fasta files")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .fasta files are found (required)")
    return parser.parse_args()

if __name__ == '__main__':
    #Parse command-line
    args = parse_cmdline()
    #Calling in files
    os.chdir(args.indirname) #change directory
    infiles = glob.glob("*.ab1")
    #creating for loop to call in .ab1 files
    count = 0
    for file in infiles:
        with open(file, "rb") as input_handle:
            #Making file names for the new fasta files
            filename = input_handle.name + '.fasta'
            filename = filename.replace('.ab1','')
            with open(filename, "w") as output_handle:
                sequences = SeqIO.parse(input_handle, "abi")
                SeqIO.write(sequences, output_handle, "fasta")
        count = count + 1
print("Converted %i .fasta files"  % count)
