#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2019
## given an fasta file with repeated names this script outputs a new fasta file with names numbered

#importing packages
import os
import re
import csv
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="fasta_rename.py", description="Takes a .fasta file and adds numbers to end of identical sequence names")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were at .fasta files are found (required)")
    parser.add_argument("-f", "--infile", dest="infile", action="store", default=None, required=True, help="Input file name were at .fasta file is found (required)")
    parser.add_argument("-o", "--outfile", dest="outfile", action="store", default=None, required=True, help="Out file name for new .fasta file with names changed (required)")
    return parser.parse_args()

def main():
    #Parse command-line
    args = parse_cmdline()
    #Calling in files
    os.chdir(args.indirname) #change directory
    file = args.infile
    open(args.outfile, 'a').close() #making empty text file for later
    with open(file) as f:
        seen = set()
        counting = list()
        count = 0
        for line in f:
            if line.startswith('>'):# grab name line
                count = count + 1
                line = re.sub(r'\n', '', line, flags=re.MULTILINE) #grab name line remove \n at end of each line
                if line in seen:
                    counting.append(line) #add line to counting list
                    num = counting.count(line) #count the number of times you see line in counting list
                    new_line = line + "_" + str(num) #make new name with the number on end
                    with open(args.outfile, "a") as f2:
                        #write line to output file
                        f2.write(new_line)
                        f2.write("\n")
                else:
                    seen.add(line)
                    counting.append(line)
                    new_line1 = line + "_1" #make new name with the number on end
                    with open(args.outfile, "a") as f3:
                        f3.write(new_line1)
                        f3.write("\n")
            else:
                with open(args.outfile, "a") as f4:
                    f4.write(line)
        print('\033[93m' + '\nThere are ' + str(count) + ' sequenes in the new file named ' + args.outfile + '\n')

if __name__ == '__main__':
    main()
