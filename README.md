# Utility_Scripts
1. abi_to_fasta.py - converts abi files in a directory from sanger sequencing to fasta files.
2. Multiplot.R - is a function a that puts together multiple plots into the same image. Use together with ggsave.
3. fasta_rename.py - if you have a fasta file with repeated sequence names this script outputs a new fasta file with repeated names numbered.
4. Unwrap_fasta.py - will convert a wrapped fasta file to a one line fasta file.
5. MaskSites.py - Given one line fasta file and a coverage file this will mask (replace with N) sites with a depth coverage below the choosen level.
6. CheckDepth.py - Checks the to see what percent of a sequence is above your depth cutoff.
