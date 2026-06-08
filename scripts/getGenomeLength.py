#!/bin/env python3

import argparse
import os
import sys
import gzip
from   collections import defaultdict


def get_arguments() -> str:
    """get the arguments"""

    d = "Get the length of a genome fasta file"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="Genome in fasta file", required=True)

    args  = parser.parse_args()
    fasta = args.fasta

    assert os.path.isfile(fasta), f"Could not locate {fasta}"

    return fasta

def countBases(fasta: str) -> None:
    """count the number of bases in the assembly"""

    fh = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')

    bases = 0

    for line in fh:
        if (len(line) == 0) or (line[0] == '>') or (line[0] == '#'):
            continue
        bases += len(line.strip('\n'))

    fh.close()

    ofh = open("GenomeSize.txt", 'w')

    ofh.write(f"Total bases: {bases}\n")

    ofh.close()

def main() -> int:
    """entry points to the small application"""

    # get arguments
    fasta = get_arguments()

    # parse the fasta file for the total bp length
    countBases(fasta)

    return 0

if __name__ == "__main__":
    main()
