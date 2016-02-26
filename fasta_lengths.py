#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
This script prints the id and length of each sequence stored in a fasta file.
"""

import sys
import argparse
from lib.fasta import read_fasta_info

# Define the command line arguments
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "fasta_file",
    help="the fasta file"
)

# Check number of arguments
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# Parse arguments
args = parser.parse_args()

# Print info
for seq_id, seq_length in read_fasta_info(args.fasta_file):
    print("{}: {}".format(seq_id, seq_length))
