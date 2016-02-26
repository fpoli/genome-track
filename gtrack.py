#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Generate WIG tracks from reads aligned on a reference genome.
"""

import sys
import argparse
import pysam
from lib.log import printlog
from lib.fasta import read_fasta_info
from lib.pipeline import *
from lib.tracks import *
import lib.tracks_lowmem as lowmem

track_classes = {
    "physical_coverage": PhysicalCoverageTrack,
    "sequence_coverage": SequenceCoverageTrack,
    "average_insert": AverageInsertTrack,
    "fr_mates": FRMatesTrack,
    "rf_mates": RFMatesTrack,
    "ff_mates": FFMatesTrack,
    "rr_mates": RRMatesTrack,
    "ffrr_mates": FFRRMatesTrack,
    "single_mate": SingleMateTrack,
    "multiple_mapping": MultipleMappingTrack,
    "hs_cigar": HSCigarTrack,
    "physical_coverage_lowmem": lowmem.PhysicalCoverageTrack,
    "sequence_coverage_lowmem": lowmem.SequenceCoverageTrack,
}

# Define the command line arguments
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "track_name",
    choices=track_classes.keys(),
    help="Name of the track to be computed."
)
parser.add_argument(
    "-i", "--input",
    type=argparse.FileType("r"),
    required=True,
    help="SAM/BAM input file."
)
parser.add_argument(
    "-o", "--output",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="WIG output file. The default is stdout."
)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "-l", "--length",
    type=int,
    help="Length of the reference genome."
)
group.add_argument(
    "-g", "--genome",
    help="Fasta file of the reference genome."
)
parser.add_argument(
    "--min",
    type=int,
    help="Inserts below this value will be ignored."
)
parser.add_argument(
    "--max",
    type=int,
    help="Inserts above this value will be ignored."
)

# Check number of arguments
if len(sys.argv) < 2:
    parser.print_help(sys.stderr)
    sys.exit(1)

# Parse arguments
args = parser.parse_args()

# Read reference genome length
if args.length is not None:
    genome_length = args.length
else:
    printlog("(*) Reading reference genome length...")
    fasta_info = read_fasta_info(args.genome)
    if len(fasta_info) == 1:
        genome_length = fasta_info[0][1]
        printlog("(i) Reference genome length is {}".format(genome_length))
    else:
        printlog("The fasta reference genome must contain exactly 1 sequence.")
        sys.exit(1)

# Open files and pipeline
Track = track_classes[args.track_name]
samfile = pysam.AlignmentFile(args.input)

printlog("(i) Input file: {}".format(args.input.name))
printlog("(i) Filter insert min: {}, max: {}".format(args.min, args.max))
printlog("(i) Selected track: {}".format(Track.name))
printlog("(i) Output file: {}".format(args.output.name))

if args.min is None and args.max is None:
    track_name = Track.name
else:
    track_name = "{} (insert min: {}, max: {})".format(
        Track.name, args.min, args.max
    )

wigfile = WigFileWriter(args.output, name=track_name)
track = Track(wigfile, genome_length)
pipeline = track

if not (args.min is None and args.max is None):
    filter_len = InsertLengthFilter(track, args.min, args.max)
    pipeline = filter_len

# Process reads
printlog("(*) Processing reads...")
for read in samfile:
    pipeline.process(read)

# Close filepipeline
printlog("(*) Closing pipeline and files...")
pipeline.close()
samfile.close()
