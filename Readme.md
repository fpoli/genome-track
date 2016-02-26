Genome-track
============

Command line program to generate WIG tracks (sequence coverage, physical coverage, average insert length, ...) from a SAM/BAM file.

Academic project for the Bioinformatics course at University of Padova. Written in Python 3 using PySAM.


Quick start
-----------

1. Check to have a Python 3 interpreter `python3 --version`
2. Install PySAM `pip3 install pysam`
3. Generate the sequence coverage
`./gtrack.py sequence_coverage -g genome.fasta -i alignments.sam -o sequence_coverage.wig`


Examples
--------

```
./gtrack.py sequence_coverage -g genome.fasta -i alignments.sam -o sequence_coverage.wig
./gtrack.py sequence_coverage -l 123456789000 -i alignments.sam -o sequence_coverage.wig
./gtrack.py physical_coverage -g genome.fasta -i alignments.sam > physical_coverage.wig
./gtrack.py physical_coverage -g genome.fasta -i alignments.sam -o physical_coverage.wig
./gtrack.py average_insert -g genome.fasta -i alignments.sam -o physical_coverage.wig
./gtrack.py average_insert --min 100 --max 50000 -g genome.fasta -i alignments.sam -o physical_coverage.wig
...
```


Usage
-----

```
usage: gtrack.py [-h] -i INPUT [-o OUTPUT] (-l LENGTH | -g GENOME) [--min MIN]
                 [--max MAX]
                 {sequence_coverage,rf_mates,hs_cigar,average_insert,ff_mates,fr_mates,single_mate,rr_mates,ffrr_mates,multiple_mapping,physical_coverage_lowmem,sequence_coverage_lowmem,physical_coverage}

Generate WIG tracks from reads aligned on a reference genome.

positional arguments:
  {sequence_coverage,rf_mates,hs_cigar,average_insert,ff_mates,fr_mates,single_mate,rr_mates,ffrr_mates,multiple_mapping,physical_coverage_lowmem,sequence_coverage_lowmem,physical_coverage}
                        Name of the track to be computed.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        SAM/BAM input file.
  -o OUTPUT, --output OUTPUT
                        WIG output file. The default is stdout.
  -l LENGTH, --length LENGTH
                        Length of the reference genome.
  -g GENOME, --genome GENOME
                        Fasta file of the reference genome.
  --min MIN             Inserts below this value will be ignored.
  --max MAX             Inserts above this value will be ignored.
```


License
-------

Genome-track generates WIG tracks from reads aligned on a reference genome.

Copyright (C) 2016 Federico Poli <federpoli@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
