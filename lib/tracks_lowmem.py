# -*- coding: UTF-8 -*-
import pysam
from lib.log import *
from lib.buffer import *
import numpy as np


class PhysicalCoverageTrack:
    name = "Physical coverage (low memory)"

    def __init__(self, _next, genome_length, buffer_length=1000000):
        self.next = _next
        self.track = PrefixSumBuffer(genome_length, buffer_length)

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads align correctly
            and read.is_paired             # SAM flag 0x1
            and read.is_proper_pair        # SAM flag 0x2
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
        ):
            for value in self.track.pop_until(read.reference_start):
                self.next.process(value)

            self.track.add_interval(
                read.reference_start,
                read.reference_start + read.template_length
            )

    def close(self):
        for value in self.track.pop_all():
            self.next.process(value)


class SequenceCoverageTrack:
    name = "Sequence coverage (low memory)"

    def __init__(self, _next, genome_length, buffer_length=1000000):
        self.next = _next
        self.track = PrefixSumBuffer(genome_length, buffer_length)

    def process(self, read):
        if (
            # The read is mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
        ):
            for value in self.track.pop_until(read.reference_start):
                self.next.process(value)

            self.track.add_interval(
                read.reference_start,
                read.reference_end
            )

    def close(self):
        for value in self.track.pop_all():
            self.next.process(value)
