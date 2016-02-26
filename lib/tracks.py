# -*- coding: UTF-8 -*-
import pysam
from lib.log import *
import numpy as np


class PhysicalCoverageTrack:
    name = "Physical coverage"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

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
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class SequenceCoverageTrack:
    name = "Sequence coverage"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            # The read is mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_end] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class AverageInsertTrack:
    name = "Average insert length"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff_sum = [0] * genome_length
        self.diff_num = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
        ):
            template_start = read.reference_start
            template_end = template_start + read.template_length
            self.diff_sum[template_start] += read.template_length
            self.diff_sum[template_end] -= read.template_length
            self.diff_num[template_start] += 1
            self.diff_num[template_end] -= 1

    def close(self):
        current_sum = 0
        current_num = 0
        for i in range(self.genome_length):
            current_sum += self.diff_sum[i]
            current_num += self.diff_num[i]
            if current_num == 0:
                average_inser_length = 0
            else:
                average_inser_length = current_sum / current_num
            self.next.process(average_inser_length)


class FRMatesTrack:
    name = "Forward-reverse mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # Forward-reverse mates
            and not read.is_reverse        # SAM flag 0x10
            and read.mate_is_reverse       # SAM flag 0x20
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class RFMatesTrack:
    name = "Reverse-forward mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # Forward-reverse mates
            and read.is_reverse            # SAM flag 0x10
            and not read.mate_is_reverse   # SAM flag 0x20
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class FFMatesTrack:
    name = "Forward-forward mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # Forward-forward mates
            and not read.is_reverse        # SAM flag 0x10
            and not read.mate_is_reverse   # SAM flag 0x20
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class RRMatesTrack:
    name = "Reverse-reverse mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # Reverse-reverse mates
            and read.is_reverse        # SAM flag 0x10
            and read.mate_is_reverse   # SAM flag 0x20
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class FFRRMatesTrack:
    name = "FF+RR mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            read.template_length > 0
            # Both reads are mapped
            and read.is_paired             # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # Forward-forward or reverse-reverse mates
            and read.is_reverse == read.mate_is_reverse  # SAM flags 0x10, 0x20
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_start + read.template_length] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class SingleMateTrack:
    name = "Single mates"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            # The read is mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # The mate is unmapped
            and read.mate_is_unmapped
        ):
            # TODO: the peak should be where the other mate was expected
            self.diff[read.reference_start] += 1
            self.diff[read.reference_end] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class MultipleMappingTrack:
    name = "Multiple mapping reads"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            # The read is mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            # This is one of multiple alignments
            and read.is_secondary          # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_end] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)


class HSCigarTrack:
    name = "Reads with H or S in the CIGAR"

    def __init__(self, _next, genome_length):
        self.next = _next
        self.genome_length = genome_length
        self.diff = [0] * genome_length

    def process(self, read):
        if (
            # The read is mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            # Remove multiple and chimeric alignments
            and not read.is_secondary      # SAM flag 0x100
            and not read.is_supplementary  # SAM flag 0x800
            # H or S in CIGAR
            and (
                "H" in read.cigarstring
                or "S" in read.cigarstring
            )
        ):
            self.diff[read.reference_start] += 1
            self.diff[read.reference_end] -= 1

    def close(self):
        current = 0
        for i in range(self.genome_length):
            current += self.diff[i]
            self.next.process(current)
