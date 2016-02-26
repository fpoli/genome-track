# -*- coding: UTF-8 -*-
import math
from lib.log import *


class WigFileWriter:
    def __init__(self, wigfile, name=None, description=None, color=None):
        self.wigfile = wigfile
        self._write_track(name, description, color)
        self._write_fixed_step_header()

    def _write_track(self, name=None, description=None, color=None):
        track = "track type=wiggle_0"
        if name is not None:
            track += " name=\"{}\"".format(name)
        if description is not None:
            track += " description=\"{}\"".format(description)
        if color is not None:
            track += " color=\"{c[0]},{c[1]},{c[2]}\"".format(c=color)
        self.wigfile.write(track + "\n")

    def _write_fixed_step_header(self, chrom="genome", start=1, step=1, span=1):
        self.wigfile.write(
            "fixedStep chrom={} start={} step={} span={}\n".format(
                chrom, start, step, span
            )
        )

    def process(self, value):
        self.wigfile.write("{}\n".format(value))

    def close(self):
        self.wigfile.close()


class InsertLengthFilter:
    def __init__(self, next, min_len=None, max_len=None):
        self.next = next
        self.min = min_len
        self.max = max_len

    def process(self, read):
        if not (
            # Both reads are mapped
            read.is_paired                 # SAM flag 0x1
            and not read.is_unmapped       # SAM flag 0x4
            and not read.mate_is_unmapped  # SAM flag 0x8
            # Length is outside range
            and (
                (self.min is not None and abs(read.template_length) < self.min)
                or
                (self.max is not None and abs(read.template_length) > self.max)
            )
        ):
            self.next.process(read)

    def close(self):
        self.next.close()
