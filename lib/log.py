# -*- coding: UTF-8 -*-
import sys


def printlog(*args):
    """Print a message to stderr.
    """
    print(*args, file=sys.stderr)
