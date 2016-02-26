# -*- coding: UTF-8 -*-


def read_fasta_info(fasta_file_path):
    """ Return id and length for each sequence read from a fasta file.

    Parameters:
        fasta_file_path (string): the fasta file

    Returns:
        [(int, int)]: a list of tuples (seq_id, seq_length) with the id and
            length of each sequence read from the fasta file.
    """
    fasta_info = []
    seq_id = None
    seq_length = 0

    with open(fasta_file_path, "r") as input_file:
        for line in input_file:
            # Check if the line is a comment
            if line[0] == ">":
                if seq_id is not None:
                    # Store the info of the previous sequence
                    fasta_info.append((seq_id, seq_length))
                # A new sequence starts
                seq_id = line[1:].split(" ", 1)[0]
                seq_length = 0
            else:
                # Add the lenghts, because the sequence is splitted over
                # several lines
                seq_length += len(line.strip())

    # Check if there was at least a sequence in the fasta file
    if seq_id is not None:
        # Store the info of the last sequence
        fasta_info.append((seq_id, seq_length))

    return fasta_info
