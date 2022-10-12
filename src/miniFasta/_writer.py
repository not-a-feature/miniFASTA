"""
miniFASTA: A simple toolbox for fasta files.

This is the writer part.

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""

from ._miniFasta import fasta_object as superFO


class fasta_object(superFO):
    def write(self, file_path: str, mode="w"):
        """
        Writes this fasta_object to a file.
        """
        write(self, file_path, mode)


def write(fasta_pairs, file_path: str, mode="w") -> None:
    """
    Writes a list of fasta_objects or a single one to a file.
    Takes fasta_objects as input.
    """

    if not isinstance(fasta_pairs, list):
        fasta_pairs = [fasta_pairs]

    with open(file_path, mode) as f:
        for fo in fasta_pairs:
            f.write(f"{fo.head}\n")
            body_len = len(fo.body)
            # Write only 70 chars per line
            for i in range(0, body_len, 70):
                f.write(f"{fo.body[i:i+70]}\n")
    return None
