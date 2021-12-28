"""
miniFASTA: A simple toolbox for fasta files

See: https://github.com/not-a-feature/miniFASTA
Or:  https://pypi.org/project/miniFasta/

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""

from os import path
from typing import List
from zipfile import ZipFile
import gzip
import tarfile

# Usual translation dictionary according to
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
translation_dict = {"TTT": "F", "TTC": "F",
                    "TTA": "L", "TTG": "L",
                    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                    "TAT": "Y", "TAC": "Y",
                    "TAA": "*", "TAG": "*", "TGA": "*",
                    "TGT": "C", "TGC": "C",
                    "TGG": "W",
                    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                    "CAT": "H", "CAC": "H",
                    "CAA": "Q", "CAG": "Q",
                    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                    "ATT": "I", "ATC": "I", "ATA": "I",
                    "ATG": "M",
                    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                    "AAT": "N", "AAC": "N",
                    "AAA": "K", "AAG": "K",
                    "AGT": "S", "AGC": "S",
                    "AGA": "R", "AGG": "R",
                    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                    "GAT": "D", "GAC": "D",
                    "GAA": "E", "GAG": "E",
                    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

# Usual reverse complement
complement_dict = {"A": "T",
                   "G": "C",
                   "C": "G",
                   "T": "A",
                   "U": "A",
                   "R": "Y",
                   "Y": "R",
                   "S": "S",
                   "W": "W",
                   "K": "M",
                   "M": "K",
                   "B": "V",
                   "V": "B",
                   "D": "H",
                   "H": "D"}


class fasta_object():
    def __init__(self, head: str, body: str):
        """
        Object to keep a valid fasta object
        """
        if head.startswith(">"):
            self.head = head
        else:
            self.head = f">{head}"

        self.body = body

    def __str__(self):
        """
        Magic method to allow fasta_object printing.
        """
        out = f'{self.head}\n'

        # Print only 70 chars per line
        for i in range(0, len(self.body), 70):
            out += f'{self.body[i:i+70]}\n'

        # Remove tailing newline
        return out[:-1]

    def __repr__(self):
        """
        Magic method to allow printing of fasta_object representation.
        """
        return f'fasta_object("{self.head}", "{self.body}")'

    def __eq__(self, o):
        """
        Magic method to allow equality check on fasta_objects.
        Does not check for header equality.
        """
        return self.body == o.body

    def __len__(self):
        """
        Magic method to allow len() on fasta_objects.
        Does not check for header-length equality.
        """
        return len(self.body)

    def toAmino(self, d=translation_dict):
        """
        Translates the dna sequence of a fasta_object to amino-acids.
        Reading frame starts at position 0, tailing bases will be ignored.
        Attention: Will replace triplet with ~ if not found.
        Input:
            d: dict, dictionary of translation.
        """
        self.body = translate_seq(self.body, d)

    def toRevComp(self,  d=complement_dict):
        """
        Reverses complement of sequence.
        If no complement was found, the nucleotide remains unchanged.
        Input:
            d: dict, dictionary of complement.
        """
        self.body = reverse_comp(self.body, d)


def read(file_path: str, upper: bool = True) -> List[fasta_object]:
    """
    Reads a compressed or non-compressed fasta file and returns a list of fasta_objects.
    Zip, tar, gz, tar.gz files are supported.
    Attention: Encoding characters (backslash) will work under certain conditions.

    Input:
        file_path: str, path to folder / file
        upper: bool, cast sequence to upper-case letters.

    Returns:
        fasta_objects: list of fasta_object
    """

    if not path.isfile(file_path):
        raise FileNotFoundError("Fasta File not found!")

    handlers = []
    file_type = file_path.split(".")[-1]

    if file_type in ["zip", "tar", "gz"]:
        # .zip file
        if file_type == "zip":
            zipHandler = ZipFile(file_path, 'r')
            # Create handler for every file in zip
            for inner_file in zipHandler.namelist():  # type:ignore
                handlers.append(zipHandler.open(inner_file, "r"))  # type:ignore
        # .tar file
        elif file_type == "tar":
            tarHandler = tarfile.open(file_path, "r")
            # Create handler for every file in tar
            for inner_file in tarHandler.getmembers():  # type:ignore
                handlers.append(tarHandler.extractfile(inner_file))  # type:ignore
        # .gz file
        elif file_type == "gz":
            # tar.gz file
            if file_path.split(".")[-2] == "tar":
                # Create handler for every file in tar.gz
                tarHandler = tarfile.open(file_path, "r")
                for inner_file in tarHandler.getmembers():  # type:ignore
                    handlers.append(tarHandler.extractfile(inner_file))  # type:ignore
            else:
                # .gz file
                handlers = [gzip.open(file_path, "r")]  # type:ignore
    else:
        handlers = [open(file_path, "r")]  # type:ignore

    fasta_objects = []
    for h in handlers:
        with h:
            head = ""
            body = ""
            newObject = True

            for maybeByteLine in h:
                # Convert byte string to string
                if isinstance(maybeByteLine, bytes):
                    line = maybeByteLine.decode('utf-8')
                else:
                    line = maybeByteLine  # type:ignore

                # Go through each line
                # First Header
                if newObject and line.startswith(">"):
                    head = line.strip()
                    body = ""
                    newObject = False
                # N-th Header
                elif line.startswith(">"):
                    fasta_objects.append(fasta_object(head, body))
                    head = line.strip()
                    body = ""
                # Sequence
                else:
                    addBody = line.strip()
                    if upper:
                        addBody = addBody.upper()

                    body += addBody
            # append last element
            fasta_objects.append(fasta_object(head, body))
    return fasta_objects


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


def print_fasta(fasta) -> None:
    """
    Prints a single or a list of fasta_objects.
    """

    if not isinstance(fasta, list):
        fasta = [fasta]

    for fo in fasta:
        print(fo)
    return None


def __maybeFind(key, d, alt):
    """
    Tries to find key in dict but has a fallback.
    Input:
        key: hashable, key to find.
        d: dict, dictionary to search.
        alt: alternative if key not found.
    Returns:
        v: found value or alt
    """
    try:
        return d[key]
    except KeyError:
        return alt


def translate_seq(seq: str, d=translation_dict) -> str:
    """
    Translates a DNA sequence to a AA sequence.
    Reading frame starts at position 0, tailing bases will be ignored.
    Attention: Will replace triplet with ~ if not found.

    To translate a fasta_object use object.toAmino()

    Input:
        seq: String, sequence to translate.
        d: dict, dictionary of translation.
    Returns:
        translated: String, translated sequence.
    """

    translated = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if not len(codon) == 3:
            break
        translated += __maybeFind(codon, d, "~")
    return translated


def reverse_comp(seq: str, d=complement_dict) -> str:
    """
    Reverses complement of sequence.
    If no complement was found, the nucleotides remains unchanged.
    Input:
        seq: String, sequence to compute the reverse complement.
        d: dict, dictionary of complement.
    Returns:
        rev: String, translated sequence
    """
    if seq == "":
        return ""
    rev = reversed(seq)
    return "".join(map(lambda b: __maybeFind(b, complement_dict, b), rev))  # type:ignore
