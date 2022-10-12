"""
miniFASTA: A simple toolbox for fasta files

See: https://github.com/not-a-feature/miniFASTA
Or:  https://pypi.org/project/miniFasta/

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""
from typing import Iterator
from dataclasses import dataclass

# Usual translation dictionary according to
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
translation_dict = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Usual reverse complement
complement_dict = {
    "A": "T",
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
    "H": "D",
}


@dataclass
class fasta_object:
    head: str
    body: str
    stype: str

    def __init__(self, head: str, body: str, stype: str = "any"):
        """
        Object to keep a fasta entry.
        Input:
            head:  str, head of fasta entry.
            body:  str, body of fasta entry.
            stype: str, type of the sequence.
                        Check with self.valid() if the body consits of valid characters.
                        One of  NA: Allows all Nucleic Acid Codes (DNA & RNA)
                               DNA: Allows all IUPAC DNA Codes
                               RNA: Allows all IUPAC DNA Codes
                              PROT: Allows all IUPAC Aminoacid Codes
                               ANY: Allows all characters [default]
        """
        if head.startswith(">"):
            self.head = head
        else:
            self.head = f">{head}"

        self.body = body

        if stype.upper() in ["NA", "DNA", "RNA", "PROT", "ANY"]:
            self.stype = stype.upper()
        else:
            raise RuntimeError("fasta object type must be one of 'dna', 'prot' or 'any'.")

    def __str__(self) -> str:
        """
        Magic method to allow fasta_object printing.
        """
        out = f"{self.head}\n"

        # Print only 70 chars per line
        for i in range(0, len(self.body), 70):
            out += f"{self.body[i:i+70]}\n"

        # Remove tailing newline
        return out[:-1]

    def __eq__(self, o) -> bool:
        """
        Magic method to allow equality check on fasta_objects.
        Does not check for header equality.
        """
        return self.body == o.body  # type:ignore

    def __len__(self) -> int:
        """
        Magic method to allow len() on fasta_objects.
        Returns the length of the body
        """
        return len(self.body)

    def __iter__(self) -> Iterator[str]:
        """
        Magic method to iterate through the sequence.
        """
        return iter(self.body)

    def getHead(self) -> str:
        """
        Getter method to return the head / sequence id.
        """
        return self.head

    def getSeq(self) -> str:
        """
        Getter method to return the sequence.
        """
        return self.body

    def valid(self, allowedChars: str = "") -> bool:
        """
        Checks if this fasta_object is valid.
        stype of fasta_object needs to be set in order to check for illegal characters in its body.
        Input:
            allowedChars: str, optional to overwrite default settings.
        """
        if 250000 <= len(self.body):
            return False

        if not allowedChars:
            if self.stype == "ANY":
                return True
            elif self.stype == "PROT":
                allowedChars = "ACDEFGHIKLMNPWRSTVWYUOBJZ*X-."
            else:
                allowedChars = "ACGRYSWKMBDHVN"
                if self.stype == "NA":
                    allowedChars += "TU"
                elif self.stype == "DNA":
                    allowedChars += "T"
                elif self.stype == "RNA":
                    allowedChars += "U"

        return all(c in allowedChars for c in self.body)

    def toAmino(self, d=translation_dict) -> None:
        """
        Translates the dna sequence of a fasta_object to amino-acids.
        Reading frame starts at position 0, tailing bases will be ignored.
        Attention: Will replace triplet with ~ if not found.
        Input:
            d: dict, dictionary of translation.
        """
        self.body = translate_seq(self.body, d)

    def toRevComp(self, d=complement_dict) -> None:
        """
        Reverses complement of sequence.
        If no complement was found, the nucleotide remains unchanged.
        Input:
            d: dict, dictionary of complement.
        """
        self.body = reverse_comp(self.body, d)


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
        codon = seq[i : i + 3]
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
