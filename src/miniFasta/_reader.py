"""
miniFASTA: A simple toolbox for fasta files.

This is the reader part.

@author: Jules Kreuer / not_a_feature
License: GPL-3.0
"""

from ._miniFasta import fasta_object

from zipfile import ZipFile
import gzip
import tarfile

from os import path
from typing import List


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
