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
from typing import Iterator, Union


def __maybeByteToStr(maybeByte) -> str:
    if isinstance(maybeByte, bytes):
        return maybeByte.decode("utf-8").rstrip()
    return str(maybeByte).rstrip()


def read(
    file_path: str, upper: bool = True, seq: bool = False
) -> Union[Iterator[fasta_object], Iterator[str]]:
    """
    Reads a compressed or non-compressed fasta file and returns a Iterator of fasta_objects.
    Zip, tar, gz, tar.gz files are supported.
    Attention: Encoding characters (backslash) will work under certain conditions.

    Input:
        file_path: str, path to folder / file.
        upper: bool, cast sequences to upper-case letters.
        seq: bool, return only the sequences.

    Returns:
        fasta_objects: Iterator of fasta_object or Iterator of strings.
    """

    if not path.isfile(file_path):
        raise FileNotFoundError("Fasta File not found!")

    handlers = []
    file_type = path.splitext(file_path)[1]

    # Compressed files
    # .zip file
    if file_type == ".zip":
        zipHandler = ZipFile(file_path, "r")
        # Create handler for every file in zip
        for inner_file in zipHandler.namelist():  # type:ignore
            handlers.append(zipHandler.open(inner_file, "r"))  # type:ignore
    # .tar file
    elif file_type == ".tar":
        tarHandler = tarfile.open(file_path, "r")
        # Create handler for every file in tar
        for inner_file in tarHandler.getmembers():  # type:ignore
            handlers.append(tarHandler.extractfile(inner_file))  # type:ignore
    # .gz file
    elif file_type == ".gz":
        # tar.gz file
        inner_file_type = path.splitext(path.splitext(file_path)[0])[1]
        if inner_file_type == ".tar":
            # Create handler for every file in tar.gz
            tarHandler = tarfile.open(file_path, "r")
            for inner_file in tarHandler.getmembers():  # type:ignore
                handlers.append(tarHandler.extractfile(inner_file))  # type:ignore
        else:
            # .gz file
            handlers = [gzip.open(file_path, "r")]  # type:ignore
    else:
        handlers = [open(file_path, "r")]  # type:ignore

    for h in handlers:
        with h:
            head = ""
            body = ""
            newObject = True

            for maybe_byte_line in h:
                # Convert byte string to string
                line = __maybeByteToStr(maybe_byte_line)

                # Go through each line

                # First Header
                if newObject and line.startswith(">"):
                    head = line.strip()
                    body = ""
                    newObject = False

                # N-th Header
                elif line.startswith(">"):
                    # Yield only sequence or complete fasta_object
                    if seq:
                        yield body
                    else:
                        yield fasta_object(head, body)
                    head = line.strip()
                    body = ""

                # Sequence
                else:
                    addBody = line.strip()
                    if upper:
                        addBody = addBody.upper()

                    body += addBody

            # Yield last element
            # Yield only sequence or complete fasta_object
            if seq:
                yield body
            else:
                yield fasta_object(head, body)
