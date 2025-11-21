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
from pathlib import Path
from typing import Iterator, Union, List, IO, Any


def _maybe_byte_to_str(maybe_byte: Union[bytes, str]) -> str:
    """
    Convert bytes or string to string and strip trailing whitespace.

    Parameters
    ----------
        maybe_byte: Union[bytes, str]
            Input that might be bytes or string.

    Returns
    -------
        str
            Decoded and stripped string.
    """
    if isinstance(maybe_byte, bytes):
        return maybe_byte.decode("utf-8").rstrip()
    return str(maybe_byte).rstrip()


def _get_file_handlers(file_path: Path) -> List[IO[Any]]:
    """
    Get appropriate file handlers for compressed or uncompressed files.

    Parameters
    ----------
        file_path: Path
            Path to the file.

    Returns
    -------
        List[IO[Any]]
            List of file handlers.
    """
    handlers: List[IO[Any]] = []
    file_suffix = file_path.suffix.lower()

    # .zip file
    if file_suffix == ".zip":
        zip_handler = ZipFile(file_path, "r")
        for inner_file in zip_handler.namelist():
            handlers.append(zip_handler.open(inner_file, "r"))

    # .tar file
    elif file_suffix == ".tar":
        tar_handler = tarfile.open(file_path, "r")
        for inner_file in tar_handler.getmembers():
            extracted = tar_handler.extractfile(inner_file)
            if extracted is not None:
                handlers.append(extracted)

    # .gz file (including .tar.gz)
    elif file_suffix == ".gz":
        # Check if it's a .tar.gz file
        if file_path.stem.endswith(".tar"):
            tar_handler = tarfile.open(file_path, "r:gz")
            for inner_file in tar_handler.getmembers():
                extracted = tar_handler.extractfile(inner_file)
                if extracted is not None:
                    handlers.append(extracted)
        else:
            handlers.append(gzip.open(file_path, "rb"))

    # Regular uncompressed file
    else:
        handlers.append(open(file_path, "r"))

    return handlers


def read(
    file_path: str, upper: bool = True, seq: bool = False
) -> Union[Iterator[fasta_object], Iterator[str]]:
    """
    Read a compressed or non-compressed FASTA file and return an Iterator of fasta_objects.

    Supports: .fasta, .fa, .zip, .tar, .gz, .tar.gz file formats.

    Parameters
    ----------
        file_path: str
            Path to the FASTA file.
        upper: bool, default: True
            Convert sequences to uppercase letters.
        seq: bool, default: False
            Return only the sequences instead of fasta_object instances.

    Returns
    -------
        Union[Iterator[fasta_object], Iterator[str]]
            Iterator of fasta_object instances or sequence strings.

    Raises
    ------
        FileNotFoundError
            If the specified file does not exist.
    """
    path_obj = Path(file_path)

    if not path_obj.is_file():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    handlers = _get_file_handlers(path_obj)

    for handler in handlers:
        with handler:
            head = ""
            body: List[str] = []
            new_object = True

            # Cache method references to avoid repeated lookups
            body_append = body.append

            for maybe_byte_line in handler:
                # Convert byte string to string if needed
                line = _maybe_byte_to_str(maybe_byte_line)

                # Skip empty lines
                if not line:
                    continue

                # Handle header lines (starting with '>')
                if line[0] == ">":  # Faster than startswith for single char
                    # If this is not the first object, yield the previous one
                    if not new_object:
                        # Apply transformations once on final sequence
                        sequence = "".join(body)
                        if upper:
                            sequence = sequence.upper()

                        if seq:
                            yield sequence
                        else:
                            yield fasta_object(head, sequence)

                    # Start new object
                    head = line.strip()
                    body = []
                    body_append = body.append
                    new_object = False

                # Handle sequence lines
                else:
                    # Strip whitespace but defer uppercasing until final join
                    body_append(line.strip())

            # Yield the last element if any data was processed
            if not new_object and (head or body):
                # Apply transformations once on final sequence
                sequence = "".join(body)
                if upper:
                    sequence = sequence.upper()

                if seq:
                    yield sequence
                else:
                    yield fasta_object(head, sequence)
