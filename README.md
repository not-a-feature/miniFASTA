![miniFASTA](https://github.com/not-a-feature/miniFASTA/raw/main/miniFASTA.png)

A simple FASTA read and write toolbox for small to medium size projects.


[![DOI](https://zenodo.org/badge/440126588.svg)](https://zenodo.org/badge/latestdoi/440126588)
![Test Badge](https://github.com/not-a-feature/miniFASTA/actions/workflows/tests.yml/badge.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)<br>
![Download Badge](https://img.shields.io/pypi/dm/miniFASTA.svg)
![Python Version Badge](https://img.shields.io/pypi/pyversions/miniFASTA)
[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](https://anaconda.org/conda-forge/minifasta)



FASTA files are text-based files for storing nucleotide or amino acid sequences.
Reading such files is not particularly difficult, yet most off the shelf packages are overloaded with strange dependencies.

miniFASTA offers an alternative to this and brings many useful functions without relying on third party packages.

## Installation
Using pip  / pip3:
```bash
pip install miniFasta
```
Or by source with pip:
```bash
git clone git@github.com:not-a-feature/miniFASTA.git
cd miniFASTA
pip install .
```
Or by conda:
```bash
conda install -c conda-forge minifasta
```

## How to use
miniFASTA offers easy to use functions for fasta handling.
The five main parts are:
- read()
- write()
- fasta_object()
    - toAmino()
    - roRevComp()
    - valid()
    - len() / str() / eq() / iter()
- translate_seq()
- reverse_comp()

## Reading FASTA files
`read()` is a fasta reader which is able to handle compressed and non-compressed files.
Following compressions are supported: zip, tar, tar.gz, gz. If multiple files are stored inside an archive, all files are read.
This function returns a Iterator of fasta_objects. If only the sequences should be returnes set the positional argument `seq=True`.
The entries are usually casted to upper case letters. Set `read("path.fasta", upper=False)` to disable casting.

```python
# Read fasta_objects
fos = mf.read("dolphin.fasta") # Iterator of fasta_objects.
fos = list(fos) # Casts the iterator to list of fasta_objects

# Read only the sequence
fasta_strings = mf.read("dolphin.fasta", seq=True) # Iterator of string.
fasta_strings = [fo.body for fo in mf.read("dolphin.fasta")] # Alternative

# Options and compressed files
fos = mf.read("mouse.fasta", upper=False) # The entries won't be casted to upper case.
fos = mf.read("reads.tar.gz") # Is able to handle compressed files.
```

## Writing FASTA files
`write()` is a basic fasta writer.
It takes a single or a list of fasta_objects and writes it to the given path.

The file is usually overwritten. Set `write(fo, "path.fasta", mode="a")` to append file.

```python
fos = mf.read("dolphin.fasta") # Iterator of fasta entries
fos = list(fos) # Materialize
mf.write(fos, "new.fasta")
```

### fasta_object()
The core component of miniFASTA is the ```fasta_object()```. This object represents an FASTA entry and consists of a head and body.

```python
import miniFasta as mf
fo = mf.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC", stype="DNA")
fo.getHead() or fo.head
# >Atlantic dolphin

fo.getSeq() or fo.body
# CGGCCTTCTATCTTCTTC

### Following functions are defined on a fasta_object():

str(fo) # will return:
# >Atlantic dolphin
# CGGCCTTCTATCTTCTTC

# Body length
len(fo) # will return 18, the length of the body

# Equality
fo == fo # True

fo_b = mf.fasta_object(">Same Body", "CGGCCTTCTATCTTCTTC")
fo == fo_b # True

fo_c = mf.fasta_object(">Different Body", "ZZZZAGCTAG")
fo == fo_c # False

for s in fo:
    # Iterates through the sequence of fo.
```

**fasta_object(...).valid()**

Checks if the body contains invalid characters.
_stype_ of fasta_object needs to be set in order to check for illegal characters in its body.

stype is one of:
- ANY : [default] Allows all characters.
- NA  : Allows all Nucleic Acid Codes (DNA & RNA).
- DNA : Allows all IUPAC DNA Codes.
- RNA : Allows all IUPAC RNA Codes.
- PROT: Allows all IUPAC Aminoacid Codes.

Optional: allowedChars can be set to overwrite default settings.

```python
# The default object allows all characters.
# True
fasta_object(">valid", "Ä'_**?.asdLLA").valid()

# Only if stype is specified, valid can check for illegal characters.
# True
fasta_object(">valid", "ACGTUAGTGU", stype="NA").valid()

# False, as W is not allowed for DNA/RNA
fasta_object(">invalid", "ACWYUOTGU", stype="NA").valid()

# True
fasta_object(">valid", "AGGATTA", stype="ANY").valid(allowedChars = "AGTC")

# True, as stype is ignored if allowedChars is set.
fasta_object(">valid", "WYU", stype="DNA").valid(allowedChars = "WYU")
```

**fasta_object(...).toAmino(translation_dict)**

Translates the body to an amino-acid sequence. See `tranlate_seq()` for more details.
```python
fo.toAmino()
fo.getBody() # Will return RPSIFF
d = {"CCG": "Z", "CTT": "A" ...}
fo.toAmino(d)
fo.getBody # Will return ZA...
```
**fasta_object(...).toRevComp(complement_dict)**

Converts the body to its reverse comlement. See `reverse_comp()` for more details.
```python
fo.toRevComp()
fo.getBody # Will return GAAGAAGATAGAAGGCCG
```

## Sequence translation
`translate_seq()` translates a sequence starting at position 0.
Unless translation_dict is provided, the standart bacterial code is used. If the codon was not found, it will be replaced by an `~`. Tailing bases that do not fit into a codon will be ignored.

```python
mf.translate_seq("CGGCCTTCTATCTTCTTC") # Will return RPSIFF

d = {"CGG": "Z", "CTT": "A"}
mf.translate_seq("CGGCTT", d) # Will return ZA.
```

## Reverse Complement
`reverse_comp()` converts a sequence to its reverse comlement.
Unless complement_dict is provided, the standart complement is used. If no complement was found, the nucleotide remains unchanged.
```python
mf.reverse_comp("CGGCCTTCTATCTTCTTC") # Will return GAAGAAGATAGAAGGCCG

d = {"C": "Z", "T": "Y"}
mf.reverse_comp("TC", d) # Will return ZY
```

## License
```
Copyright (C) 2022 by Jules Kreuer - @not_a_feature
This piece of software is published unter the GNU General Public License v3.0
TLDR:

| Permissions      | Conditions                   | Limitations |
| ---------------- | ---------------------------- | ----------- |
| ✓ Commercial use | Disclose source              | ✕ Liability |
| ✓ Distribution   | License and copyright notice | ✕ Warranty  |
| ✓ Modification   | Same license                 |             |
| ✓ Patent use     | State changes                |             |
| ✓ Private use    |                              |             |
```
Go to [LICENSE.md](https://github.com/not-a-feature/miniFASTA/blob/main/LICENSE) to see the full version.

## Dependencies
In addition to packages included in Python 3, this piece of software uses 3rd-party software packages for development purposes that are not required in the published version.
Go to [DEPENDENCIES.md](https://github.com/not-a-feature/miniFASTA/blob/main/DEPENDENCIES.md) to see all dependencies and licenses.
