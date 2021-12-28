![miniFASTA](https://github.com/not-a-feature/miniFASTA/raw/main/miniFASTA.png)

An easy FASTA object handler, reader, writer and translator for small to medium size projects without dependencies.

![Test Badge](https://github.com/not-a-feature/miniFASTA/actions/workflows/tests.yml/badge.svg)
![Download Badge](https://img.shields.io/pypi/dm/miniFASTA.svg)
## Installation
Using pip  / pip3:
```bash
pip install miniFasta
```
Or by source:
```bash
git clone git@github.com:not-a-feature/miniFASTA.git
cd miniFASTA
pip install .
```

## How to use
miniFASTA offers easy to use functions for fasta handling.
The five main parts are:
- fasta_object()
    - toAmino()
    - roRevComp()
    - len() / str() / eq()
- read()
- write()
- translate_seq()
- reverse_comp()


### fasta_object()
The core component of miniFASTA is the ```fasta_object()```. This object represents an entry in a FASTA file and consists of a head and body.

```python 
import miniFasta as mf
fo = mf.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")
print(fo.head) # >Atlantic dolphin
print(fo.body) # CGGCCTTCTATCTTCTTC

### Following functions are defined on a fasta_object():

str(fo) # will return:
# >Atlantic dolphin
# CGGCCTTCTATCTTCTTC

# Body length
len(fo) # will return 18, the length of the body

# Equality 
print(fo == fo) # True

fo_b = mf.fasta_object(">Same Body", "CGGCCTTCTATCTTCTTC")
print(fo == fo_b) # True

fo_c = mf.fasta_object(">Different Body", "ZZZZAGCTAG")
print(fo == fo_c) # False
```

**fasta_object(...).toAmino(translation_dict)**

Translates the body to an amino-acid sequence. See `tranlate_seq()` for more details.
```python 
fo.toAmino() 
print(fo.body) # Will return RPSIFF
d = {"CCG": "Z", "CTT": "A" ...}
fo.toAmino(d) 
print(fo.body) # Will return ZA...
```
**fasta_object(...).toRevComp(complement_dict)**

Converts the body to its reverse comlement. See `reverse_comp()` for more details.
```python 
fo.toRevComp() 
print(fo.body) # Will return GAAGAAGATAGAAGGCCG
```
## Reading FASTA files
`read()` is a fasta reader which is able to handle compressed and non-compressed files.
Following compressions are supported: zip, tar, tar.gz, gz. If multiple files are stored inside an archive, all files are read. 
This function returns a list of fasta_objects. 
The entries are usually casted to upper case letters. Set `read("path.fasta", upper=False)` to disable casting.

```python
fos = mf.read("dolphin.fasta") # List of fasta entries.
fos = mf.read("mouse.fasta", upper=False) # The entries won't be casted to upper case.
fos = mf.read("reads.tar.gz") # Is able to handle compressed files.
```

## Writing FASTA files
`write()` is a basic fasta reader.
It takes a single or a list of fasta_objects and writes it to the given path. 

The file is usually overwritten. Set `write(fo, "path.fasta", mode="a")` to append file.

```python
fos = mf.read("dolphin.fasta") # List of fasta entries
mf.write(fos, "new.fasta")
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
Copyright (C) 2021 by Jules Kreuer - @not_a_feature
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