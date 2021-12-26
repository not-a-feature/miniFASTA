# miniFASTA
 An easy FASTA object handler, reader, writer and translator for small to medium size projects without dependencies.

![Test Badge](https://github.com/not-a-feature/miniFASTA/actions/workflows/tests.yml/badge.svg)
![Download Badge](https://img.shields.io/pypi/dm/miniFASTA.svg)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fnot-a-feature%2FminiFASTA.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2Fnot-a-feature%2FminiFASTA?ref=badge_shield)
## Installation
Using pip  / pip3:
```bash
pip install miniFasta
```
Or by source:
```bash
git clone git@github.com:not-a-feature/miniFASTA.git
cd miniFASTA
pip3 install .
```

## How to use
miniFASTA offers easy to use functions for fasta handling.
The five main parts are:
- fasta_object()
- read_fasta()
- write_fasta()
- translate_seq()
- reverse_comp()


### fasta_object()
The core component of miniFASTA is the ```fasta_object()```. This object represents an entry in a FASTA file and consists of a head and body.

```python 
from miniFasta import miniFasta as fasta
fo = fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")
print(fo.head) # >Atlantic dolphin
print(fo.body) # CGGCCTTCTATCTTCTTC
```

### Following functions are defined on a fasta_object():

**str()**
```python 
str(fo) # will return:
# >Atlantic dolphin
# CGGCCTTCTATCTTCTTC
```
**len()**
```python 
len(fo) # will return 18, the length of the body
```
**==**
```python 
print(fo == fo) # True

fo_b = fasta.fasta_object(">Same Body", "CGGCCTTCTATCTTCTTC")
print(fo == fo_b) # True

fo_c = fasta.fasta_object(">Different Body", "ZZZZAGCTAG")
print(fo == fo_c) # False
```
**toAmino(translation_dict)**

Translates the body to an amino-acid sequence. See `tranlate_seq()` for more details.
```python 
fo.toAmino() 
print(fo.body) # Will return RPSIFF
d = {"CCG": "Z", "CTT": "A" ...}
fo.toAmino(d) 
print(fo.body) # Will return ZA...
```
**toRevComp(complement_dict)**

Converts the body to its reverse comlement. See `reverse_comp()` for more details.
```python 
fo.toRevComp() 
print(fo.body) # Will return GAAGAAGATAGAAGGCCG
```
## Reading FASTA files
`read_fasta()` is a basic fasta reader.
It reads a fasta-style file and returns a list of fasta_objects.
The entries are usually casted to upper case letters. Set `read_fasta("path.fasta", upper=False)` to disable casting.

```python
fos = fasta.read_fasta("dolphin.fasta") # List of fasta entries
fos = fasta.read_fasta("cat.fasta", upper=False)
```

## Writing FASTA files
`write_fasta()` is a basic fasta reader.
It takes a single or a list of fasta_objects and writes it to the given path. 

The file is usually overwritten. Set `write_fasta(fo, "path.fasta", mode="a")` to append file.

```python
fos = fasta.read_fasta("dolphin.fasta") # List of fasta entries
fasta.write_fasta(fos, "new.fasta")
```
## Sequence translation
`translate_seq()` translates a sequence starting at position 0.
Unless translation_dict is provided, the standart bacterial code is used. If the codon was not found, it will be replaced by an `~`. Tailing bases that do not fit into a codon will be ignored.

```python 
fasta.translate_seq("CGGCCTTCTATCTTCTTC") # Will return RPSIFF

d = {"CGG": "Z", "CTT": "A"}
fasta.translate_seq("CGGCTT", d) # Will return ZA.
```

## Reverse Complement
`reverse_comp()` converts a sequence to its reverse comlement.
Unless complement_dict is provided, the standart complement is used. If no complement was found, the nucleotide remains unchanged.
```python 
fasta.reverse_comp("CGGCCTTCTATCTTCTTC") # Will return GAAGAAGATAGAAGGCCG

d = {"C": "Z", "T": "Y"}
fasta.reverse_comp("TC", d) # Will return ZY
```

## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fnot-a-feature%2FminiFASTA.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2Fnot-a-feature%2FminiFASTA?ref=badge_large)