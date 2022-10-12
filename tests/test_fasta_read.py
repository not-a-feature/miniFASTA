import miniFasta as mf

from os import path
import pytest

dolphin = [
    mf.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
    mf.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT"),
]

multi = dolphin.copy()
multi.extend([mf.fasta_object(">RANDOM", "ASGDTASGDTASGD"), mf.fasta_object(">R2", "ASDZASDJ")])


@pytest.mark.parametrize(
    "file_path, expected",
    [
        (path.join(path.dirname(__file__), "test_data/test0.fasta"), dolphin),
        (path.join(path.dirname(__file__), "test_data/test.fasta.zip"), dolphin),
        (path.join(path.dirname(__file__), "test_data/test.fasta.tar"), dolphin),
        (path.join(path.dirname(__file__), "test_data/test.fasta.tar.gz"), dolphin),
        (path.join(path.dirname(__file__), "test_data/test.fasta.gz"), dolphin),
        (path.join(path.dirname(__file__), "test_data/test.multi.zip"), multi),
        (path.join(path.dirname(__file__), "test_data/test.multi.tar"), multi),
        (path.join(path.dirname(__file__), "test_data/test.multi.tar.gz"), multi),
    ],
)
def test_read(file_path, expected):
    assert list(mf.read(file_path)) == expected


def test_read_seq_only():
    file_path = path.join(path.dirname(__file__), "test_data/test0.fasta")
    assert list(mf.read(file_path, seq=True)) == [d.body for d in dolphin]


def test_read_strange_header():
    file_path = path.join(path.dirname(__file__), "test_data/test1.fasta")
    fos_heads = [fo.head for fo in mf.read(file_path)]
    assert fos_heads == [""">FirstTestFASTA !"$%_-&/()=?'#'"",. :""", ">", "> Third >"]


def test_read_long_body():
    file_path = path.join(path.dirname(__file__), "test_data/test2.fasta")
    b = "".join(
        [
            "CCAGTCTGGTCTCTCTCTAGATCAATTTTAACGGGCAAATTGTTGCTATTGCTCCAAATTCAATGGGACGCGCCTCT",
            "TAGTCAGAAAGTACACGCAAACTACATGCTAGTAATGCCTCTGGGCAGGTATGTGGAGGACCGGCATTAGCGTGGCGC",
            "GTTCGCAGGGGGATTGGCGAACCGAAGGTGACCTTAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCTGAACC",
            "ACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCTGAAC",
            "CACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCTGAA",
            "CCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCTGA",
            "ACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCTG",
            "AACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTCT",
            "GAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGTC",
            "TGAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAGT",
            "CTGAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACAG",
            "TCTGAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAATAGGGATATTGTGTTGGGAGCGCGGTCCCCCTACA",
            "GTCTGAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAA",
        ]
    )
    assert next(mf.read(file_path)).body == b


def test_read_empty_lines():
    file_path = path.join(path.dirname(__file__), "test_data/test3.fasta")
    b0 = "IVPRFRKIDPRFSIDMMRLVGSFLKDREAEIIDGYGAQRSLNSVESADDTTRHFPSTVGVALGNAAIKELHRDENFDNCL"
    b1 = "MLPTRLANLGLESADHYPNPIQLNADDWDIPFEFELTHQVPTSVAVQYGSLSRAAATLERFVGS"
    fos = list(mf.read(file_path))

    assert fos[0].body == b0
    assert fos[1].body == b1
