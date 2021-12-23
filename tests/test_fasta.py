from miniFasta import miniFasta as fasta
from os import path, remove


def test_read_fasta():
    file_path = path.join(path.dirname(__file__), "test_data/test0.fasta")

    assert fasta.read_fasta(file_path) == [
        fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_read_fasta():

    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [fasta.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
          fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")]

    fasta.write_fasta(fo, file_path)

    fo_read = fasta.read_fasta(file_path)
    remove(file_path)
    assert fo == fo_read


def test_read_strange_header():
    file_path = path.join(path.dirname(__file__), "test_data/test1.fasta")
    fos_heads = [fo.head for fo in fasta.read_fasta(file_path)]
    assert fos_heads == [""">FirstTestFASTA !"$%_-&/()=?'#'"",. :""", ">", "> Third >"]


def test_read_long_body():
    file_path = path.join(path.dirname(__file__), "test_data/test2.fasta")
    b = "".join(["CCAGTCTGGTCTCTCTCTAGATCAATTTTAACGGGCAAATTGTTGCTATTGCTCCAAATTCAATGGGACGCGCCTCT",
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
                 "GTCTGAACCACCCCTTCAAGAGTCGCTAGGAAGCTCTTGGTCAA"])
    assert fasta.read_fasta(file_path)[0].body == b


def test_read_empty_lines():
    file_path = path.join(path.dirname(__file__), "test_data/test3.fasta")
    b0 = "IVPRFRKIDPRFSIDMMRLVGSFLKDREAEIIDGYGAQRSLNSVESADDTTRHFPSTVGVALGNAAIKELHRDENFDNCL"
    b1 = "MLPTRLANLGLESADHYPNPIQLNADDWDIPFEFELTHQVPTSVAVQYGSLSRAAATLERFVGS"
    fos = fasta.read_fasta(file_path)

    assert fos[0].body == b0
    assert fos[1].body == b1


def test_translate_seq():
    assert fasta.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert fasta.translate_seq("HELLO") == "~"


def test_translate_toAmino():
    fo = fasta.fasta_object("test", "CGGCCTTCTATCTTCTTC")
    fo.toAmino()
    assert fo.body == "RPSIFF"


def test_reverse_comp_seq():
    assert fasta.reverse_comp("CGGCCTTCTATCTTCTTC") == "GAAGAAGATAGAAGGCCG"


def test_reverse_comp_toRevComp():
    fo = fasta.fasta_object("test", "CGGCCTTCTATCTTCTTC")
    fo.toRevComp()
    assert fo.body == "GAAGAAGATAGAAGGCCG"


def test_print_fasta_by_func(capsys):
    file_path = path.join(path.dirname(__file__), "test_data/test0.fasta")

    assert fasta.print_fasta(fasta.read_fasta(file_path)) is None

    # check if it prints the sequences correctly
    c = capsys.readouterr()
    assert c.out == ">Atlantic dolphin\nCGGCCTTCTATCTTCTTC\n>Pacific dolphin\nCTTTCTATCTCTTTCCTCT\n"


def test_len_fasta():
    assert len(fasta.fasta_object("test", "abc")) == 3


def test_str_fasta():
    assert str(fasta.fasta_object("test", "abc")) == ">test\nabc"


def test_eq_fasta():
    foa = fasta.fasta_object("test", "abc")
    fob = fasta.fasta_object("different header", "abc")
    foc = fasta.fasta_object("different body", "zzz")

    assert foa == foa
    assert fob == foa
    assert foa == fob
    assert not foa == foc
    assert not foc == foa
