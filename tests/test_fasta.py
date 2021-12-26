import miniFasta as mf
from os import path, remove


def test_read():
    file_path = path.join(path.dirname(__file__), "test_data/test0.fasta")

    assert mf.read(file_path) == [
        mf.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        mf.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_read():

    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [mf.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
          mf.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")]

    mf.write(fo, file_path)

    fo_read = mf.read(file_path)
    remove(file_path)
    assert fo == fo_read


def test_read_strange_header():
    file_path = path.join(path.dirname(__file__), "test_data/test1.fasta")
    fos_heads = [fo.head for fo in mf.read(file_path)]
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
    assert mf.read(file_path)[0].body == b


def test_read_empty_lines():
    file_path = path.join(path.dirname(__file__), "test_data/test3.fasta")
    b0 = "IVPRFRKIDPRFSIDMMRLVGSFLKDREAEIIDGYGAQRSLNSVESADDTTRHFPSTVGVALGNAAIKELHRDENFDNCL"
    b1 = "MLPTRLANLGLESADHYPNPIQLNADDWDIPFEFELTHQVPTSVAVQYGSLSRAAATLERFVGS"
    fos = mf.read(file_path)

    assert fos[0].body == b0
    assert fos[1].body == b1


def test_translate_seq():
    assert mf.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert mf.translate_seq("HELLO") == "~"


def test_translate_toAmino():
    fo = mf.fasta_object("test", "CGGCCTTCTATCTTCTTC")
    fo.toAmino()
    assert fo.body == "RPSIFF"


def test_reverse_comp_seq():
    assert mf.reverse_comp("CGGCCTTCTATCTTCTTC") == "GAAGAAGATAGAAGGCCG"


def test_reverse_comp_toRevComp():
    fo = mf.fasta_object("test", "CGGCCTTCTATCTTCTTC")
    fo.toRevComp()
    assert fo.body == "GAAGAAGATAGAAGGCCG"


def test_print_fasta_by_func(capsys):
    file_path = path.join(path.dirname(__file__), "test_data/test0.fasta")
    mf.print_fasta(mf.read(file_path))

    # check if it prints the sequences correctly
    c = capsys.readouterr()
    assert c.out == ">Atlantic dolphin\nCGGCCTTCTATCTTCTTC\n>Pacific dolphin\nCTTTCTATCTCTTTCCTCT\n"


def test_len_fasta():
    assert len(mf.fasta_object("test", "abc")) == 3


def test_str_fasta():
    assert str(mf.fasta_object("test", "abc")) == ">test\nabc"


def test_eq_fasta():
    foa = mf.fasta_object("test", "abc")
    fob = mf.fasta_object("different header", "abc")
    foc = mf.fasta_object("different body", "zzz")

    assert foa == foa
    assert fob == foa
    assert foa == fob
    assert not foa == foc
    assert not foc == foa
