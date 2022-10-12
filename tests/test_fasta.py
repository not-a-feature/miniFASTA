import miniFasta as mf

from os import path, remove
import pytest


def test_write_read():
    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [
        mf.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
        mf.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT"),
    ]

    mf.write(fo, file_path)

    fo_read = list(mf.read(file_path))
    remove(file_path)
    assert fo == fo_read


def test_write_byself():
    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = mf.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")

    fo.write(file_path)

    fo_read = list(mf.read(file_path))
    remove(file_path)
    assert [fo] == fo_read


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
    mf.print_fasta(list(mf.read(file_path)))

    # check if it prints the sequences correctly
    c = capsys.readouterr()
    assert c.out == ">Atlantic dolphin\nCGGCCTTCTATCTTCTTC\n>Pacific dolphin\nCTTTCTATCTCTTTCCTCT\n"


def test_len_fasta():
    assert len(mf.fasta_object(">test", "abc")) == 3


def test_str_fasta():
    assert str(mf.fasta_object("test", "abc")) == ">test\nabc"


def test_eq_fasta():
    foa = mf.fasta_object(">test", "abc")
    fob = mf.fasta_object(">different header", "abc")
    foc = mf.fasta_object(">different body", "zzz")

    assert foa == foa
    assert fob == foa
    assert foa == fob
    assert not foa == foc
    assert not foc == foa


def test_getter():
    fo = mf.fasta_object(">test", "abc")
    assert fo.body == fo.getSeq()
    assert fo.head == fo.getHead()


def test_iter_fasta():
    fo = mf.fasta_object(">test", "abc")
    body = ""
    for s in fo:
        body = body + s

    assert body == fo.body


def test_valid():
    fo = mf.fasta_object(">valid", "Ä'_**?.asdLLA")
    assert fo.valid()

    fo = mf.fasta_object(">valid", "ACGTUAGTGU", stype="NA")
    assert fo.valid()

    fo = mf.fasta_object(">invalid", "ACWYUOTGU", stype="NA")
    assert not fo.valid()

    fo = mf.fasta_object(">valid", "ACGTAGGT", stype="DNA")
    assert fo.valid()

    fo = mf.fasta_object(">invalid", "ACG TAGGTZ", stype="DNA")
    assert not fo.valid()

    fo = mf.fasta_object(">valid", "ACGUAGGU", stype="RNA")
    assert fo.valid()

    fo = mf.fasta_object(">invalid", "ACGTAGGT", stype="RNA")
    assert not fo.valid()

    fo = mf.fasta_object(">valid", "WYUOBJZ*X", stype="PROT")
    assert fo.valid()

    fo = mf.fasta_object(">invalid", "WY#X", stype="PROT")
    assert not fo.valid()

    fo = mf.fasta_object(">invalid", "ACGTAA", stype="ANY")
    assert fo.valid(allowedChars="ACGT")

    fo = mf.fasta_object(">invalid", "ACGTAA", stype="DNA")
    assert fo.valid(allowedChars="ACGT")

    fo = mf.fasta_object(">invalid", "FACG", stype="ANY")
    assert not fo.valid(allowedChars="ACGT")


@pytest.mark.xfail
def test_invalid_type_exception():
    fo = mf.fasta_object(">valid", "ACGTAGGT", stype="AA")
    assert fo.valid()
