import edlib

from conftest import UNI12, random_seq
from flucut_posthoc.primers import (OLIGOS, contains_primer, search_left,
                                    search_read_start)
from flucut_posthoc.seqio import revcomp


def test_primer_at_the_boundary():
    read = random_seq(40, 2) + UNI12 + random_seq(100, 3)
    hit = search_left(read, 52, OLIGOS, 10, 0.2)
    assert (hit.oligo.name, hit.start, hit.end, hit.errors) == \
        ("Uni12", 40, 52, 0)


def test_primer_inside_the_alignment():
    # the aligner aligned the primer bases too: the boundary is before it
    read = random_seq(40, 2) + UNI12 + random_seq(100, 3)
    hit = search_left(read, 40, OLIGOS, 10, 0.2)
    assert (hit.start, hit.end) == (40, 52)


def test_primer_with_errors_and_degenerate_base():
    primer = "AGCGAAAGCTGG"          # G variant, one substitution
    read = random_seq(30, 4) + primer + random_seq(100, 5)
    hit = search_left(read, 42, OLIGOS, 10, 0.2)
    assert hit.oligo.name == "Uni12" and hit.errors == 1


def test_primer_far_from_the_boundary_is_not_taken():
    read = UNI12 + random_seq(200, 6)
    assert search_left(read, 150, OLIGOS, 10, 0.2) is None


def test_read_beginning_inside_the_primer():
    read = UNI12[3:] + random_seq(100, 7)
    hit = search_left(read, 9, OLIGOS, 10, 0.2)
    assert hit.end == 9 and hit.missing == 3


def test_no_primer_in_a_random_read():
    read = random_seq(300, 8)
    assert search_left(read, 60, OLIGOS, 10, 0.2) is None


def test_primer_at_the_start_of_a_read_not_mapped():
    read = random_seq(25, 9) + "AGTAGAAACAAGG" + random_seq(80, 10)
    hit = search_read_start(read, OLIGOS, 100, 0.2)
    assert hit.oligo.name == "Uni13" and hit.end == 38


def test_contains_primer_on_either_strand():
    assert contains_primer(random_seq(50, 11) + revcomp(UNI12), OLIGOS)
    assert not contains_primer(random_seq(200, 12), OLIGOS)


def test_controls_are_far_from_the_primers():
    real = ["AGCAAAAGCAGG", "AGCGAAAGCAGG", "AGTAGAAACAAGG",
            "ACGCGTGATCAGTAGAAACAAGG", "CGGGTTATTAGTAGAAACAAGG"]
    real += [revcomp(s) for s in real]
    for oligo in OLIGOS:
        assert len(set(oligo.controls)) == len(oligo.controls)
        for control in oligo.controls:
            assert sorted(control) == sorted(oligo.seq)
            for seq in real:
                distance = edlib.align(
                    control, seq, mode="HW",
                    additionalEqualities=[("R", "A"), ("R", "G")]
                )["editDistance"]
                assert distance >= 5
