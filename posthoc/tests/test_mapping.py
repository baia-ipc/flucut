from flucut_posthoc.mapping import parse_record, read_groups


def record(cigar: str, flag: int = 0, pos: int = 101) -> str:
    seq = "A" * 60
    return f"r1\t{flag}\tref\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t*\tNM:i:7\n"


def test_coordinates_of_an_alignment():
    _, aln, _ = parse_record(record("10S20M5I15M10S"))
    assert (aln.qs, aln.qe, aln.rs, aln.re) == (10, 50, 100, 135)
    assert aln.mismatches == 7 and not aln.reverse
    assert aln.read_start == 10


def test_reference_position_of_a_read_position():
    _, aln, _ = parse_record(record("10S20M5I15M10S"))
    assert aln.ref_of(10, True) == 100
    assert aln.ref_of(29, True) == 119
    assert aln.ref_of(32, True) == 120     # inside the insertion
    assert aln.ref_of(32, False) == 119
    assert aln.ref_of(6, True) == 96       # before the alignment
    assert aln.ref_of(52, False) == 137    # after it


def test_pieces_at_a_long_deletion():
    _, aln, _ = parse_record(record("20M300D20M2D20M"))
    assert aln.pieces(100) == [(0, 20, 100, 120), (20, 60, 420, 462)]


def test_reverse_strand_and_unmapped():
    _, aln, _ = parse_record(record("5S55M", flag=16))
    assert aln.reverse and aln.read_start == 0
    assert parse_record(record("*", flag=4))[1] is None


def test_records_grouped_by_read():
    lines = ["@SQ\tSN:ref\n", "a\t0\n", "a\t2048\n", "b\t4\n"]
    assert [(n, len(g)) for n, g in read_groups(lines)] == \
        [("a", 2), ("b", 1)]
