from conftest import UNI12, UNI13_RC, exact_alignment, random_seq
from flucut_posthoc.amplicons import build_units, classify_unmapped
from flucut_posthoc.seqio import revcomp

FLANK = random_seq(30, 20)


def units_of(reference, params, alignments):
    return build_units(alignments, {reference.name: reference}, params)


def test_complete_amplicon(reference, params):
    read = FLANK + reference.seq + revcomp(FLANK)
    aln = exact_alignment(reference, read, 30, 0, 600)
    (unit,), _, parts = units_of(reference, params, [aln])
    assert unit.category == "terminal+"
    assert (unit.left.kind, unit.left.origin, unit.left.pos) == \
        ("12", "terminal", 12)
    assert (unit.right.kind, unit.right.origin, unit.right.pos) == \
        ("13", "terminal", 586)
    assert unit.template_len == 575 and parts == [unit]


def test_same_amplicon_read_on_the_minus_strand(reference, params):
    read = FLANK + reference.seq + revcomp(FLANK)
    aln = exact_alignment(reference, read, 30, 0, 600, reverse=True)
    (unit,), _, _ = units_of(reference, params, [aln])
    assert unit.category == "terminal-"
    assert unit.read_span == (43, 618)


def test_internally_primed_amplicon(reference, params):
    # the primer replaced the 12 bases before position 300
    read = FLANK + UNI12 + reference.seq[300:] + revcomp(FLANK)
    aln = exact_alignment(reference, read, 42, 300, 300)
    (unit,), _, _ = units_of(reference, params, [aln])
    assert unit.category == "internal+"
    assert (unit.left.origin, unit.left.pos) == ("internal", 300)


def test_fragment_without_primers(reference, params):
    read = reference.seq[200:400]
    aln = exact_alignment(reference, read, 0, 200, 200)
    (unit,), _, _ = units_of(reference, params, [aln])
    assert unit.category == "noprimer_fragment" and unit.n_primers == 0


def test_internal_deletion_from_two_alignments(reference, params):
    read = FLANK + reference.seq[:150] + reference.seq[450:] + FLANK
    first = exact_alignment(reference, read, 30, 0, 150)
    second = exact_alignment(reference, read, 180, 450, 150)
    (unit,), _, _ = units_of(reference, params, [second, first])
    assert unit.category == "deletion"
    assert unit.deletions == [(150, 450)]
    assert unit.blocks == [(0, 150), (450, 600)]


def test_stretch_not_aligned_is_not_a_deletion(reference, params):
    # 120 bases of the read between the alignments, 120 of the reference
    read = reference.seq[:200] + random_seq(120, 21) + reference.seq[320:]
    first = exact_alignment(reference, read, 0, 0, 200)
    second = exact_alignment(reference, read, 320, 320, 280)
    (unit,), _, _ = units_of(reference, params, [first, second])
    assert unit.deletions == [] and unit.blocks == [(0, 600)]


def test_two_amplicons_joined_end_to_end(reference, params):
    amplicon = reference.seq
    read = FLANK + amplicon + amplicon + FLANK
    first = exact_alignment(reference, read, 30, 0, 600)
    second = exact_alignment(reference, read, 630, 0, 600)
    units, junctions, parts = units_of(reference, params, [first, second])
    assert len(units) == 2 and len(parts) == 2
    assert junctions == ["amplicons joined end to end"]


def test_amplicon_joined_to_one_which_does_not_map(reference, params):
    other = UNI12 + random_seq(200, 22) + UNI13_RC
    read = FLANK + reference.seq + other
    aln = exact_alignment(reference, read, 30, 0, 600)
    units, _, parts = units_of(reference, params, [aln])
    assert len(units) == 1 and parts == [units[0], None]


def test_reads_which_do_not_map(params):
    template = random_seq(150, 23)
    assert classify_unmapped(UNI12 + UNI13_RC, params).kind == "dimer"
    found = classify_unmapped(FLANK + UNI12 + template + UNI13_RC, params)
    assert found.kind == "offtarget" and found.template == template
    assert classify_unmapped(UNI12 + template, params).kind == "one_primer"
    assert classify_unmapped(template, params).kind == "no_primer"


def test_flanks_of_a_deletion_which_overlap_in_the_read(reference, params):
    # 40 bases of the read fit both flanks
    read = reference.seq[:200] + reference.seq[400:]
    first = exact_alignment(reference, read, 0, 0, 200)
    second = exact_alignment(reference, read, 200, 400, 200)
    second.qs, second.rs = 160, 360          # as if aligned 40 bases back
    second.bq, second.br, second.bl = [160], [360], [240]
    (unit,), _, _ = units_of(reference, params, [first, second])
    assert unit.category == "deletion"
    # the 40 bases are counted once, with the second flank
    assert unit.blocks == [(0, 160), (360, 600)]
    assert unit.deletions == [(160, 360)]


def test_deletion_with_a_flank_too_short_to_map(reference, params):
    # the primer, 20 bases of the start of the segment, then its end
    read = FLANK + reference.seq[:32] + reference.seq[400:] + FLANK
    aln = exact_alignment(reference, read, 62, 400, 200)
    (unit,), _, _ = units_of(reference, params, [aln])
    assert unit.category == "deletion"
    assert (unit.left.kind, unit.left.origin, unit.left.pos) == \
        ("12", "terminal", 12)
    assert unit.blocks == [(12, 32), (400, 600)]
    assert unit.deletions == [(32, 400)]
    assert unit.template_len == 207


def test_same_at_the_end_of_the_segment(reference, params):
    read = FLANK + reference.seq[:200] + reference.seq[568:] + FLANK
    aln = exact_alignment(reference, read, 30, 0, 200)
    (unit,), _, _ = units_of(reference, params, [aln])
    assert unit.category == "deletion"
    assert (unit.right.kind, unit.right.origin, unit.right.pos) == \
        ("13", "terminal", 586)
    assert unit.blocks == [(0, 200), (568, 587)]
    assert unit.deletions == [(200, 568)]


def test_short_flank_at_a_reference_without_its_terminus(params):
    from flucut_posthoc.references import Reference
    body = random_seq(600, 31)
    ref = Reference("ref", "segment 4", body[4:], 4, 0, 0)  # lacks 4 bases
    read = FLANK + UNI12 + body[:20] + body[404:] + FLANK
    aln = exact_alignment(ref, read, 62, 400, 196)
    (unit,), _, _ = build_units([aln], {"ref": ref}, params)
    assert unit.category == "deletion" and unit.left.kind == "12"
    assert unit.blocks == [(0, 16), (400, 596)]
