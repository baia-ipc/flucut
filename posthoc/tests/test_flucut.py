import json
from collections import Counter

import numpy as np
import pytest

from conftest import exact_alignment, random_seq
from flucut_posthoc.amplicons import build_units
from flucut_posthoc.agreement import origin_scores
from flucut_posthoc.flucut_summary import (load_group, on_references,
                                           read_coverage)
from flucut_posthoc.flucut_tags import (compare_read, parse_part,
                                        predicted_spans)

LENGTHS = {"ref": 600}


def test_tags_of_a_part():
    part = parse_part("575\tfs:Z:+\tfl:Z:uni12r\tol:Z:terminal"
                      "\tml:Z:ref:13:+\ttl:i:30\ttr:i:23")
    assert part["len"] == 575 and part["ml"] == ("ref", 12, False)
    assert part["tl"] == 30 and part["ol"] == "terminal"


def test_position_on_another_reference():
    lift = {"other": ("ref", np.arange(100) + 12)}
    part = parse_part("100\tfl:Z:uni12r\tml:Z:other:1:+", lift)
    assert part["ml"] == ("ref", 12, False)


def test_predicted_span_of_a_terminal_read():
    part = parse_part("575\tfs:Z:+\tfl:Z:u\tol:Z:terminal\tml:Z:ref:13:+")
    assert predicted_spans(part, LENGTHS) == ("terminal+", "ref",
                                              [(12, 587)])


def test_predicted_span_of_a_deletion():
    part = parse_part("300\tfs:Z:-\tfl:Z:u\tol:Z:terminal\tml:Z:ref:1:+"
                      "\tfr:Z:v\tor:Z:terminal\tmr:Z:ref:600:-")
    assert predicted_spans(part, LENGTHS) == (
        "deletion", "ref", [(0, 150), (450, 600)])


def test_read_without_located_end():
    part = parse_part("300\tfs:Z:+\tfl:Z:u\tol:Z:unknown")
    assert predicted_spans(part, LENGTHS)[0] == "not_located"


def test_category_tag_is_used_when_present():
    # no end located: the category is known only from the tag
    part = parse_part("300\tfs:Z:-\tfl:Z:u\tol:Z:unknown\tct:Z:internal")
    assert predicted_spans(part, LENGTHS) == ("internal-", None, [])
    part = parse_part("300\tfs:Z:.\tct:Z:unknown")
    assert predicted_spans(part, LENGTHS) == ("not_located", None, [])
    # an end located: the span from its position, the category from the tag
    part = parse_part("575\tfs:Z:+\tfl:Z:u\tol:Z:terminal\tml:Z:ref:13:+"
                      "\tct:Z:full")
    assert predicted_spans(part, LENGTHS) == ("terminal+", "ref",
                                              [(12, 587)])
    part = parse_part("300\tfs:Z:+\tml:Z:ref:13:+\tct:Z:fragment")
    assert predicted_spans(part, LENGTHS) == ("no_primer", "ref",
                                              [(12, 312)])
    part = parse_part("300\tfs:Z:-\tfl:Z:u\tol:Z:terminal\tml:Z:ref:1:+"
                      "\tfr:Z:v\tor:Z:terminal\tmr:Z:ref:600:-"
                      "\tct:Z:deletion")
    assert predicted_spans(part, LENGTHS) == (
        "deletion", "ref", [(0, 150), (450, 600)])


def test_agreement_of_a_read(reference, params):
    flank = random_seq(30, 30)
    read = flank + reference.seq + flank
    aln = exact_alignment(reference, read, 30, 0, 600)
    units, _, layout = build_units([aln], {"ref": reference}, params)
    counts: Counter = Counter()
    tags = ("599\tfs:Z:+\tfl:Z:uni12r\tol:Z:terminal\tml:Z:ref:13:+"
            "\tfr:Z:mbtuni-13\tor:Z:internal\tmr:Z:ref:587:-"
            "\ttl:i:30\ttr:i:43")
    predicted = compare_read([tags], layout, None, LENGTHS, {}, counts)
    assert predicted == [("ref", "terminal+", 12, 587)]
    assert counts["origin", "left", "terminal", "terminal"] == 1
    assert counts["origin", "right", "internal", "terminal"] == 1
    assert counts["position", "left", "0"] == 1
    assert counts["trim", "terminal", "terminus kept"] == 1
    assert counts["trim", "terminal", "whole primer removed"] == 1


def test_coverage_of_the_summary_of_flucut(tmp_path, reference):
    # two classes at 4 positions of a reference of 600 bases
    group = {"name": "sample", "processed": True, "coverage": [{
        "reference": "ref", "length": 600, "positions": [0, 200, 400, 599],
        "depth": {"terminal_plus": [10, 10, 20, 20],
                  "internal_plus": [0, 0, 10, 10],
                  "not_located": [1, 1, 1, 1]}}]}
    path = tmp_path / "flucut_summary.json"
    path.write_text(json.dumps({"min_depth": 50, "groups": [
        group, {"name": "other", "processed": True, "coverage": []}]}))
    section = load_group(str(path), None, "run1_sample")
    assert section["name"] == "sample" and section["min_depth"] == 50
    with pytest.raises(ValueError):
        load_group(str(path), "missing", "sample")
    coverage = read_coverage(section)
    assert list(coverage["ref"]["terminal+"]) == [10, 10, 20, 20]
    assert list(coverage["ref"]["not_located"]) == [1, 1, 1, 1]
    profiles = on_references(coverage, [reference], {})
    assert len(profiles["ref"]["terminal+"]) == 600
    assert profiles["ref"]["terminal+"][0] == 10
    assert profiles["ref"]["terminal+"][300] == 15
    assert profiles["ref"]["internal+"][599] == 10
    assert read_coverage({"coverage": None}) == {}


def test_coverage_on_another_reference(reference):
    # the reference of flucut has 12 bases more at its start
    coverage = {"other": {"positions": np.array([0., 12., 611.]),
                          "terminal+": np.array([4., 4., 4.])}}
    lift = {"other": ("ref", np.arange(612) - 12)}
    profiles = on_references(coverage, [reference], lift)
    assert profiles["ref"]["terminal+"][0] == 4
    assert profiles["ref"]["terminal+"][599] == 4


def test_coverage_of_several_references_is_added(reference):
    # flucut located reads in the chosen reference and in another one
    # of the segment, and names a third without any depth
    coverage = {"empty": {"positions": np.array([0., 599.])},
                "other": {"positions": np.array([0., 12., 611.]),
                          "terminal+": np.array([4., 4., 4.])},
                "ref": {"positions": np.array([0., 599.]),
                        "terminal+": np.array([10., 10.]),
                        "deletion": np.array([1., 1.])}}
    lift = {"empty": ("ref", np.arange(600)),
            "other": ("ref", np.arange(612) - 12)}
    profiles = on_references(coverage, [reference], lift)
    assert profiles["ref"]["terminal+"][300] == 14
    assert profiles["ref"]["deletion"][300] == 1


def test_scores_count_missed_and_false_primers():
    origin = {"left": {
        "internal": {"internal": 8, "terminal": 1, "no primer": 1},
        "terminal": {"terminal": 90, "internal": 2, "no primer": 10},
        "unknown": {"internal": 3, "terminal": 4},
        "no primer": {"internal": 7, "terminal": 6}}}
    scores = origin_scores(origin)["left"]
    assert scores["called"] == 101
    assert scores["agree"] == 98
    assert scores["internal_precision"] == 8 / 10
    assert scores["internal_recall"] == 8 / 20
    assert scores["terminal_precision"] == 90 / 102
    assert scores["terminal_recall"] == 90 / 101
    assert scores["internal_called"] == 10
    assert scores["internal_measured"] == 20
    assert origin_scores({"right": {}})["right"]["internal_recall"] is None
