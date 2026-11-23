"""The agreement checks and the pages written from a summary."""

import json

import numpy as np

from flucut_posthoc.agreement import (grade, origin_table, sample_checks,
                                      segment_rows, verdict)
from flucut_posthoc.origins import (MEASURED_ORIGIN, depth_stats,
                                    low_intervals)
from flucut_posthoc.overview import write_overview
from flucut_posthoc.plots import coverage_plot, length_plot
from flucut_posthoc.report import write_report
from flucut_posthoc.segment_coverage import write_segment_coverage
from report_example import example_subtype, example_summary


def test_states_of_a_measure():
    assert grade("same", 0.995, 100) == "ok"
    assert grade("same", 0.96, 100) == "warn"
    assert grade("same", 0.5, 100) == "bad"
    assert grade("same", 0.5, 5) == "na"       # too few cases
    assert grade("primer_left", 0.01, 100) == "ok"
    assert grade("primer_left", 0.5, 100) == "bad"
    assert grade("depth_ratio", 1.05) == "ok"
    assert grade("depth_ratio", 0.8) == "warn"
    assert grade("depth_ratio", None) == "na"


def test_depth_of_a_segment():
    profiles = {"terminal+": [60] * 80 + [10] * 20,
                "terminal-": [20] * 100, "internal+": [0] * 100,
                "deletion": [20] * 50 + [0] * 50}
    stats = depth_stats(profiles, MEASURED_ORIGIN, 50, 0.9)
    assert stats["state"] == "partial" and stats["well_covered"] == 0.8
    assert stats["covered"] == 1.0 and stats["internal"] == 0.0
    assert stats["deletion"] == 1000 / 8000
    assert stats["minus_strand"] == 2000 / 7000
    assert low_intervals(profiles, 50) == [(80, 100)]
    assert depth_stats({"terminal+": [0] * 10}, MEASURED_ORIGIN)["state"] \
        == "missing"
    assert depth_stats({"terminal+": [50] * 10}, MEASURED_ORIGIN)["state"] \
        == "complete"


def test_checks_of_a_sample():
    s = example_summary()
    checks = {(c.label, c.end): c for c in sample_checks(s)}
    same = checks["same origin", "left"]
    assert same.cases == 997 and same.share == 995 / 997
    found = checks["internally primed ends found", "left"]
    assert found.cases == 102 and found.status == "ok"
    assert checks["primer bases left", None].status == "ok"
    origin = checks["same read origin", None]
    assert origin.cases == 1020 and origin.share == 990 / 1020
    table = origin_table(s["flucut"]["category"])
    assert table["full"]["del"] == 20 and table["oth"]["full"] == 10
    split = checks["same number of amplicons in the read", None]
    assert split.cases == 1092 and split.share == 1090 / 1092
    rows = segment_rows(s)
    assert rows[0]["same_state"] and rows[0]["status"] == "ok"
    assert verdict(list(checks.values()), rows)[0] == "ok"
    assert verdict(list(checks.values()), rows,
                   example_subtype("H7N1"))[0] == "bad"


def test_depth_which_differs():
    s = example_summary(scale=0.5)
    rows = segment_rows(s)
    assert rows[0]["status"] == "bad" and rows[0]["ratio"] < 0.6
    status, text = verdict(sample_checks(s), rows)
    assert status == "bad" and "outside their limits" in text


def test_without_flucut():
    s = example_summary()
    s["flucut"] = None
    for cover in s["coverage"].values():
        cover["predicted"] = cover["agreement"] = None
    for entry in s["complete_amplicons"]:
        entry["flucut"] = None
    assert sample_checks(s) == []
    assert segment_rows(s)[0]["flucut"] is None


def test_plots():
    depth = {"full": np.full(400, 100.0), "del": np.full(400, 20.0)}
    svg = coverage_plot(depth, 400, 150, 50, [(10, 30)], {"IRMA": [90] * 400},
                        site=201)
    assert svg.count("<polygon") == 2 and 'class="c-low"' in svg
    assert ">201</text>" in svg and 'class="depth d0"' in svg
    svg = length_plot([0] * 19 + [5] + [0] * 19 + [20], [0] * 39 + [18],
                      10, 400)
    assert svg.count('class="c-short"') == 1
    assert svg.count('class="c-full"') == 1 and "polyline" in svg


def test_report_of_a_sample(tmp_path):
    path = tmp_path / "posthoc_report.html"
    write_report(example_summary(), str(path))
    page = path.read_text()
    for part in ("flucut and the mapping", "Agreement by measure",
                 "Agreement by segment", "The sample, measured",
                 "Internal priming sites", "measured (mapping)",
                 "predicted (flucut)", "All numbers", "Terms used"):
        assert part in page
    assert "✓ flucut agrees with the mapping" in page
    assert "✓ Subtype H5N1: confirmed by the mapping" in page
    assert page.count('class="seg"') == 1


def test_report_without_flucut(tmp_path):
    s = example_summary()
    s["flucut"] = None
    for cover in s["coverage"].values():
        cover["predicted"] = cover["agreement"] = None
    for entry in s["complete_amplicons"]:
        entry["flucut"] = None
    path = tmp_path / "posthoc_report.html"
    write_report(s, str(path))
    page = path.read_text()
    assert "Agreement by measure" not in page
    assert "predicted (flucut)" not in page
    assert "The sample, measured" in page


def _summaries(tmp_path) -> list[str]:
    paths = []
    for name, scale in (("one", 1.0), ("two", 0.5)):
        directory = tmp_path / name
        directory.mkdir()
        path = directory / "posthoc_summary.json"
        path.write_text(json.dumps(example_summary(scale, name)))
        paths.append(str(path))
    return paths


def test_overview(tmp_path):
    write_overview(_summaries(tmp_path), str(tmp_path), reports=True)
    lines = (tmp_path / "overview.tsv").read_text().splitlines()
    assert lines[0].split("\t")[:2] == ["sample", "reads"]
    assert len(lines) == 3 and lines[1].startswith("one\t1200\t")
    page = (tmp_path / "overview.html").read_text()
    assert '<a href="one/posthoc_report.html">one</a>' in page
    assert "✓ agrees" in page and "✕ differs" in page
    assert lines[0].split("\t")[-4:] == [
        "subtype_flucut", "subtype_mapping", "subtype_irma",
        "subtype_status"]
    assert lines[1].split("\t")[-4:] == ["H5N1", "H5N1", "H5N1", "same"]
    assert page.count("✓ same") == 2
    assert (tmp_path / "two" / "posthoc_report.html").exists()


def test_segment_coverage_page(tmp_path):
    """The page of one segment has a strip per sample."""
    paths = _summaries(tmp_path)
    page = open(write_segment_coverage(paths, str(tmp_path), "HA")).read()
    assert page.count('class="seg"') == 2
    assert "<b>one</b>" in page and "<b>two</b>" in page
    assert page.count("predicted (flucut)") == 2
    page = open(write_segment_coverage(paths, str(tmp_path), "NA")).read()
    assert page.count('class="seg"') == 0
    assert "No reads map to this segment in: one, two." in page
