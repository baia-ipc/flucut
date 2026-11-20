from flucut_posthoc.irma import _deletions, _primers_inside, gene_of
from flucut_posthoc.offtarget import cluster_offtarget
from flucut_posthoc.summary import _footprint, cluster_sites

from conftest import UNI12, random_seq


def test_positions_grouped_into_sites():
    sites = cluster_sites({100: 3, 102: 10, 107: 1, 113: 4})
    assert [(s["position"], s["from"], s["to"], s["reads"])
            for s in sites] == [(102, 100, 107, 14), (113, 113, 113, 4)]


def test_template_under_the_primer(reference):
    seq = reference.seq
    # left end at 1-based position 301: the 12 bases before it
    template, mismatches = _footprint(reference, "left", 301,
                                      "AGCRAAAGCAGG")
    assert template == seq[288:300]
    # the terminus itself has no mismatch with the primer
    assert _footprint(reference, "left", 13, "AGCRAAAGCAGG") == (UNI12, 0)
    assert _footprint(reference, "right", 587, "AGTAGAAACAAGG")[1] == 0


def test_irma_names_and_consensus_checks():
    assert gene_of("A_HA_H5") == "HA" and gene_of("A_MP") == "M"
    assert _deletions([(0, 200, 30, 230), (200, 440, 2100, 2340)]) == \
        [[231, 2100]]
    assert _deletions([(0, 200, 30, 230), (350, 500, 380, 530)]) == []
    joined = random_seq(200, 40) + "CCTTGTTTCTACT" + random_seq(30, 41) \
        + UNI12 + random_seq(200, 42)
    assert _primers_inside(joined) == 2
    assert _primers_inside(random_seq(400, 43)) == 0


def test_offtarget_groups():
    one, two = random_seq(200, 50), random_seq(150, 51)
    samples = [(one[:24], 200, True, one)] * 8 + \
        [(two[:24], 150, True, two)] * 2
    clusters = cluster_offtarget(samples)
    assert [c["sampled_reads"] for c in clusters] == [8, 2]
    assert clusters[0]["representative"] == one
    assert clusters[0]["template_length"] == 200


def test_segment_coverage_page(tmp_path):
    """The page of one segment has a pair of plots per sample."""
    import json
    from flucut_posthoc.segment_coverage import write_segment_coverage
    paths = []
    for name in ("one", "two"):
        cover = {"measured": {"terminal+": [3, 4, 4, 2]},
                 "predicted": {"terminal+": [3.0, 3.0, 4.0, 2.0]},
                 "irma": {}, "agreement": {
                     "measured_mean_depth": 3.25,
                     "predicted_mean_depth": 3.0, "correlation": 0.9,
                     "mean_absolute_difference": 0.25}}
        summary = {"sample": name, "internal_sites": [],
                   "references": [{"name": "r", "gene": "HA"}],
                   "coverage": {"r": cover}}
        path = tmp_path / f"{name}.json"
        path.write_text(json.dumps(summary))
        paths.append(str(path))
    page = open(write_segment_coverage(paths, str(tmp_path), "HA")).read()
    assert page.count('class="pair"') == 2
    assert "one measured" in page and "two as flucut predicts" in page
    assert open(write_segment_coverage(paths, str(tmp_path), "NA")).read() \
        .count('class="pair"') == 0
