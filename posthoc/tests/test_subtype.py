"""The subtype of the mapping, of flucut and of IRMA, and the verdict."""

from flucut_posthoc.report_subtype import subtype_card
from flucut_posthoc.subtype import irma_subtype, part, subtype_of
from report_example import example_subtype


def test_subtype_named_by_a_reference():
    assert subtype_of("Influenza A virus (A/swan/Germany/1/2025(H5N1)) "
                      "segment 4 hemagglutinin (HA) gene") == "H5N1"
    assert subtype_of("A/x/1/2020(H10N8) segment 6") == "H10N8"
    assert subtype_of("test virus segment 4") is None
    assert part("H10N8", "H") == "H10" and part("H5N?", "N") is None
    assert part(None, "H") is None


def test_subtype_of_irma():
    assert irma_subtype([["4-A_HA_H5", 400], ["4-A_NA_N1", 250],
                         ["4-A_HA_H3", 20], ["5-A_NA_N8", 900]]) == "H5N1"
    assert irma_subtype([["4-A_HA_H5", 400]]) == "H5N?"
    assert irma_subtype([["1-initial", 1000]]) is None


def test_subtype_confirmed():
    subtype = example_subtype("H5N1", other_h=12)
    assert subtype["status"] == "same"
    assert subtype["text"] == "Subtype H5N1: confirmed by the mapping"
    assert subtype["measured"] == "H5N1" and subtype["irma"] == {
        "IRMA": "H5N1"}
    ha = subtype["genes"]["HA"]
    assert ha["reads"] == 500 and ha["other"] == {"type": "H7", "reads": 12}
    # another reference for an internal gene is not another subtype
    pb2 = subtype["segments"][0]
    assert pb2["same_reference"] is False
    assert pb2["reference_subtype"] == "H7N9"
    assert pb2["flucut_reference_subtype"] == "H5N1"


def test_subtype_different():
    subtype = example_subtype("H7N1")
    assert subtype["status"] == "different"
    assert subtype["text"] == "Subtype: flucut H7N1, mapping H5N1"
    assert subtype["genes"]["HA"]["status"] == "different"
    assert subtype["genes"]["NA"]["status"] == "same"


def test_subtype_cannot_be_told():
    # too few reads on HA: not confirmed, although the H is the same
    subtype = example_subtype("H5N1", ha_reads=5)
    assert subtype["status"] == "unknown"
    assert subtype["text"] == \
        "Subtype: cannot be told (HA 5 mapped reads)"
    assert subtype["measured"] == "H?N1"
    assert subtype["genes"]["HA"]["closest"] == "H5"
    # flucut assigned none, or no N
    assert example_subtype(None)["status"] == "unknown"
    assert example_subtype("H5N?")["status"] == "unknown"
    # a wrong N is different, whatever the reads on HA
    assert example_subtype("H5N2", ha_reads=5)["status"] == "different"


def test_subtype_card():
    page = subtype_card({"subtype": example_subtype("H7N1", other_h=12)})
    assert "✕ Subtype: flucut H7N1, mapping H5N1" in page
    assert "✕ different" in page and "✓ same" in page
    assert "H7: 12" in page and "· other reference" in page
    assert "ha5 (H5N1)" in page
    assert subtype_card({}) == ""
