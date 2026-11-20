from flucut_posthoc.html import pct
from flucut_posthoc.plots import bars


def test_zero_and_all_only_when_exact():
    assert pct(0, 2000) == "0.0%"
    assert pct(1, 2000) == "<0.1%"
    assert pct(1, 1000) == "0.1%"
    assert pct(423, 1000) == "42.3%"
    assert pct(999, 1000) == "99.9%"
    assert pct(1999, 2000) == ">99.9%"
    assert pct(2000, 2000) == "100.0%"
    assert pct(5, 0) == "–"


def test_other_numbers_of_decimals():
    assert pct(1, 100000, 2) == "<0.01%"
    assert pct(99999, 100000, 2) == ">99.99%"
    assert pct(0.99999, 1, 2) == ">99.99%"
    assert pct(1.0, 1, 2) == "100.00%"
    assert pct(1, 1000, 0) == "<1%"
    assert pct(999, 1000, 0) == ">99%"
    assert pct(1000, 1000, 0) == "100%"


def test_ratio_above_the_total():
    assert pct(100001, 100000) == ">100.0%"
    assert pct(3, 2) == "150.0%"


def test_bar_labels():
    html = bars([("most", 9999, "#000"), ("few", 1, "#fff")], 10000)
    assert "<em>&gt;99.9%</em>" in html
    assert "<em>&lt;0.1%</em>" in html
