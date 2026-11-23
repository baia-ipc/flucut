from flucut_posthoc.html import (Markup, badge, fold, pct, stacked_bar,
                                 table)


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
    html = stacked_bar("Reads", "reads", [("c-full", "most", 9999),
                                          ("c-oth", "few", 1)])
    assert "of 10,000 reads" in html
    assert "most <b>&gt;99.9%</b>" in html
    assert "few <b>&lt;0.1%</b>" in html
    assert stacked_bar("Reads", "reads", [("c-full", "none", 0)]) == ""


def test_badge_has_a_symbol_and_a_word():
    assert badge("ok", "agrees") == '<span class="st ok">✓ agrees</span>'
    assert badge("warn", "differs").count("! differs") == 1


def test_table_cells():
    html = table(["a", "b"], [["x<y", 1234], [Markup("<b>z</b>"), None]],
                 groups=[("both", 2)])
    assert "x&lt;y" in html and "1,234" in html
    assert "<b>z</b>" in html and "<td>–</td>" in html
    assert '<th colspan="2" class="g">both</th>' in html
    assert fold("Title", "") == ""
