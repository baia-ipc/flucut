"""One segment in several samples: the strip of the segment of each
sample (numbers, coverage measured and predicted, amplicon lengths) on
one page."""

import json
import os

from .agreement import depth_limits, segment_rows
from .html import (Markup, badge, card, heading, note, page, pct, table,
                   terms, write_page)
from .report_segments import segment_strip, strips_head


def _agreement_line(sample: str, row: dict) -> list:
    r = row["correlation"]
    return [sample,
            Markup(badge(row["status"], "agrees" if row["status"] == "ok"
                         else "not compared" if row["status"] == "na"
                         else "differs")),
            f"{row['measured']['mean_depth']:,.0f}×",
            f"{row['flucut']['mean_depth']:,.0f}×"
            if row["flucut"].get("mean_depth") is not None else None,
            pct(row["ratio"], 1) if row["ratio"] is not None else None,
            pct(row["difference"], 1) if row["difference"] is not None
            else None,
            f"{r:.3f}" if r is not None else None]


def write_segment_coverage(summaries: list[str], outdir: str,
                           gene: str) -> str:
    """Write coverage_<gene>.html: the strip of one segment of each
    sample, from the summaries given. Returns the path."""
    strips, lines, missing = [], [], []
    with_flucut, min_depth = False, 50
    for path in summaries:
        with open(path) as handle:
            s = json.load(handle)
        ref = next((r for r in s["references"]
                    if r["gene"].upper() == gene.upper()), None)
        row = ref and next(row for row in segment_rows(s)
                           if row["reference"] == ref["name"])
        strip = ref and segment_strip(s, ref, row, s["sample"],
                                      with_lines=False)
        if not strip:
            missing.append(s["sample"])
            continue
        strips.append(strip)
        min_depth = depth_limits(s)[0]
        if row["flucut"] is not None:
            with_flucut = True
            lines.append(_agreement_line(s["sample"], row))
    title = f"flucut post-hoc: {gene} in all samples"
    parts = [f"<h1><small>flucut post-hoc</small>{gene} in all samples</h1>"]
    if lines:
        parts += [heading("agreement", "flucut and the mapping"),
                  card("Depth of flucut and depth measured", table(
                      ["sample", "", "measured", "flucut",
                       "flucut / measured", "difference per position",
                       "correlation"], lines, left=2,
                      groups=[("", 2), ("mean depth", 2),
                              ("depth along the segment", 3)]))]
    parts.append(card(f"{gene} by sample",
                      strips_head(with_flucut, min_depth), *strips))
    if missing:
        parts.append(note("No reads map to this segment in: "
                          + ", ".join(missing) + "."))
    parts.append(terms([
        ("measured", "From the base-by-base mapping of every read."),
        ("predicted, flucut", "The depth flucut estimates from where it "
         "locates each read, without mapping."),
        ("difference per position", "Mean absolute difference between "
         "the two depths, as a share of the mean depth measured."),
        ("agrees", "Same state of the segment, and depth within the "
         "limits given in the report of the sample.")]))
    target = os.path.join(outdir, f"coverage_{gene}.html")
    write_page(target, page(title, *parts))
    return target
