"""The coverage of a segment, measured and as flucut predicts it: the
pair of plots of a sample, and a page with the pairs of several
samples."""

import json
import os
from html import escape

from .html import STYLE, para, pct, table
from .plots import MEASURED_BANDS, PREDICTED_BANDS, legend, stacked_coverage


def largest_sites(s: dict) -> dict[str, int]:
    """The position of the largest internal priming site of each
    reference."""
    largest = {}
    for site in s["internal_sites"]:
        largest.setdefault(site["reference"], site["position"])
    return largest


def coverage_pair(s: dict, ref: dict, title: str,
                  irma_lines: bool = True) -> str | None:
    """The measured coverage of a reference and, at its right and on the
    same scale, the one flucut predicts (None: nothing maps there)."""
    cover = s["coverage"][ref["name"]]
    ymax = max(sum(v) for v in zip(*cover["measured"].values()))
    predicted = cover["predicted"]
    if predicted:
        ymax = max(ymax, max(sum(v) for v in zip(*predicted.values())))
    if ymax == 0:
        return None
    lines = cover["irma"] if irma_lines else {}
    for depth in lines.values():
        ymax = max(ymax, max(depth))
    plot = (f'<div class="plot"><b>{escape(title)} measured</b>'
            + stacked_coverage(cover["measured"], MEASURED_BANDS, ymax,
                               lines, largest_sites(s).get(ref["name"]))
            + "</div>")
    if predicted:
        plot += (f'<div class="plot"><b>{escape(title)} as flucut '
                 f"predicts</b>"
                 + stacked_coverage(predicted, PREDICTED_BANDS, ymax)
                 + "</div>")
    return f'<div class="pair">{plot}</div>'


def agreement_row(label: str, agreement: dict) -> list:
    """The line of a segment in a table of predicted and measured
    depth."""
    r = agreement["correlation"]
    return [label, agreement["measured_mean_depth"],
            agreement["predicted_mean_depth"],
            pct(agreement["predicted_mean_depth"],
                agreement["measured_mean_depth"]),
            f"{r:.3f}" if r is not None else "–",
            pct(agreement["mean_absolute_difference"],
                agreement["measured_mean_depth"])]


AGREEMENT_COLUMNS = ["measured depth", "predicted depth",
                     "predicted / measured", "correlation",
                     "difference per position"]


def write_segment_coverage(summaries: list[str], outdir: str,
                           gene: str) -> str:
    """Write coverage_<gene>.html: the coverage of one segment in each
    sample, from the summaries given. Returns the path."""
    pairs, rows, missing = [], [], []
    for path in summaries:
        with open(path) as handle:
            s = json.load(handle)
        ref = next((r for r in s["references"]
                    if r["gene"].upper() == gene.upper()), None)
        pair = ref and coverage_pair(s, ref, s["sample"], irma_lines=False)
        if not pair:
            missing.append(s["sample"])
            continue
        pairs.append(pair)
        agreement = s["coverage"][ref["name"]]["agreement"]
        if agreement:
            rows.append(agreement_row(s["sample"], agreement))
    out = [para(
        f"The depth at each position of {gene} in each sample. At the "
        "left, measured from the mappings (the bases of a read aligned "
        "there; a deletion of a read is not covered), stacked by what "
        "the amplicon is: where its 5′ (left) end is and whether a "
        "primer is there. The dashed red line is the largest internal "
        "priming site. At the right, on the same scale, the coverage "
        "which flucut predicts."),
        "<b>Measured</b>", legend(MEASURED_BANDS),
        "<b>As flucut predicts</b>", legend(PREDICTED_BANDS)]
    out.extend(pairs)
    if missing:
        out.append(para("No reads map to this segment in: "
                        + ", ".join(missing) + ".", "note"))
    if rows:
        out.append("<h3>Predicted and measured depth</h3>")
        out.append(table(["sample"] + AGREEMENT_COLUMNS, rows))
        out.append(para(
            "The depths are means over the segment. The difference per "
            "position is the mean of the absolute differences between "
            "the two depths, as a share of the measured mean depth.",
            "note"))
    title = f"flucut post-hoc: coverage of {gene}"
    page = ('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{escape(title)}</title><style>{STYLE}</style></head>"
            f"<body><main><h1>{escape(title)}</h1>{''.join(out)}"
            "</main></body></html>")
    target = os.path.join(outdir, f"coverage_{gene}.html")
    with open(target, "w") as handle:
        handle.write(page)
    return target
