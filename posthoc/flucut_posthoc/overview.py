"""An overview of several samples: one line per sample, from the
summaries of their runs (posthoc_summary.json)."""

import json
import os
from html import escape

from .html import STYLE, num, para, pct
from .report import write_report

COLUMNS = [
    ("sample", "sample"),
    ("reads", "reads"),
    ("influenza", "map to influenza"),
    ("offtarget", "other amplicons"),
    ("amplicons", "influenza amplicons"),
    ("deletion", "with an internal deletion"),
    ("internal", "5′ ends internally primed"),
    ("sites", "internal priming sites"),
    ("origin_agree", "flucut: origin calls confirmed"),
    ("internal_confirmed", "flucut: internal calls confirmed"),
    ("internal_found", "flucut: internal ends found"),
    ("terminal_found", "flucut: terminal ends found"),
    ("primer_left", "flucut: internal ends with primer bases left"),
    ("depth_ratio", "coverage: predicted / measured"),
    ("depth_difference", "coverage: mean difference per position"),
]


def _ratio(part: float, total: float) -> float | None:
    return part / total if total else None


def sample_row(s: dict) -> dict:
    """The line of a sample from its summary."""
    reads, classes = s["reads"]["total"], s["reads"]["classes"]
    amplicons = s["amplicons"]
    left = s["ends"].get("left", {})
    primed = sum(sum(v.values()) for k, v in left.items() if k != "none")
    internal = sum(v.get("internal", 0) for k, v in left.items()
                   if k != "none")
    row = {"sample": s["sample"], "reads": reads,
           "influenza": _ratio(classes.get("influenza", 0), reads),
           "offtarget": _ratio(classes.get("offtarget", 0), reads),
           "amplicons": amplicons["total"],
           "deletion": _ratio(amplicons["categories"].get("deletion", 0),
                              amplicons["total"]),
           "internal": _ratio(internal, primed),
           "sites": s["internal_sites_total"]}
    f = s["flucut"]
    if f:
        scores = [sc for sc in f["scores"].values() if sc["called"]]
        row["origin_agree"] = _ratio(sum(sc["agree"] for sc in scores),
                                     sum(sc["called"] for sc in scores))
        left_scores = f["scores"].get("left", {})
        row["internal_confirmed"] = left_scores.get("internal_precision")
        row["internal_found"] = left_scores.get("internal_recall")
        row["terminal_found"] = left_scores.get("terminal_recall")
        trim = f["trim"].get("internal", {})
        row["primer_left"] = _ratio(trim.get("primer bases left", 0),
                                    sum(trim.values()))
    agreements = [(r["length"], s["coverage"][r["name"]]["agreement"])
                  for r in s["references"]]
    agreements = [(n, a) for n, a in agreements if a]
    if agreements:
        measured = sum(n * a["measured_mean_depth"] for n, a in agreements)
        row["depth_ratio"] = _ratio(
            sum(n * a["predicted_mean_depth"] for n, a in agreements),
            measured)
        row["depth_difference"] = _ratio(
            sum(n * a["mean_absolute_difference"] for n, a in agreements),
            measured)
    return row


def _cell(key: str, value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_overview(summaries: list[str], outdir: str,
                   reports: bool = False) -> None:
    """Write overview.tsv and overview.html for the summaries given (the
    paths of posthoc_summary.json files); with `reports`, also write the
    report of each sample again from its summary."""
    rows, links = [], []
    for path in summaries:
        with open(path) as handle:
            summary = json.load(handle)
        rows.append(sample_row(summary))
        report = os.path.join(os.path.dirname(path), "posthoc_report.html")
        if reports:
            write_report(summary, report)
        links.append(os.path.relpath(report, outdir))
    with open(os.path.join(outdir, "overview.tsv"), "w") as out:
        out.write("\t".join(key for key, _ in COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(_cell(key, row.get(key))
                                for key, _ in COLUMNS) + "\n")
    body = []
    for row, link in zip(rows, links):
        cells = []
        for key, _ in COLUMNS:
            value = row.get(key)
            if key == "sample":
                cells.append(f'<td class="l"><a href="{escape(link)}">'
                             f"{escape(value)}</a></td>")
            elif isinstance(value, float):
                share = pct(value, 1, 2 if value > 0.9 else 1)
                cells.append(f"<td>{escape(share)}</td>")
            else:
                cells.append(f"<td>{num(value)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    head = "".join(f"<th>{escape(label)}</th>" for _, label in COLUMNS)
    page = (
        '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
        "<title>flucut post-hoc: overview</title>"
        f"<style>{STYLE} th {{ white-space:normal; max-width:110px; "
        "vertical-align:bottom; }} main {{ max-width:1500px; }}</style>"
        "</head><body><main><h1>flucut post-hoc: overview</h1>"
        + para("One line per sample, measured from the mapping of every "
               "read; the name of a sample links to its report. The "
               "columns of flucut compare its predictions with the "
               "mapping: origin calls confirmed (of the ends flucut "
               "calls terminal or internal, those the mapping places "
               "there), internal calls confirmed, internal and terminal "
               "ends found (of those of the mapping, the ones flucut "
               "calls so), internally primed ends at which flucut left "
               "primer bases, and its coverage plot against the depth "
               "measured (all segments together).")
        + f'<div class="scroll"><table><tr>{head}</tr>{"".join(body)}'
        "</table></div></main></body></html>")
    with open(os.path.join(outdir, "overview.html"), "w") as out:
        out.write(page)
