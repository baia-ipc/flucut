"""An overview of several samples: one line per sample, from the
summaries of their runs (posthoc_summary.json)."""

import json
import os
from html import escape

from .agreement import (depth_agreement, grade, limit_text, sample_checks,
                        segment_rows, verdict, worst)
from .html import (Markup, badge, card, chip, chips, heading, legend, num,
                   page, pct, table, terms, tile, tiles, write_page)
from .report import write_report
from .subtype import STATUS

# (key, label, the measure whose limits give the state of the cell)
COLUMNS = [
    ("sample", "sample", None),
    ("reads", "reads", None),
    ("influenza", "map to influenza", None),
    ("offtarget", "other amplicons", None),
    ("amplicons", "influenza amplicons", None),
    ("deletion", "with an internal deletion", None),
    ("internal", "5′ ends internally primed", None),
    ("sites", "internal priming sites", None),
    ("origin_agree", "flucut: origin calls confirmed", "same"),
    ("internal_confirmed", "flucut: internal calls confirmed", "confirmed"),
    ("internal_found", "flucut: internal ends found", "found"),
    ("terminal_found", "flucut: terminal ends found", "found"),
    ("primer_left", "flucut: internal ends with primer bases left",
     "primer_left"),
    ("depth_ratio", "coverage: predicted / measured", "depth_ratio"),
    ("depth_difference", "coverage: mean difference per position",
     "depth_difference"),
    ("subtype_flucut", "subtype: flucut", None),
    ("subtype_mapping", "subtype: mapping", None),
    ("subtype_irma", "subtype: IRMA", None),
    ("subtype_status", "subtype of flucut: same, different, unknown", None),
]
N_SHARED = 15     # columns shown as they are in the table of the page
MEASURED = ["reads", "map to influenza", "other amplicons", "amplicons",
            "internal deletion", "5′ ends int. primed", "sites"]
SUBTYPE_MARKS = {"same": "same", "different": "different",
                 "unknown": "cannot be told"}
COMPARED = ["origins the same", "internal confirmed", "internal found",
            "terminal found", "primer bases left", "flucut / measured",
            "difference"]


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
    depth = depth_agreement(s)
    if depth:
        row["depth_ratio"] = depth["ratio"]
        row["depth_difference"] = depth["difference"]
    subtype = s.get("subtype")
    if subtype:
        irma = sorted({found for found in subtype["irma"].values() if found})
        row.update(subtype_flucut=subtype["flucut"],
                   subtype_mapping=subtype["measured"],
                   subtype_irma=", ".join(irma) or None,
                   subtype_status=subtype["status"])
    return row


def _cell(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _write_tsv(rows: list[dict], path: str) -> None:
    with open(path, "w") as out:
        out.write("\t".join(key for key, _, _ in COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(_cell(row.get(key))
                                for key, _, _ in COLUMNS) + "\n")


def _line(row: dict, link: str, s: dict) -> tuple[list, str]:
    """The line of a sample in the table, and the state of its
    agreement."""
    segments = segment_rows(s)
    subtype = s.get("subtype")
    status, _ = verdict(sample_checks(s), segments, subtype)
    compared = [seg for seg in segments if seg["flucut"] is not None]
    same = sum(seg["status"] == "ok" for seg in compared)
    label = {"ok": "agrees", "na": "not compared"}.get(status, "differs")
    cells = [Markup(f'<a href="{escape(link)}">{escape(row["sample"])}</a>'),
             Markup(badge(status, label)),
             Markup(badge(STATUS[subtype["status"]],
                          SUBTYPE_MARKS[subtype["status"]]))
             if subtype else None,
             row.get("subtype_flucut"), row.get("subtype_mapping"),
             row.get("subtype_irma"),
             f"{same}/{len(compared)}" if compared else None]
    for key, _, kind in COLUMNS[1:N_SHARED]:
        value = row.get(key)
        if isinstance(value, float):
            text = pct(value, 1, 2 if value > 0.9 and kind else 1)
        else:
            text = num(value)
        cells.append((text, f"sc {grade(kind, value)}")
                     if kind and value is not None else text)
    return cells, status


def _page(lines: list[list], states: list[str]) -> str:
    counts = {state: states.count(state) for state in set(states)}
    differ = counts.get("warn", 0) + counts.get("bad", 0)
    summary = [f'<b class="count">{len(states)} samples</b> ']
    if counts.get("ok"):
        summary.append(badge("ok", f"{counts['ok']} agree with flucut"))
    if differ:
        summary.append(badge(worst(states),
                             f"{differ} with checks outside their limits"))
    if counts.get("na"):
        summary.append(chip(f"{counts['na']} not compared with flucut"))
    headers = ["sample", "", "", "flucut", "mapping", "IRMA",
               "segments agree"] + MEASURED + COMPARED
    limits = [(label, limit_text(kind)) for label, (_, _, kind) in
              zip(COMPARED, COLUMNS[8:N_SHARED])]
    return page(
        "flucut post-hoc: overview",
        "<h1><small>flucut post-hoc</small>Overview</h1>",
        heading("samples", "flucut and the mapping, by sample"),
        chips(*summary),
        tiles(tile(str(len(states)), "samples"),
              tile(str(counts.get("ok", 0)), "agree in all checks",
                   status="ok" if counts.get("ok") else ""),
              tile(str(differ), "with checks outside their limits",
                   status=worst(states) if differ else "")),
        card("Samples", table(
            headers, lines, left=2, cls="samples",
            groups=[("", 2), ("subtype", 4), ("", 1),
                    ("measured from the mapping", len(MEASURED)),
                    ("flucut read ends", 5), ("flucut depth", 2)]),
            legend([("t-ok", "within the limit"),
                    ("t-warn", "within the wider limit"),
                    ("t-bad", "outside"), ("t-na", "too few cases")])),
        terms([("sample", "Links to the report of the sample, which has "
                "every check and its limits."),
               ("subtype", "Of flucut, of the mapping (the H of the HA "
                "reference and the N of the NA reference the reads map "
                "to best) and of the IRMA runs. ✓ the same for flucut "
                "and the mapping, ✕ different, · cannot be told (too few "
                "reads on HA or NA, or none assigned)."),
               ("segments agree", "Segments in the same state for flucut "
                "and the mapping, with the depth within the limits."),
               ("measured", "Shares of the reads (map to influenza, "
                "other amplicons), of the influenza amplicons (internal "
                "deletion), of the 5′ ends with a primer.")]
              + [(label, f"✓ {text}") for label, text in limits]
              + [("read ends", "Origins the same: of the ends both call "
                  "terminal or internal. Confirmed: of the 5′ ends flucut "
                  "calls internal. Found: of the 5′ ends the mapping "
                  "places there. Primer bases left: of the internally "
                  "primed ends trimmed."),
                 ("depth", "All segments together; difference: mean "
                  "absolute difference per position, of the depth "
                  "measured.")]), wide=True)


def write_overview(summaries: list[str], outdir: str,
                   reports: bool = False) -> None:
    """Write overview.tsv and overview.html for the summaries given (the
    paths of posthoc_summary.json files); with `reports`, also write the
    report of each sample again from its summary."""
    rows, lines, states = [], [], []
    for path in summaries:
        with open(path) as handle:
            summary = json.load(handle)
        report = os.path.join(os.path.dirname(path), "posthoc_report.html")
        if reports:
            write_report(summary, report)
        row = sample_row(summary)
        line, status = _line(row, os.path.relpath(report, outdir), summary)
        rows.append(row)
        lines.append(line)
        states.append(status)
    _write_tsv(rows, os.path.join(outdir, "overview.tsv"))
    write_page(os.path.join(outdir, "overview.html"), _page(lines, states))
