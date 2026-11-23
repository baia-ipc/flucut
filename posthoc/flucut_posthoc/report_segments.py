"""The segments of a sample, one strip for each: the numbers measured
next to those of flucut, the coverage by read origin measured above the
one flucut predicts, and the lengths of the complete amplicons."""

from html import escape

from .agreement import depth_limits, segment_rows
from .html import Markup, badge, card, legend, note, num, pct
from .origins import (MEASURED_ORIGIN, ORIGINS, PREDICTED_ORIGIN,
                      STATE_LABELS, STATE_STATUS, by_origin, low_intervals)
from .plots import coverage_plot, length_plot

AGREEMENT_LABELS = {"ok": "agrees", "warn": "differs", "bad": "differs",
                    "na": "not compared"}
LINE_CLASSES = ["c-line", "c-dash", "c-dash"]


def largest_sites(s: dict) -> dict[str, int]:
    """The position of the largest internal priming site of each
    reference."""
    largest = {}
    for site in s["internal_sites"]:
        largest.setdefault(site["reference"], site["position"])
    return largest


def _state(stats: dict | None) -> Markup | None:
    if stats is None:
        return None
    return Markup(badge(STATE_STATUS[stats["state"]],
                        STATE_LABELS[stats["state"]]))


def _share(stats: dict | None, key: str) -> str | None:
    if stats is None or stats.get(key) is None:
        return None
    return pct(stats[key], 1)


def _depth(stats: dict | None) -> str | None:
    if stats is None or stats.get("median_depth") is None:
        return None
    return f"{stats['median_depth']:,.0f}×"


def _numbers(row: dict, complete: dict, min_depth: int) -> str:
    """The numbers of a segment: measured and, if given, of flucut."""
    measured, flucut = row["measured"], row["flucut"]
    said = complete.get("flucut") or {}
    both = flucut is not None or bool(said)
    lines = [
        ("state", _state(measured), _state(flucut)),
        ("median depth", _depth(measured), _depth(flucut)),
        (f"at ≥{min_depth}×", _share(measured, "well_covered"),
         _share(flucut, "well_covered")),
        ("covered", _share(measured, "covered"), _share(flucut, "covered")),
        ("share of the depth",),
        ("internally primed", _share(measured, "internal"),
         _share(flucut, "internal")),
        ("internal deletion", _share(measured, "deletion"),
         _share(flucut, "deletion")),
        ("minus strand", _share(measured, "minus_strand"),
         _share(flucut, "minus_strand")),
        ("complete amplicons",),
        ("number", num(complete["amplicons"]),
         num(said["amplicons"]) if said else None),
        ("short", pct(complete["short"], complete["amplicons"]),
         pct(said["short"], said["amplicons"]) if said else None)]
    out = ['<table><tr><th class="l"></th><th>measured</th>'
           + ("<th>flucut</th>" if both else "") + "</tr>"]
    for line in lines:
        if len(line) == 1:
            out.append(f'<tr><td class="kh" colspan="{3 if both else 2}">'
                       f"{line[0]}</td></tr>")
            continue
        cells = "".join(f"<td>{value if value is not None else '–'}</td>"
                        for value in line[1:(3 if both else 2)])
        out.append(f'<tr><td class="l">{escape(line[0])}</td>{cells}</tr>')
    return "".join(out) + "</table>"


def _notes(ref: dict, row: dict) -> str:
    """One line on the agreement of the depth, one for each consensus of
    IRMA."""
    out = []
    if row["ratio"] is not None:
        r = row["correlation"]
        out.append(note(
            f"depth: flucut {pct(row['ratio'], 1)} of measured, "
            f"{pct(row['difference'], 1)} difference per position"
            + (f", r {r:.3f}" if r is not None else "")))
    for label, info in ref.get("irma", {}).items():
        text = (f"{label}: consensus {info['consensus_length']:,} b, "
                f"{pct(info['reference_covered'], ref['length'], 0)} of "
                f"the reference")
        if info["deletions"]:
            text += "; deleted " + ", ".join(
                f"{a}–{b}" for a, b in info["deletions"])
        flagged = info["reference_covered"] < 0.9 * ref["length"]
        if info["primers_inside"]:
            text += f"; {info['primers_inside']} primers inside"
            flagged = True
        out.append(f'<p class="note">{"! " if flagged else ""}'
                   f"{escape(text)}</p>")
    return "".join(out)


def _coverage(s: dict, ref: dict, with_lines: bool) -> str | None:
    """The coverage measured and, below it and on the same scale, the
    one flucut predicts (None: nothing maps to the reference)."""
    cover = s["coverage"][ref["name"]]
    min_depth, _ = depth_limits(s)
    measured = by_origin(cover["measured"], MEASURED_ORIGIN)
    ymax = float(sum(measured.values()).max())
    predicted = None
    if cover.get("predicted"):
        predicted = by_origin(cover["predicted"], PREDICTED_ORIGIN)
        ymax = max(ymax, float(sum(predicted.values()).max()))
    if ymax == 0:
        return None
    lines = cover.get("irma", {}) if with_lines else {}
    for depth in lines.values():
        ymax = max(ymax, max(depth))
    out = coverage_plot(measured, ref["length"], ymax, min_depth,
                        low_intervals(cover["measured"], min_depth), lines,
                        largest_sites(s).get(ref["name"]))
    if predicted is None:
        return out
    return ('<p class="note">measured (mapping)</p>' + out
            + '<p class="note">predicted (flucut)</p>'
            + coverage_plot(predicted, ref["length"], ymax, min_depth,
                            low_intervals(cover["predicted"], min_depth)))


def _lengths(s: dict, complete: dict) -> str:
    said = complete.get("flucut") or {}
    settings = s.get("settings", {})
    text = (f"{num(complete['amplicons'])} complete amplicons, "
            f"{pct(complete['short'], complete['amplicons'])} short")
    if said:
        text += (f" · flucut {num(said['amplicons'])}, "
                 f"{pct(said['short'], said['amplicons'])} short")
    if not complete["amplicons"] and not said.get("amplicons"):
        return note("no complete amplicons")
    return length_plot(complete["histogram"], said.get("histogram"),
                       settings.get("length_bin", 10),
                       complete["reference_length"],
                       settings.get("short_fraction", 0.9)) + note(text)


def segment_strip(s: dict, ref: dict, row: dict, title: str | None = None,
                  with_lines: bool = True) -> str | None:
    """The strip of a reference (None: nothing maps to it); `row`: its
    entry of agreement.segment_rows(); `title`: shown instead of the
    name of the segment."""
    coverage = _coverage(s, ref, with_lines)
    if coverage is None:
        return None
    complete = next(c for c in s["complete_amplicons"]
                    if c["reference"] == ref["name"])
    min_depth, _ = depth_limits(s)
    state = row["measured"]["state"]
    mark = badge(row["status"], AGREEMENT_LABELS[row["status"]]) \
        if row["flucut"] is not None else \
        badge(STATE_STATUS[state], STATE_LABELS[state])
    identity = f" · identity {pct(ref['identity'], 1, 2)}" \
        if ref.get("identity") else ""
    return (
        f'<div class="seg" id="seg-{escape(ref["gene"])}">'
        f'<div class="sinfo"><div class="sname">'
        f'<b>{escape(title or ref["gene"])}</b>{mark}</div>'
        + note(f"{ref['name']} · {ref['length']:,} b{identity}")
        + _numbers(row, complete, min_depth) + _notes(ref, row)
        + f'</div><div class="scov">{coverage}</div>'
        f'<div class="slen">{_lengths(s, complete)}</div></div>')


def strips_head(with_flucut: bool, min_depth: int,
                lines: list[str] = ()) -> str:
    """The titles and the legends above the strips."""
    cover = [(f"c-{key}", label) for key, label in ORIGINS[1:] + ORIGINS[:1]]
    cover.append(("c-low", f"below {min_depth}×"))
    cover.append(("c-site", "largest internal priming site"))
    cover += [(LINE_CLASSES[k % 3], f"depth of {label}")
              for k, label in enumerate(lines)]
    lengths = [("c-full", "full length"), ("c-short", "short")]
    if with_flucut:
        lengths.append(("c-line", "flucut"))
    return ('<div class="seghead"><span></span><div>'
            f"<b>Coverage by read origin</b>{legend(cover)}</div><div>"
            f"<b>Complete amplicons by length (bases)</b>{legend(lengths)}"
            "</div></div>")


def segments_card(s: dict) -> str:
    """The strips of all segments of a sample."""
    rows = {row["reference"]: row for row in segment_rows(s)}
    strips, missing = [], []
    for ref in s["references"]:
        strip = segment_strip(s, ref, rows[ref["name"]])
        if strip:
            strips.append(strip)
        else:
            missing.append(ref["gene"])
    with_flucut = any(row["flucut"] is not None for row in rows.values())
    parts = [strips_head(with_flucut, depth_limits(s)[0],
                         [run["label"] for run in s["irma"]]), *strips]
    if missing:
        parts.append(note("No reads map to: " + ", ".join(missing) + "."))
    return card("Segments", *parts)
