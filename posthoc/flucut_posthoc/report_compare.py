"""The parts of the report which compare flucut and IRMA with the
mapping."""

from html import escape

from .html import cross_table, num, para, pct, section, table
from .plots import MEASURED_BANDS

CATEGORY_LABELS = {key: label for key, label, _ in MEASURED_BANDS}
CATEGORY_LABELS.update({
    "several": "several amplicons in the read",
    "dimer": "primer dimer", "offtarget": "other amplicon, not mapped",
    "one_primer": "not mapped, one primer",
    "no_primer": "not mapped, no primer"})
ORIGINS = ["terminal", "internal", "no primer"]
BUCKETS = ["0", "1-2", "3-5", "6-20", ">20", "other reference"]


def _share(value: float | None) -> str:
    return pct(value, 1, 2) if value is not None else "–"


def _scores(f: dict) -> str:
    rows = []
    for end, label in (("left", "5′ (Uni12)"), ("right", "3′ (Uni13)")):
        sc = f["scores"].get(end)
        if not sc or not sc["called"]:
            continue
        rows.append([label, sc["called"], pct(sc["agree"], sc["called"], 2),
                     _share(sc["internal_precision"]),
                     _share(sc["internal_recall"]),
                     _share(sc["terminal_precision"]),
                     _share(sc["terminal_recall"])])
    return table(["end", "ends called by both", "same origin",
                  "internal: confirmed", "internal: found by flucut",
                  "terminal: confirmed", "terminal: found by flucut"], rows)


def _origin(f: dict) -> str:
    out = []
    agree = {"terminal": "terminal", "internal": "internal"}
    for end, label in (("left", "5′ end (Uni12 primer)"),
                       ("right", "3′ end (Uni13 primer)")):
        counts = f["origin"].get(end)
        if not counts:
            continue
        out.append(f"<h3>{label}</h3>")
        out.append(cross_table("flucut \\ mapping", counts,
                               ["terminal", "internal", "unknown",
                                "no primer"], ORIGINS, agree))
    return "".join(out)


def _positions(f: dict) -> str:
    rows = []
    for end in ("left", "right"):
        counts = f["position"].get(end, {})
        total = sum(counts.values())
        if total:
            rows.append([end] + [pct(counts.get(b, 0), total, 2)
                                 for b in BUCKETS] + [total])
    return table(["end"] + [b + (" bases" if b[0].isdigit() or b[0] == ">"
                                 else "") for b in BUCKETS] + ["ends"],
                 rows)


def _trimming(f: dict) -> str:
    out = []
    labels = {"internal": "Internally primed ends",
              "terminal": "Terminal ends"}
    for origin in ("internal", "terminal"):
        counts = f["trim"].get(origin)
        if not counts:
            continue
        total = sum(counts.values())
        rows = [[k, v, pct(v, total, 2)] for k, v in
                sorted(counts.items(), key=lambda kv: -kv[1])]
        out.append(f"<h3>{labels[origin]}</h3>")
        out.append(table(["what flucut removed", "ends", "share"], rows))
    left = f["primer_bases_left"]
    if left:
        total = sum(left.values())
        mean = sum(int(k) * v for k, v in left.items()) / total
        out.append(para(
            f"Internally primed ends at which primer bases were left: "
            f"{num(total)}, with {mean:.1f} primer bases left on "
            f"average.", "note"))
    return "".join(out)


def _steps(s: dict) -> str:
    said = (s["flucut"]["summary"] or {}).get("reads")
    if not said:
        return ""
    classes, units = s["reads"]["classes"], s["amplicons"]
    primers = units["primers"]
    off = classes.get("offtarget", 0) + classes.get("dimer", 0)
    ends = s["ends"]
    said_ends = s["flucut"]["summary"].get("ends") or {}

    def measured(kind: str, origin: str) -> int:
        return sum(ends.get(side, {}).get(kind, {}).get(origin, 0)
                   for side in ("left", "right"))

    rows = [
        ["reads", s["reads"]["total"], said["in"]],
        ["amplicons with both primers", primers.get("2", 0) + off,
         said["both_primers"]],
        ["amplicons with one primer",
         primers.get("1", 0) + classes.get("one_primer", 0),
         said["one_primer"]],
        ["amplicons with no primer",
         primers.get("0", 0) + classes.get("no_primer", 0),
         said["no_primer"]],
        ["primer dimers", classes.get("dimer", 0), said["dimers"]],
        ["reads with several amplicons",
         sum(v for k, v in s["reads"]["amplicons_per_read"].items()
             if k != "1"), said["split"]]]
    for end, kind in (("left", "12"), ("right", "13")):
        for origin in ("terminal", "internal"):
            rows.append([f"Uni{kind} ends, {origin}",
                         measured(kind, origin),
                         said_ends.get(end, {}).get(origin)])
    for row in rows:
        a, b = row[1], row[2]
        row.append(pct(b, a) if a and b is not None else "–")
    return "<h3>Counts of the run</h3>" + para(
        "The counts of the summary of flucut next to those measured. "
        "flucut counts the unknown ends apart; the mapping places them, "
        "so its terminal and internal counts are larger.", "note") + \
        table(["", "measured", "flucut", "flucut / measured"], rows)


def _split(f: dict) -> str:
    counts: dict = {}
    for key, n in f["split"].items():
        said, truth = key.split(",")
        counts.setdefault(said, {})[truth] = n
    label = lambda k: "4 or more" if k == "4" else k
    rows = {label(k): {label(t): n for t, n in v.items()}
            for k, v in counts.items()}
    return cross_table("parts by flucut \\ amplicons by mapping", rows,
                       ["1", "2", "3", "4 or more"],
                       ["0", "1", "2", "3", "4 or more"])


def _categories(f: dict) -> str:
    cols = [k for k, _, _ in MEASURED_BANDS][::-1]
    counts = {said: {CATEGORY_LABELS.get(t, t): n for t, n in row.items()}
              for said, row in f["category"].items()}
    return cross_table("flucut \\ mapping", counts,
                       ["terminal+", "terminal-", "internal+", "internal-",
                        "no_primer", "deletion", "not_located"],
                       [CATEGORY_LABELS[c] for c in cols])


def flucut_sections(s: dict) -> str:
    f = s["flucut"]
    written = f["written"].get("not written", {})
    parts = [
        para("flucut classifies each read end from a probe of 32 bases "
             "after the primer; here the same reads are compared with "
             "their complete mapping. A read is compared when flucut and "
             "the mapping find the same number of amplicons in it."),
        _scores(f),
        para("<b>Confirmed</b>: of the ends flucut calls internal (or "
             "terminal), the share which the mapping places there. "
             "<b>Found by flucut</b>: of the ends the mapping places "
             "there, the share flucut calls so (the others it calls "
             "unknown, or the opposite).", "note"),
        "<h3>Origin of the ends, read by read</h3>", _origin(f),
        para("Other: reads which do not map, or in which flucut and the "
             "mapping find different numbers of amplicons.", "note"),
        "<h3>Position of the ends</h3>",
        para("Distance between the position flucut reports for an end "
             "(the tags <code>ml</code>, <code>mr</code>) and that of "
             "the mapping."), _positions(f)]
    if f["reads_are_flucut_output"]:
        parts.append(para(
            "The reads analysed are those written by flucut: their "
            "primers are already removed. An end at which flucut found "
            "a primer is therefore judged by its position alone (at a "
            "terminus or inside the segment), and the bases removed "
            "cannot be checked. Give the untrimmed reads to check the "
            "primers and the trimming.",
            "note"))
    else:
        parts += [
            "<h3>Bases removed</h3>",
            para("For the reads with one amplicon: the bases flucut "
                 "removed at an end with a primer, compared with the "
                 "bases before the template found by the mapping. At an "
                 "internally primed end the whole primer must go; at a "
                 "terminal end flucut keeps the conserved terminus when "
                 "it is sure of the terminus, and removes the whole "
                 "primer otherwise."), _trimming(f)]
    parts += [_steps(s),
              "<h3>Amplicons in a read</h3>", _split(f),
              "<h3>Category in the coverage plot</h3>",
              para("The category flucut gives a read in its coverage "
                   "plot (rebuilt from the tags), by the category "
                   "measured."), _categories(f)]
    only = f["sites_only_flucut"]
    if only:
        rows = [[x["reference"], x["end"], x["position"],
                 f"{x['from']}–{x['to']}", x["reads"]] for x in only[:15]]
        parts += ["<h3>Sites of flucut which the mapping does not "
                  "show</h3>",
                  table(["reference", "end", "position", "range", "reads"],
                        rows, left=2)]
    if written:
        rows = [[k, v] for k, v in sorted(written.items(),
                                          key=lambda kv: -kv[1])]
        parts += ["<h3>Reads not written by flucut</h3>",
                  table(["what the read is", "reads"], rows)]
    return section("flucut", "flucut compared with the mapping", *parts)


def irma_section(s: dict) -> str:
    parts = [para("What IRMA did with the reads: the reads at each of "
                  "its steps, and which kinds of reads it assembled (the "
                  "reads in its BAM files).")]
    for run in s["irma"]:
        parts.append(f"<h3>{escape(run['label'])}"
                     f"{' (' + escape(run['module']) + ')' if run['module'] else ''}"
                     f"</h3>")
        rows = [[k, v] for k, v in run["read_counts"]]
        kinds = []
        for kind, counts in sorted(
                run["reads_by_kind"].items(),
                key=lambda kv: -sum(kv[1].values())):
            used, unused = counts.get("True", 0), counts.get("False", 0)
            kinds.append([CATEGORY_LABELS.get(kind, kind), used, unused,
                          pct(used, used + unused)])
        same = run["same_segment"]
        note = ""
        if same:
            total = sum(same.values())
            note = para(f"Of the reads assembled by IRMA which map, "
                        f"{escape(pct(same.get('True', 0), total, 2))} are "
                        f"assigned by IRMA to the segment they map to.",
                        "note")
        parts.append('<div class="pair"><div>'
                     + table(["step of IRMA", "reads"], rows) + "</div><div>"
                     + table(["what the read is", "assembled",
                              "not assembled", "assembled"], kinds)
                     + note + "</div></div>")
    return section("irma", "IRMA", *parts)
