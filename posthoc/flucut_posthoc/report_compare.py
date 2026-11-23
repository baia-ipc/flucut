"""The parts of the report on flucut read by read and on IRMA: the
tables behind the checks, each shown on request, and what IRMA did with
each kind of read."""

from .agreement import origin_scores, origin_table
from .html import card, cross_table, fold, more, note, num, pct, table
from .origins import ORIGIN_LABELS, ORIGINS
from .report_details import CATEGORY_LABELS

READ_LABELS = dict(CATEGORY_LABELS)
READ_LABELS.update({
    "several": "several amplicons in the read",
    "dimer": "primer dimer", "offtarget": "other amplicon, not mapped",
    "one_primer": "not mapped, one primer",
    "no_primer": "not mapped, no primer"})
FLUCUT_LABELS = {
    "terminal+": "full template, plus strand",
    "terminal-": "full template, minus strand",
    "internal+": "internally primed, plus strand",
    "internal-": "internally primed, minus strand",
    "no_primer": "other: no primer at the 5′ end",
    "deletion": "internal deletion", "not_located": "other: not located"}
# the kinds of reads, in few words
SHORT_LABELS = dict(READ_LABELS)
SHORT_LABELS.update({
    "terminal+": "full template, plus strand",
    "terminal-": "full template, minus strand",
    "internal+": "internally primed, plus strand",
    "internal-": "internally primed, minus strand",
    "internal3": "internally primed at the 3′ end",
    "noprimer_terminal": "other, from the terminus",
    "noprimer_fragment": "other, from inside"})
END_ORIGINS = ["terminal", "internal", "no primer"]
BUCKETS = ["0", "1-2", "3-5", "6-20", ">20", "other reference"]
END_TITLES = (("left", "5′ end (Uni12 primer)"),
              ("right", "3′ end (Uni13 primer)"))


def _share(value: float | None) -> str:
    return pct(value, 1, 2) if value is not None else "–"


def _scores(f: dict) -> str:
    rows = []
    for end, sc in origin_scores(f["origin"]).items():
        label = dict(END_TITLES).get(end)
        if not label or not sc["called"]:
            continue
        rows.append([label, sc["called"], pct(sc["agree"], sc["called"], 2),
                     _share(sc["internal_precision"]),
                     _share(sc["internal_recall"]),
                     _share(sc["terminal_precision"]),
                     _share(sc["terminal_recall"])])
    return table(["end", "ends called by both", "same origin",
                  "confirmed", "found by flucut", "confirmed",
                  "found by flucut"], rows,
                 groups=[("", 3), ("internal", 2), ("terminal", 2)])


def _origin(f: dict) -> str:
    out = [_scores(f)]
    agree = {"terminal": "terminal", "internal": "internal"}
    for end, label in END_TITLES:
        counts = f["origin"].get(end)
        if counts:
            out.append(f"<h4>{label}</h4>" + cross_table(
                "flucut \\ mapping", counts,
                ["terminal", "internal", "unknown", "no primer"],
                END_ORIGINS, agree))
    out.append(note("Other: reads which do not map, or in which flucut "
                    "and the mapping find different numbers of "
                    "amplicons."))
    return fold("Origin of the ends, read by read", *out,
                hint="the call of flucut for each end against the mapping")


def _positions(f: dict) -> str:
    rows = []
    for end in ("left", "right"):
        counts = f["position"].get(end, {})
        total = sum(counts.values())
        if total:
            rows.append([end] + [pct(counts.get(b, 0), total, 2)
                                 for b in BUCKETS] + [total])
    if not rows:
        return ""
    return fold(
        "Position of the ends",
        table(["end"] + BUCKETS + ["ends"], rows,
              groups=[("", 1), ("bases between flucut and the mapping, "
                                "of the ends located by both", 6),
                      ("", 1)]),
        hint="distance between the tags ml, mr and the mapping")


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
        out.append(f"<h4>{labels[origin]}</h4>" + table(
            ["what flucut removed", "ends", "of these ends"], rows))
    left = f["primer_bases_left"]
    if left:
        total = sum(left.values())
        mean = sum(int(k) * v for k, v in left.items()) / total
        out.append(note(f"Internally primed ends with primer bases left: "
                        f"{num(total)}, {mean:.1f} bases on average."))
    return fold("Bases removed", *out,
                hint="at the ends with a primer, reads with one amplicon")


def _counts(s: dict) -> str:
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
    return fold(
        "Counts of the run",
        table(["", "measured", "flucut", "flucut / measured"], rows),
        note("flucut counts the ends it cannot place apart (unknown); "
             "the mapping places them, so its counts are larger."),
        hint="the counts of the summary of flucut next to those measured")


def _split(f: dict) -> str:
    counts: dict = {}
    for key, n in f["split"].items():
        said, truth = key.split(",")
        counts.setdefault(said, {})[truth] = n
    label = lambda k: "4 or more" if k == "4" else k
    rows = {label(k): {label(t): n for t, n in v.items()}
            for k, v in counts.items()}
    if not rows:
        return ""
    return fold("Amplicons in a read", cross_table(
        "parts by flucut \\ amplicons by mapping", rows,
        ["1", "2", "3", "4 or more"], ["0", "1", "2", "3", "4 or more"]),
        hint="reads split by flucut against the amplicons mapped")


def _categories(f: dict) -> str:
    if not f["category"]:
        return ""
    by_origin = origin_table(f["category"])
    labels = [label for _, label in ORIGINS[1:] + ORIGINS[:1]]
    origins = {ORIGIN_LABELS[a]: {ORIGIN_LABELS[b]: n for b, n in
                                  row.items()}
               for a, row in by_origin.items()}
    cols = list(CATEGORY_LABELS)
    counts = {FLUCUT_LABELS.get(said, said):
              {READ_LABELS.get(t, t): n for t, n in row.items()}
              for said, row in f["category"].items()}
    return fold(
        "Read origin, read by read",
        cross_table("flucut \\ mapping", origins, labels, labels,
                    {label: label for label in labels}, rest=False),
        "<h4>By category</h4>",
        cross_table("flucut \\ mapping", counts,
                    list(FLUCUT_LABELS.values()),
                    [READ_LABELS[c] for c in cols]),
        note("Other: reads which do not map, or split differently."),
        hint="the origin flucut gives a read against the one measured")


def _rest(f: dict) -> str:
    out = ""
    only = f["sites_only_flucut"]
    if only:
        rows = [[x["reference"], x["end"], x["position"],
                 f"{x['from']}–{x['to']}", x["reads"]] for x in only[:15]]
        out += fold("Sites of flucut which the mapping does not show",
                    table(["reference", "end", "position", "range",
                           "reads"], rows, left=2),
                    hint=f"{len(only)} sites")
    written = f["written"].get("not written", {})
    if written:
        rows = [[READ_LABELS.get(k.replace("unmapped: ", ""), k), v]
                for k, v in sorted(written.items(), key=lambda kv: -kv[1])]
        out += fold("Reads not written by flucut",
                    table(["what the read is", "reads"], rows),
                    hint=f"{sum(written.values()):,} reads")
    return out


def flucut_details(s: dict) -> str:
    """The tables of the comparison with flucut, each shown on
    request."""
    f = s["flucut"]
    out = ""
    if f["origin"]:
        out += _origin(f) + _positions(f)
        if not f["reads_are_flucut_output"]:
            out += _trimming(f)
    return out + _counts(s) + _split(f) + _categories(f) + _rest(f)


def _kinds(s: dict) -> list[str]:
    """The kinds of reads which an IRMA run has, most frequent first."""
    totals: dict[str, int] = {}
    for run in s["irma"]:
        for kind, counts in run["reads_by_kind"].items():
            totals[kind] = totals.get(kind, 0) + sum(counts.values())
    return sorted(totals, key=lambda kind: -totals[kind])


def irma_card(s: dict) -> str:
    """For each kind of read, the share which each IRMA run
    assembled."""
    if not s["irma"]:
        return ""
    rows = []
    for kind in _kinds(s):
        row = [SHORT_LABELS.get(kind, kind)]
        for run in s["irma"]:
            counts = run["reads_by_kind"].get(kind, {})
            used, unused = counts.get("True", 0), counts.get("False", 0)
            row += [pct(used, used + unused), used + unused]
        rows.append(row)
    headers = ["what the read is"] + ["assembled", "of reads"] * \
        len(s["irma"])
    groups = [("", 1)] + [(run["label"], 2) for run in s["irma"]]
    lines = []
    for run in s["irma"]:
        same = run["same_segment"]
        if same:
            lines.append(f"{run['label']}: "
                         f"{pct(same.get('True', 0), sum(same.values()), 2)}"
                         f" of the reads assembled are assigned to the "
                         f"segment they map to")
    return card("Reads assembled by IRMA",
                table(headers, rows[:6], groups=groups),
                more(f"{len(rows) - 6} more kinds of reads",
                     table(headers, rows[6:], groups=groups))
                if len(rows) > 6 else "",
                note(". ".join(lines) + ".") if lines else "")


def irma_details(s: dict) -> str:
    out = ""
    for run in s["irma"]:
        kinds = []
        for kind, counts in sorted(run["reads_by_kind"].items(),
                                   key=lambda kv: -sum(kv[1].values())):
            used, unused = counts.get("True", 0), counts.get("False", 0)
            kinds.append([SHORT_LABELS.get(kind, kind), used, unused,
                          pct(used, used + unused)])
        module = f" ({run['module']})" if run["module"] else ""
        out += fold(
            f"IRMA run {run['label']}{module}",
            table(["step of IRMA", "reads"],
                  [[k, v] for k, v in run["read_counts"]]),
            table(["what the read is", "assembled", "not assembled",
                   "assembled, of these reads"], kinds),
            hint="reads at each step, reads assembled by kind")
    return out
