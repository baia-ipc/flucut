"""The tables of the report with every number measured, each shown on
request, and the terms of the report."""

from .agreement import LIMITS, MIN_CASES, depth_limits
from .html import fold, note, num, pct, table, terms
from .origins import MEASURED_ORIGIN, ORIGIN_LABELS

CATEGORY_LABELS = {
    "terminal+": "primer at the 5′ end, at the terminus, plus strand",
    "terminal-": "primer at the 5′ end, at the terminus, minus strand",
    "internal+": "primer at the 5′ end, inside the segment, plus strand",
    "internal-": "primer at the 5′ end, inside the segment, minus strand",
    "internal3": "5′ end at the terminus, primer at the 3′ end inside",
    "noprimer_terminal": "no primer at the 5′ end, starts at the terminus",
    "noprimer_fragment": "no primer at the 5′ end, starts inside",
    "deletion": "internal deletion"}
SITES_LISTED = 40


def _references(s: dict) -> str:
    labels = [run["label"] for run in s["irma"]]
    rows = []
    for r in s["references"]:
        row = [r["gene"], r["name"], r["length"],
               f"{r['terminus_12_bases']} / {r['terminus_13_bases']}",
               r["adapted_bases"], r["amplicons"],
               pct(r["identity"], 1, 2) if r["identity"] else "–"]
        for label in labels:
            info = r["irma"].get(label)
            if info is None:
                row.append("! no consensus")
                continue
            text = (f"{info['consensus_length']:,} "
                    f"({pct(info['reference_covered'], r['length'], 0)})")
            if info["deletions"]:
                text += "; deleted " + ", ".join(
                    f"{a}–{b}" for a, b in info["deletions"])
            if info["primers_inside"]:
                text += f"; {info['primers_inside']} primers inside"
            flagged = info["reference_covered"] < 0.9 * r["length"] \
                or info["primers_inside"]
            row.append(("! " if flagged else "") + text)
        rows.append(row)
    headers = ["segment", "reference", "length", "terminus bases 12 / 13",
               "bases of the sample", "amplicons", "read identity"] + \
        [f"consensus of {label}" for label in labels]
    return fold(
        "References", table(headers, rows, left=2),
        note("Bases of the sample: bases of the reference replaced by "
             "those of most reads before the mapping. Consensus: its "
             "length and the share of the reference it covers; ! marks "
             "one which covers under 90% or has primers inside."),
        hint="the reference of each segment and the consensus of IRMA")


def _ends(s: dict) -> str:
    rows = []
    for side in ("left", "right"):
        for kind, label in (("12", "Uni12 primer"), ("13", "Uni13 primer"),
                            ("none", "no primer")):
            counts = s["ends"].get(side, {}).get(kind, {})
            total = sum(counts.values())
            if total:
                rows.append([f"{side} ({'start' if side == 'left' else 'end'}"
                             f" of the segment)", label,
                             counts.get("terminal", 0),
                             counts.get("internal", 0), total])
    offsets = s["terminal_offsets"]
    exact = {side: pct(sum(v for k, v in offsets.get(side, {}).items()
                           if abs(int(k)) <= 2),
                       sum(offsets.get(side, {}).values()))
             for side in ("left", "right")}
    return fold(
        "Ends of the amplicons",
        table(["end", "primer", "terminal", "internal", "total"], rows,
              left=2),
        note(f"Terminal ends within 2 bases of the terminus: left "
             f"{exact['left']}, right {exact['right']}."),
        hint="terminal or internal, by primer")


def _primer_search(s: dict) -> str:
    search = s["primer_search"]
    error_rows = []
    for name, errors in search["errors"].items():
        total = sum(errors.values())
        error_rows.append([name] + [pct(errors.get(str(e), 0), total)
                                    for e in range(5)] +
                          [search["truncated"].get(name, 0), total])
    control_rows = [[context, c["searched"],
                     ", ".join(f"{n:,}" for n in c["found_each"]),
                     pct(c["found"], c["searched"], 2)]
                    for context, c in search["control"].items()]
    before = []
    for name, tails in search.get("before_primer", {}).items():
        total = tails[-1][1]
        before += [[name, (seq, "seq"), pct(n, total)]
                   for seq, n in tails[:-1][:3]]
    return fold(
        "Primer search",
        table(["primer", "0", "1", "2", "3", "4", "truncated", "found"],
              error_rows, groups=[("", 1), ("errors, of the primers found",
                                            5), ("", 2)]),
        table(["control search", "searched", "found by each control",
               "rate"], control_rows),
        note("Controls: the bases of each primer in another order, "
             "searched at the same places; their hits are chance "
             "matches. The rate is that of the middle one of three."),
        table(["primer", "bases before it", "of the primers found"],
              before),
        hint="errors, chance matches, bases before the primers")


def _sites(s: dict) -> str:
    sites = s["internal_sites"]
    if not sites:
        return ""
    with_flucut = bool((s["flucut"] or {}).get("origin"))
    rows = []
    for site in sites[:SITES_LISTED]:
        row = [site["gene"], site["reference"], site["side"],
               "Uni" + site["primer"], site["position"],
               f"{site['from']}–{site['to']}", site["reads"],
               (site["template"], "seq"),
               site["mismatches"] if site["mismatches"] >= 0 else None]
        rows.append(row + [site["flucut_reads"]] if with_flucut else row)
    headers = ["segment", "reference", "end", "primer", "position",
               "range", "reads", "template under the primer",
               "mismatches"] + (["flucut reads"] if with_flucut else [])
    return fold(
        "Internal priming sites", table(headers, rows, left=4),
        note("Ends with a primer inside a segment, at most 5 bases apart. "
             "Template: the reference where the primer annealed; "
             "mismatches: with the conserved terminus of the primer."),
        hint=f"the largest of {s['internal_sites_total']:,}, with the "
             f"template under the primer")


def _complete(s: dict) -> str:
    with_flucut = any(c["flucut"] for c in s["complete_amplicons"])
    rows = []
    for c in s["complete_amplicons"]:
        f = c["flucut"] or {}
        row = [c["gene"], c["reference_length"], c["amplicons"],
               c["main_length"], c["short"], c["short_main_length"]]
        if with_flucut:
            row += [f.get("amplicons"), f.get("main_length"),
                    f.get("short"), f.get("short_main_length")]
        rows.append(row)
    headers = ["segment", "reference length"] + \
        ["amplicons", "usual length", "short", "usual short length"] \
        * (2 if with_flucut else 1)
    groups = [("", 2), ("measured", 4)] + \
        ([("flucut", 4)] if with_flucut else [])
    return fold("Complete amplicons",
                table(headers, rows, groups=groups),
                note("Amplicons with Uni12 at the left and Uni13 at the "
                     "right. Short: under 90% of the reference."),
                hint="number and usual length by segment")


def _categories(s: dict) -> str:
    counts = s["amplicons"]["categories"]
    total = s["amplicons"]["total"]
    rows = [[CATEGORY_LABELS[key], ORIGIN_LABELS[MEASURED_ORIGIN[key]],
             counts.get(key, 0), pct(counts.get(key, 0), total)]
            for key in CATEGORY_LABELS]
    return fold("Amplicons by category",
                table(["category", "origin", "amplicons", "of amplicons"],
                      rows, left=2),
                hint="the categories behind the four origins")


def _other(s: dict) -> str:
    genes = {r["name"]: r["gene"] for r in s["references"]}
    out = []
    if s["deletions"]:
        rows = [[genes[d["reference"]], f"{d['from'] + 1}–{d['to']}",
                 d["to"] - d["from"], d["reads"]]
                for d in s["deletions"][:20]]
        out.append(fold(
            "Internal deletions",
            table(["segment", "deleted", "bases", "amplicons"], rows),
            hint="the most frequent, positions in bins of 10"))
    if s["fragment_ends"]:
        rows = [[genes[d["reference"]], d["side"], d["position"],
                 d["reads"]] for d in s["fragment_ends"][:15]]
        out.append(fold(
            "Ends without a primer",
            table(["segment", "end", "position", "amplicons"], rows),
            hint="the positions where most amplicons end without one"))
    if s["offtarget"]:
        rows = [[i, pct(c["share_of_sample"], 1), c["template_length"],
                 (c["representative"][:90] + ("…" if len(
                     c["representative"]) > 90 else ""), "seq")]
                for i, c in enumerate(s["offtarget"], 1)]
        out.append(fold(
            "Amplicons which are not influenza",
            table(["group", "of these reads", "template length",
                   "sequence"], rows),
            note("Reads which do not map but begin with Uni12, grouped "
                 "by their first bases; complete sequences in "
                 "offtarget_amplicons.fasta."),
            hint="grouped, with a sequence to search in a database"))
    return "".join(out)


def measured_details(s: dict) -> str:
    return (_references(s) + _ends(s) + _primer_search(s) + _sites(s)
            + _complete(s) + _categories(s) + _other(s))


def terms_table(s: dict) -> str:
    min_depth, fraction = depth_limits(s)
    window = s["settings"]["terminal_window"]
    rows = [
        ("measured", "From the base-by-base mapping of every read to the "
         "reference of each segment."),
        ("predicted, flucut", "What flucut reports without mapping: the "
         "tags of its reads and its summary. Its depth is estimated from "
         "where it locates each read."),
        ("amplicon", "A part of a read which maps to a segment; a read "
         "can hold several, joined end to end."),
        ("left, right end", "The start and the end of the plus strand of "
         "the segment; Uni12 is expected at the left (5′), Uni13 at the "
         "right (3′)."),
        ("terminal, internal", f"An end within {window} bases of the "
         "terminus of the segment, or anywhere else. Internal with a "
         "primer: internally primed."),
        ("full template", "A primer at the 5′ end, at the terminus."),
        ("internally primed", "A primer inside the segment, at the 5′ "
         "end or, with the 5′ end at the terminus, at the 3′ end."),
        ("internal deletion", f"At least {s['settings']['min_deletion']} "
         "bases of the reference missing in the amplicon; only the "
         "bases present count as depth."),
        ("other", "No primer at the 5′ end; for flucut also reads it "
         "does not locate."),
        ("complete segment", f"At least {fraction:.0%} of its positions "
         f"at a depth of {min_depth}× or more; partial: at least half; "
         "else low depth."),
        ("confirmed, found", "Confirmed: of the ends flucut calls so, "
         "those the mapping places there. Found: of the ends the mapping "
         "places there, those flucut calls so."),
        ("difference per position", "Mean absolute difference between "
         "the depth of flucut and the depth measured, as a share of the "
         "mean depth measured."),
        ("✓ ! ✕ ·", "Within the limit of the measure, within its wider "
         "limit, outside both; · fewer than "
         f"{MIN_CASES} cases, or nothing to compare. A segment agrees "
         f"when its state is the same and its depth is within "
         f"{LIMITS['depth_ratio'][0]:.0%}."),
        ("<0.1%, >99.9%", "Not zero, not all: 0.0% and 100.0% are "
         "written only for exactly none and exactly all."),
    ]
    return terms(rows)
