"""The HTML report of a run: what the mapping shows about the reads."""

import os
from html import escape

from .html import STYLE, num, para, pct, section, table, tile
from .segment_coverage import (AGREEMENT_COLUMNS, agreement_row,
                               coverage_pair)
from .plots import (MEASURED_BANDS, PREDICTED_BANDS, bars, legend,
                    length_histogram)
from .report_compare import flucut_sections, irma_section

CLASS_LABELS = [
    ("influenza", "influenza amplicons (mapped)", "#0b9e73"),
    ("offtarget", "other amplicons (both primers, not mapped)", "#e8a013"),
    ("one_primer", "not mapped, one primer", "#c9a45c"),
    ("no_primer", "not mapped, no primer", "#9e9e9e"),
    ("dimer", "primer dimers", "#d7263d"),
]


def _head(s: dict) -> str:
    inputs = s["inputs"]
    primers = ", ".join(f"{p['name']} ({p['sequence']})"
                        for p in s["settings"]["primers"])
    rows = [["reads", ", ".join(inputs["reads"])],
            ["IRMA results", ", ".join(inputs["irma"]) or "none"],
            ["flucut reads", inputs["flucut_reads"] or (
                "the reads themselves" if s["flucut"] and
                s["flucut"]["reads_are_flucut_output"] else "none")],
            ["flucut summary", inputs["flucut_summary"] or "none"],
            ["flucut report", inputs["flucut_report"] or "none"],
            ["references", inputs["references"]],
            ["primers", primers],
            ["settings", f"primer errors up to "
             f"{s['settings']['error_rate']:.0%} of the primer; terminal: "
             f"within {s['settings']['terminal_window']} bases of the "
             f"terminus; deletion: at least "
             f"{s['settings']['min_deletion']} bases"]]
    cells = "".join(f"<tr><td>{escape(k)}</td><td>{escape(v)}</td></tr>"
                    for k, v in rows)
    return (f"<h1>flucut post-hoc: {escape(s['sample'])}</h1>"
            f'<div class="sub">measured from the base-by-base mapping of '
            f"every read · flucut-posthoc {s['version']} · "
            f"{escape(s['updated'])}</div>"
            f'<table class="params">{cells}</table>')


def _tiles(s: dict) -> str:
    reads, classes = s["reads"]["total"], s["reads"]["classes"]
    amplicons = s["amplicons"]
    complete = sum(c["amplicons"] for c in s["complete_amplicons"])
    left = s["ends"].get("left", {})
    primed = sum(sum(v.values()) for k, v in left.items() if k != "none")
    internal = sum(v.get("internal", 0) for k, v in left.items()
                   if k != "none")
    items = [tile("reads", num(reads)),
             tile("map to influenza",
                  pct(classes.get("influenza", 0), reads)),
             tile("amplicons of another template",
                  pct(classes.get("offtarget", 0), reads)),
             tile("influenza amplicons", num(amplicons["total"])),
             tile("with both primers", pct(complete, amplicons["total"])),
             tile("5′ ends internally primed", pct(internal, primed)),
             tile("with an internal deletion",
                  pct(amplicons["categories"].get("deletion", 0),
                      amplicons["total"]))]
    if s["flucut"]:
        scores = s["flucut"]["scores"].get("left", {})
        if scores.get("called"):
            items.append(tile("flucut 5′ origin calls confirmed",
                              pct(scores["agree"], scores["called"])))
    return f'<div class="tiles">{"".join(items)}</div>'


def _reads(s: dict) -> str:
    reads = s["reads"]
    items = [(label, reads["classes"].get(key, 0), colour)
             for key, label, colour in CLASS_LABELS]
    per_read = reads["amplicons_per_read"]
    several = sum(v for k, v in per_read.items() if k != "1")
    text = (f"{num(several)} of the influenza reads "
            f"({escape(pct(several, sum(per_read.values())))}) contain more "
            f"than "
            f"one influenza amplicon; "
            f"{num(reads['joined_not_influenza'])} are joined to an "
            f"amplicon which is not influenza.")
    junctions = "".join(f" {escape(k)}: {num(v)}."
                        for k, v in sorted(reads["junctions"].items()))
    return section(
        "reads", "What the reads are",
        para("Every read is mapped to the segment references. A read "
             "which does not map is looked at for the primers at its two "
             "ends: with both and nothing between them it is a primer "
             "dimer, with both and a template it is an amplicon of "
             "something which is not influenza."),
        bars(items, reads["total"]),
        para(text + (" Junctions between consecutive amplicons of a "
                     "read:" + junctions if junctions else ""), "note"))


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
                row.append(("no consensus", "bad"))
                continue
            covered = info["reference_covered"] / r["length"]
            covered_text = pct(info["reference_covered"], r["length"], 0)
            text = f"{info['consensus_length']:,} ({covered_text})"
            if info["deletions"]:
                text += "; deleted " + ", ".join(
                    f"{a}–{b}" for a, b in info["deletions"])
            if info["primers_inside"]:
                text += f"; {info['primers_inside']} primers inside"
            row.append((text, "bad" if covered < 0.9
                        or info["primers_inside"] else ""))
        rows.append(row)
    headers = ["segment", "reference", "length", "terminus bases 12 / 13",
               "bases of the sample", "amplicons", "read identity"] + \
        [f"consensus of {label}" for label in labels]
    return section(
        "references", "References",
        para("For each segment, the reference of the set to which the "
             "reads map best. Where most reads of the sample have "
             "another base than the reference, the reference is given "
             "that base before the mapping (bases of the sample); "
             "insertions and deletions are not applied, so the positions "
             "in this report are those of the references. The read "
             "identity is that of the reads to the reference so "
             "adapted. A database record often lacks "
             "the conserved termini (0 terminus bases): position 1 is "
             "then the first base after the Uni12 primer. For each IRMA "
             "run: the length of its consensus, the share of the "
             "reference it covers, its deletions and the primers found "
             "inside it (a consensus built from amplicons joined end to "
             "end)."),
        table(headers, rows, left=2))


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
    offsets = s["terminal_offsets"]
    exact = {side: escape(pct(sum(v for k, v in offsets.get(side, {}).items()
                                  if abs(int(k)) <= 2),
                              sum(offsets.get(side, {}).values())))
             for side in ("left", "right")}
    return section(
        "ends", "Ends of the amplicons",
        para("Each end of an amplicon is where its template begins or "
             "ends on the reference. With a primer before it, it is "
             "<b>terminal</b> (at the terminus of the segment: a normal "
             "amplicon) or <b>internal</b> (the primer annealed inside "
             "the segment). Without a primer it is the end of a fragment "
             "(a broken molecule, an incomplete read)."),
        table(["end", "primer", "terminal", "internal", "total"], rows,
              left=2),
        para(f"Terminal ends exactly at the terminus (within 2 bases): "
             f"left {exact['left']}, right {exact['right']}.", "note"),
        "<h3>How well the primers were found</h3>",
        para("Errors in the primers found (share of the primers with 0 "
             "to 4 errors), and primers of which the read lacks the "
             "first bases. A control is the search, at the same places, of "
             "the bases of each primer in another order: its hits are "
             "chance matches, and tell how many of the primers found can "
             "be expected to be chance matches as well. Three controls "
             "are searched; the rate is that of the middle one (one "
             "control alone can resemble an adapter or a barcode)."),
        table(["primer", "0", "1", "2", "3", "4", "truncated", "found"],
              error_rows),
        table(["control search", "searched", "found by each control",
               "rate"],
              control_rows),
        "<h3>What precedes the primers</h3>",
        para("The most frequent bases of the reads just before each "
             "primer, as shares of the primers found: the tail of the "
             "primer used in the laboratory, or the barcode or adapter "
             "next to an untailed primer."),
        table(["primer", "bases before it", "share"], _before(search)))


def _before(search: dict) -> list[list]:
    rows = []
    for name, tails in search.get("before_primer", {}).items():
        total = tails[-1][1]
        rows += [[name, (seq, "seq"), pct(n, total)]
                 for seq, n in tails[:-1][:3]]
    return rows


def _sites(s: dict) -> str:
    sites = s["internal_sites"]
    if not sites:
        return section("sites", "Internal priming sites",
                       para("No internally primed end was found."))
    with_flucut = s["flucut"] is not None
    rows = []
    for site in sites[:40]:
        row = [site["gene"], site["reference"], site["side"],
               "Uni" + site["primer"], site["position"],
               f"{site['from']}–{site['to']}", site["reads"],
               (site["template"], "seq"),
               site["mismatches"] if site["mismatches"] >= 0 else None]
        if with_flucut:
            row.append(site["flucut_reads"])
        rows.append(row)
    headers = ["segment", "reference", "end", "primer", "position",
               "range", "reads", "template under the primer",
               "mismatches"] + (["flucut reads"] if with_flucut else [])
    return section(
        "sites", "Internal priming sites",
        para("The ends with a primer which are inside a segment, grouped "
             "into sites (positions at most 5 bases apart), by decreasing "
             "number of reads. The position is that of the first base "
             "after the primer (left end) or of the last base before it "
             "(right end), as in flucut. The template under the primer is "
             "the sequence of the reference where the primer annealed, "
             "with its mismatches to the conserved terminus of the "
             "primer: few mismatches explain why the primer anneals "
             f"there. {num(s['internal_sites_total'])} sites in all; the "
             "largest are shown."),
        table(headers, rows, left=4))


def _complete(s: dict) -> str:
    rows, plots = [], []
    with_flucut = any(c["flucut"] for c in s["complete_amplicons"])
    bin_size = s["settings"]["length_bin"]
    for c in s["complete_amplicons"]:
        f = c["flucut"] or {}
        row = [c["gene"], c["reference_length"], c["amplicons"],
               c["main_length"], c["short"], c["short_main_length"]]
        if with_flucut:
            row += [f.get("amplicons"), f.get("main_length"),
                    f.get("short"), f.get("short_main_length")]
        rows.append(row)
        if c["amplicons"] or f.get("amplicons"):
            plots.append(
                f'<div class="plot"><b>{escape(c["gene"])}</b>'
                + length_histogram(c["histogram"], f.get("histogram"),
                                   bin_size, c["reference_length"])
                + "</div>")
    headers = ["segment", "reference length", "amplicons", "usual length",
               "short", "usual short length"]
    if with_flucut:
        headers += ["flucut: amplicons", "usual length", "short",
                    "usual short length"]
    return section(
        "complete", "Complete amplicons",
        para("The amplicons with the Uni12 primer at their left end and "
             "the Uni13 primer at their right end, by segment: the length "
             "of the template between the primers, and those shorter "
             "than 90% of the reference (internal priming, deletions). "
             "Bars: measured; red line: flucut (from its summary); dashed "
             "line: the length of the reference."),
        table(headers, rows), f'<div class="plots">{"".join(plots)}</div>')


def _coverage(s: dict) -> str:
    labels = [run["label"] for run in s["irma"]]
    with_flucut = any(c["predicted"] for c in s["coverage"].values())
    parts, rows = [], []
    for ref in s["references"]:
        pair = coverage_pair(s, ref, ref["gene"])
        if not pair:
            continue
        parts.append(pair)
        agreement = s["coverage"][ref["name"]]["agreement"]
        if agreement:
            rows.append(agreement_row(ref["gene"], agreement))
    text = ("The depth at each position of each segment, from the "
            "mappings (the bases of a read aligned there; a deletion of "
            "a read is not covered), stacked by what the amplicon is: "
            "where its 5′ (left) end is and whether a primer is there. "
            "The dashed red line is the largest internal priming site. ")
    if labels:
        text += ("The lines are the depths reported by IRMA, at the "
                 "positions of its consensus on the reference. ")
    out = [para(text), legend(MEASURED_BANDS, labels)]
    if s["predicted_coverage_from"] == "report":
        out.append(para(
            "At the right, on the same scale, the coverage which flucut "
            "predicts: the plot of its own report, at the positions of "
            "the reference used here.", "note"))
    elif with_flucut:
        out.append(para(
            "At the right, on the same scale, the coverage which flucut "
            "predicts, rebuilt from the tags of its reads as flucut "
            "builds its plot. The tags do not say where flucut placed a "
            "read by a probe from inside it (such reads are missing "
            "here) nor where it located the junction of a deletion (the "
            "read is drawn as two halves). Give the report of flucut "
            "(<code>--flucut-report</code>) to compare with the plot it "
            "really shows.", "note"))
        out.append(legend(PREDICTED_BANDS))
    out.extend(parts)
    if rows:
        out.append("<h3>Predicted and measured depth</h3>")
        out.append(table(["segment"] + AGREEMENT_COLUMNS, rows))
        out.append(para(
            "The depths are means over the segment. The difference per "
            "position is the mean of the absolute differences between "
            "the two depths, as a share of the measured mean depth. The "
            "correlation along the segment says whether the two "
            "profiles have the same shape; it is low for a flat profile "
            "even when the depths agree.", "note"))
    return section("coverage", "Coverage by read origin", *out)


def _other(s: dict) -> str:
    genes = {r["name"]: r["gene"] for r in s["references"]}
    out = []
    if s["deletions"]:
        rows = [[genes[d["reference"]], f"{d['from'] + 1}–{d['to']}",
                 d["to"] - d["from"], d["reads"]]
                for d in s["deletions"][:20]]
        out.append(section(
            "deletions", "Internal deletions",
            para("The most frequent deletions of the amplicons (bases of "
                 "the reference missing in the read; positions in bins "
                 "of 10): genomes with a large internal deletion, such "
                 "as defective interfering genomes."),
            table(["segment", "deleted", "bases", "amplicons"], rows)))
    if s["fragment_ends"]:
        rows = [[genes[d["reference"]], d["side"], d["position"],
                 d["reads"]] for d in s["fragment_ends"][:15]]
        out.append(section(
            "fragments", "Ends without a primer",
            para("The positions inside the segments at which most "
                 "amplicons end without a primer. Ends of broken "
                 "molecules are spread along a segment; many at one "
                 "position are a step of the coverage which no primer "
                 "explains."),
            table(["segment", "end", "position", "amplicons"], rows)))
    if s["offtarget"]:
        rows = [[i, pct(c["share_of_sample"], 1), c["template_length"],
                 (c["representative"][:90] + ("…" if len(
                     c["representative"]) > 90 else ""), "seq")]
                for i, c in enumerate(s["offtarget"], 1)]
        out.append(section(
            "offtarget", "Amplicons which are not influenza",
            para("The reads which do not map but begin with the Uni12 "
                 "primer, grouped by their first bases after it: share "
                 "of these reads, usual length of the template, and a "
                 "read of the group (complete in "
                 "<code>offtarget_amplicons.fasta</code>, for a database "
                 "search)."),
            table(["group", "share", "template length", "sequence"],
                  rows)))
    return "".join(out)


def write_report(s: dict, path: str) -> None:
    """Write the HTML report of the summary `s`."""
    nav = [("reads", "Reads"), ("references", "References"),
           ("ends", "Ends"), ("sites", "Internal priming"),
           ("complete", "Complete amplicons"), ("coverage", "Coverage")]
    if s["flucut"]:
        nav.append(("flucut", "flucut compared"))
    if s["irma"]:
        nav.append(("irma", "IRMA"))
    links = "".join(f'<a href="#{a}">{t}</a>' for a, t in nav)
    body = [_head(s), f"<nav>{links}</nav>", _tiles(s), _reads(s),
            _references(s), _ends(s), _sites(s), _complete(s),
            _coverage(s)]
    if s["flucut"]:
        body.append(flucut_sections(s))
    if s["irma"]:
        body.append(irma_section(s))
    body.append(_other(s))
    page = (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,'
            f'initial-scale=1"><title>flucut post-hoc: '
            f'{escape(s["sample"])}</title><style>{STYLE}</style></head>'
            f'<body><main>{"".join(body)}</main></body></html>')
    with open(path + ".tmp", "w") as handle:
        handle.write(page)
    os.replace(path + ".tmp", path)
