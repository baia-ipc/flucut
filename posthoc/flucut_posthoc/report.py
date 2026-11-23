"""The HTML report of a sample: first whether flucut agrees with the
mapping, then what the mapping shows about the reads and the segments,
then the tables with every number, shown on request."""

from html import escape

from .agreement import depth_limits, sample_checks, segment_rows
from .html import (badge, card, chip, chips, heading, key_values, more,
                   note, num, page, pct, stacked_bar, table, tile, tiles,
                   write_page)
from .origins import MEASURED_ORIGIN, ORIGINS, origin_counts
from .report_agreement import agreement_section
from .report_compare import flucut_details, irma_card, irma_details
from .report_details import measured_details, terms_table
from .report_segments import segments_card
from .report_subtype import subtype_card

READ_CLASSES = [
    ("influenza", "influenza", "c-full"),
    ("offtarget", "other amplicons", "c-p1"),
    ("one_primer", "not mapped, one primer", "c-short"),
    ("no_primer", "not mapped, no primer", "c-oth"),
    ("dimer", "primer dimers", "c-nw")]
FLAG_INTERNAL, FLAG_DELETION, FLAG_OFFTARGET = 0.1, 0.5, 0.05
SITES_SHOWN = 5


def _head(s: dict) -> str:
    inputs = s["inputs"]
    f = s["flucut"]
    primers = ", ".join(f"{p['name']} ({p['sequence']})"
                        for p in s["settings"]["primers"])
    rows = [["reads", ", ".join(inputs["reads"])],
            ["IRMA results", ", ".join(inputs["irma"]) or "none"],
            ["flucut reads", inputs["flucut_reads"] or (
                "the reads themselves" if f and
                f["reads_are_flucut_output"] else "none")],
            ["flucut summary", inputs["flucut_summary"] or "none"],
            ["references", inputs["references"]],
            ["primers", primers],
            ["settings", f"primer errors up to "
             f"{s['settings']['error_rate']:.0%} of the primer; terminal: "
             f"within {s['settings']['terminal_window']} bases of the "
             f"terminus; deletion: at least "
             f"{s['settings']['min_deletion']} bases"]]
    sequencing = ((f or {}).get("summary") or {}).get("sequencing") or {}
    rows += [[f"sequencing: {key.replace('_', ' ')}", str(value)]
             for key, value in sequencing.items() if value]
    return (f"<h1><small>flucut post-hoc</small>{escape(s['sample'])}</h1>"
            f'<p class="sub">Measured from the base-by-base mapping of '
            f"every read · flucut-posthoc {s['version']} · "
            f"{escape(s['updated'])}</p>"
            + more("Inputs and settings",
                   table(["", ""], rows, cls="params")))


def _flags(s: dict, rows: list[dict]) -> str:
    """The state of the segments as measured, and what deserves
    attention."""
    found = [row for row in rows if row["measured"]["state"] != "missing"]
    complete = sum(row["measured"]["state"] == "complete" for row in found)
    min_depth, fraction = depth_limits(s)
    out = [badge("ok" if found and complete == len(rows) else "warn",
                 f"{complete}/{len(rows)} segments complete")]
    reads, classes = s["reads"]["total"], s["reads"]["classes"]
    if reads and classes.get("offtarget", 0) >= FLAG_OFFTARGET * reads:
        out.append(chip(f"{pct(classes['offtarget'], reads, 0)} of reads "
                        f"are amplicons of another template"))
    for key, limit, what in (
            ("internal", FLAG_INTERNAL, "of depth internally primed"),
            ("deletion", FLAG_DELETION,
             "of depth from reads with an internal deletion")):
        named = [f"{row['gene']} {row['measured'][key]:.0%}"
                 for row in found if (row["measured"][key] or 0) >= limit]
        if named:
            out.append(chip(f"{', '.join(named)} {what}"))
    return chips(*out) + note(
        f"Complete: at least {fraction:.0%} of the positions at a depth "
        f"of {min_depth}× or more, measured.")


def _read_lengths(s: dict) -> str:
    """The bases of the reads and, from the summary of flucut, their
    usual length and quality."""
    text = f"{s['reads']['bases'] / 1e6:,.1f} Mb"
    summary = ((s["flucut"] or {}).get("summary") or {})
    lengths = (summary.get("read_lengths") or {}).get("input")
    quality = (summary.get("read_quality") or {}).get("input")
    if lengths:
        text += f", median {lengths['median']:,} b, N50 {lengths['n50']:,}"
    if quality and quality.get("median") is not None:
        text += f", median quality Q{quality['median']:.0f}"
    return text


def _tiles(s: dict) -> str:
    reads, classes = s["reads"]["total"], s["reads"]["classes"]
    amplicons = s["amplicons"]
    total = amplicons["total"]
    complete = sum(c["amplicons"] for c in s["complete_amplicons"])
    left = s["ends"].get("left", {})
    primed = sum(sum(v.values()) for k, v in left.items() if k != "none")
    internal = sum(v.get("internal", 0) for k, v in left.items()
                   if k != "none")
    deletion = amplicons["categories"].get("deletion", 0)
    out = [tile(num(reads), "reads", _read_lengths(s))]
    subtype = s.get("subtype")
    if subtype:
        genes = subtype["genes"]
        out.append(tile(subtype["measured"] or "–", "subtype, mapping",
                        f"HA {genes['HA']['reads']:,}, NA "
                        f"{genes['NA']['reads']:,} of the first "
                        f"{subtype['sampled_reads']:,} reads"))
    out += [
        tile(pct(classes.get("influenza", 0), reads), "map to influenza",
             f"{classes.get('influenza', 0):,} of {reads:,} reads"),
        tile(num(total), "influenza amplicons", "a read can hold several"),
        tile(pct(complete, total), "with both primers",
             f"{complete:,} of {total:,} amplicons"),
        tile(pct(internal, primed), "5′ ends internally primed",
             f"{internal:,} of {primed:,} 5′ ends with a primer"),
        tile(pct(deletion, total), "with an internal deletion",
             f"{deletion:,} of {total:,} amplicons")]
    return tiles(*out)


def _reads_card(s: dict) -> str:
    reads, amplicons = s["reads"], s["amplicons"]
    classes = reads["classes"]
    summary = ((s["flucut"] or {}).get("summary") or {})
    said = summary.get("reads") or {}
    parts = [stacked_bar("What the reads are", "reads",
                         [(cls, label, classes.get(key, 0))
                          for key, label, cls in READ_CLASSES])]
    counts = origin_counts(amplicons["categories"], MEASURED_ORIGIN)
    order = ORIGINS[1:] + ORIGINS[:1]
    parts.append(stacked_bar(
        "Amplicons by origin, measured", "influenza amplicons",
        [(f"c-{key}", label, counts[key]) for key, label in order]))
    written = summary.get("written_by_category")
    if written:
        names = {"full": "full", "int": "internal", "del": "deletion",
                 "oth": "other"}
        parts.append(stacked_bar(
            "Reads by origin, flucut", "reads after splitting",
            [(f"c-{key}", label, written.get(names[key], 0))
             for key, label in order]
            + [("c-nw", "not written", said.get("too_short", 0))]))
    primers = amplicons["primers"]
    off = classes.get("offtarget", 0) + classes.get("dimer", 0)
    parts.append(stacked_bar(
        "Primers found, measured", "amplicons and reads not mapped",
        [("c-full", "both", primers.get("2", 0) + off),
         ("c-p1", "one", primers.get("1", 0) + classes.get("one_primer", 0)),
         ("c-oth", "none", primers.get("0", 0)
          + classes.get("no_primer", 0))]))
    if said:
        parts.append(stacked_bar(
            "Primers found, flucut", "reads after splitting",
            [("c-full", "both", said["both_primers"]),
             ("c-p1", "one", said["one_primer"]),
             ("c-oth", "none", said["no_primer"])]))
    per_read = reads["amplicons_per_read"]
    mapped = sum(per_read.values())
    several = sum(v for k, v in per_read.items() if k != "1")
    minus = amplicons.get("strand", {}).get("-", 0)
    items = [("several amplicons in the read", num(several),
              pct(several, mapped)),
             ("joined to another template",
              num(reads["joined_not_influenza"]),
              pct(reads["joined_not_influenza"], mapped)),
             ("amplicons on the minus strand", num(minus),
              pct(minus, amplicons["total"]))]
    items += [(kind, num(n), "") for kind, n in
              sorted(reads["junctions"].items())]
    parts.append(key_values(items))
    parts.append(note("Shares in the list: of the reads which map, of "
                      "the amplicons. The last lines count the junctions "
                      "between the amplicons of a read."))
    return card("Reads", *parts)


def _sites_card(s: dict) -> str:
    sites = s["internal_sites"]
    if not sites:
        return card("Internal priming sites",
                    note("No internally primed end was found."))
    with_flucut = bool((s["flucut"] or {}).get("origin"))
    total = s["amplicons"]["total"]

    def line(site: dict) -> list:
        row = [site["gene"], site["side"], site["position"], site["reads"],
               pct(site["reads"], total)]
        return row + [site["flucut_reads"]] if with_flucut else row

    headers = ["segment", "end", "position", "reads", "% of amplicons"] \
        + (["flucut reads"] if with_flucut else [])
    parts = [table(headers, [line(x) for x in sites[:SITES_SHOWN]])]
    if len(sites) > SITES_SHOWN:
        parts.append(more(
            f"{min(len(sites), 40) - SITES_SHOWN} smaller sites, of "
            f"{s['internal_sites_total']:,}",
            table(headers, [line(x) for x in sites[SITES_SHOWN:40]])))
    parts.append(note("Position: first base after the primer (left end), "
                      "last base before it (right end)."))
    return card("Internal priming sites", *parts)


def write_report(s: dict, path: str) -> None:
    """Write the HTML report of the summary `s`."""
    checks = sample_checks(s)
    body = [_head(s)]
    if s["flucut"]:
        body.append(agreement_section(s, checks))
    body += [heading("sample", "The sample, measured"),
             _flags(s, segment_rows(s)), _tiles(s),
             "" if s["flucut"] else subtype_card(s),
             f'<div class="cols">{_reads_card(s)}<div>{_sites_card(s)}'
             f'{irma_card(s)}</div></div>', segments_card(s),
             heading("details", "All numbers"), measured_details(s)]
    if s["flucut"]:
        body.append(flucut_details(s))
    body += [irma_details(s), terms_table(s)]
    write_page(path, page(f"flucut post-hoc: {s['sample']}", *body))
