"""The summary of a run: from the counts to the document written as JSON
(posthoc_summary.json), from which the HTML report is made."""

import time
from collections import Counter

import numpy as np

from . import __version__
from .amplicons import CATEGORIES
from .analysis import LENGTH_BIN, LENGTH_MAX, SHORT_FRACTION
from .flucut_tags import FLUCUT_CATEGORIES
from .references import Reference
from .seqio import revcomp

SITE_DISTANCE = 5     # positions at most this far apart are one site
MAX_SITES = 200
_IUPAC = {"R": "AG", "Y": "CT", "N": "ACGT"}


def _select(counts: Counter, name: str) -> dict:
    """The counts whose key begins with `name`, by the rest of the key."""
    out = {}
    for key, value in counts.items():
        if isinstance(key, tuple) and key[0] == name:
            out[key[1:] if len(key) > 2 else key[1]] = value
    return out


def _nested(counts: Counter, name: str) -> dict:
    """The counts of `name` as nested dictionaries with string keys."""
    out: dict = {}
    for key, value in _select(counts, name).items():
        keys = key if isinstance(key, tuple) else (key,)
        node = out
        for k in keys[:-1]:
            node = node.setdefault(str(k), {})
        node[str(keys[-1])] = value
    return out


def cluster_sites(positions: dict[int, int]) -> list[dict]:
    """Group positions (position: reads) at most SITE_DISTANCE apart into
    sites; the position of a site is the one with most reads."""
    sites, current = [], []
    for pos in sorted(positions):
        if current and pos - current[-1] > SITE_DISTANCE:
            sites.append(current)
            current = []
        current.append(pos)
    if current:
        sites.append(current)
    return [{"position": max(group, key=lambda p: positions[p]),
             "from": group[0], "to": group[-1],
             "reads": sum(positions[p] for p in group)} for group in sites]


def _footprint(ref: Reference, side: str, position: int, terminus: str
               ) -> tuple[str, int]:
    """The bases of the reference where the terminus part of the primer
    annealed at a site, written as the primer reads, and their mismatches
    with it."""
    size = len(terminus)
    if side == "left":
        template = ref.seq[max(0, position - 1 - size):position - 1]
    else:
        template = revcomp(ref.seq[position:position + size])
    if len(template) < size:
        return template, -1
    mismatches = sum(b not in _IUPAC.get(p, p)
                     for p, b in zip(terminus, template))
    return template, mismatches


def _sites(counts: Counter, references: list[Reference], oligos
           ) -> list[dict]:
    termini = {o.kind: o.seq[-o.terminus_len:] for o in oligos}
    groups: dict[tuple, dict[int, int]] = {}
    for (ref, side, kind, pos), n in _select(counts, "site").items():
        groups.setdefault((ref, side, kind), {})[pos] = n
    said: dict[tuple, dict[int, int]] = {}
    for (ref, end, pos), n in _select(counts, "flucut_site").items():
        said.setdefault((ref, end), {})[pos] = n
    by_name = {r.name: r for r in references}
    sites = []
    for (ref, side, kind), positions in groups.items():
        flucut = said.get((ref, "left" if kind == "12" else "right"), {})
        for site in cluster_sites(positions):
            template, mismatches = _footprint(by_name[ref], side,
                                              site["position"],
                                              termini[kind])
            site.update(reference=ref, gene=by_name[ref].gene, side=side,
                        primer=kind, template=template,
                        mismatches=mismatches,
                        flucut_reads=sum(
                            n for p, n in flucut.items()
                            if site["from"] - SITE_DISTANCE <= p
                            <= site["to"] + SITE_DISTANCE))
            sites.append(site)
    sites.sort(key=lambda s: -s["reads"])
    return sites


def _flucut_only_sites(counts: Counter, sites: list[dict]) -> list[dict]:
    """The sites of flucut with no site of the mapping within
    SITE_DISTANCE."""
    groups: dict[tuple, dict[int, int]] = {}
    for (ref, end, pos), n in _select(counts, "flucut_site").items():
        groups.setdefault((ref, end), {})[pos] = n
    measured: dict[tuple, list[tuple[int, int]]] = {}
    for s in sites:
        end = "left" if s["primer"] == "12" else "right"
        measured.setdefault((s["reference"], end), []).append(
            (s["from"] - SITE_DISTANCE, s["to"] + SITE_DISTANCE))
    out = []
    for (ref, end), positions in groups.items():
        spans = measured.get((ref, end), [])
        for site in cluster_sites(positions):
            if not any(a <= site["position"] <= b for a, b in spans):
                site.update(reference=ref, end=end)
                out.append(site)
    out.sort(key=lambda s: -s["reads"])
    return out


def _histogram(bins: dict[int, int]) -> list[int]:
    if not bins:
        return []
    out = [0] * (max(bins) // LENGTH_BIN + 1)
    for start, n in bins.items():
        out[start // LENGTH_BIN] = n
    return out


def _main_length(histogram: list[int], below: int | None = None
                 ) -> int | None:
    part = histogram if below is None else \
        histogram[:max(0, -(-below // LENGTH_BIN))]
    if not part or max(part) == 0:
        return None
    return int(np.argmax(part)) * LENGTH_BIN


def _complete(counts: Counter, references: list[Reference],
              group: dict | None, lift: dict) -> list[dict]:
    said = {lift[g["group"]][0] if g["group"] in lift else g["group"]: g
            for g in (group or {}).get("complete_amplicons", [])}
    out = []
    for ref in references:
        bins = {k[1]: v for k, v in _select(counts, "length").items()
                if k[0] == ref.name}
        histogram = _histogram(bins)
        limit = int(SHORT_FRACTION * ref.length)
        entry = {"reference": ref.name, "gene": ref.gene,
                 "reference_length": ref.length,
                 "amplicons": counts["complete", ref.name],
                 "main_length": _main_length(histogram),
                 "short": counts["short", ref.name],
                 "short_main_length": _main_length(histogram, limit),
                 "histogram": histogram, "flucut": None}
        if ref.name in said:
            g = said[ref.name]
            entry["flucut"] = {k: g[k] for k in (
                "amplicons", "main_length", "short", "short_main_length",
                "histogram")}
        out.append(entry)
    return out


def _coverage(references, measured, predicted, shown, irma, lifted) -> dict:
    out = {}
    for i, ref in enumerate(references):
        depth = measured.depth(i)
        entry = {"measured": {c: depth[j].tolist()
                              for j, c in enumerate(CATEGORIES)},
                 "predicted": None, "irma": {}, "agreement": None}
        pred = None
        if shown is not None:
            if ref.name in shown:
                pred = np.array([shown[ref.name].get(c, [0] * ref.length)
                                 for c in FLUCUT_CATEGORIES])
        elif predicted is not None:
            pred = predicted.depth(i)
        if pred is not None:
            entry["predicted"] = {c: pred[j].tolist()
                                  for j, c in enumerate(FLUCUT_CATEGORIES)}
            a, b = depth.sum(axis=0), pred.sum(axis=0)
            r = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and \
                b.std() > 0 else None
            entry["agreement"] = {
                "measured_mean_depth": float(a.mean()),
                "predicted_mean_depth": float(b.mean()),
                "correlation": r,
                "mean_absolute_difference": float(np.abs(a - b).mean())}
        for run, lift in zip(irma, lifted):
            if ref.name in lift and lift[ref.name]["depth"] is not None:
                entry["irma"][run.label] = lift[ref.name]["depth"]
        out[ref.name] = entry
    return out


def _agreement(origin: dict) -> dict:
    """Scores of the origin calls of flucut, over the ends which both
    flucut and the mapping assign to a terminus or to the interior."""
    scores = {}
    for end, table in origin.items():
        get = lambda said, truth: table.get(said, {}).get(truth, 0)
        tt, ti = get("terminal", "terminal"), get("terminal", "internal")
        it, ii = get("internal", "terminal"), get("internal", "internal")
        ut, ui = get("unknown", "terminal"), get("unknown", "internal")
        scores[end] = {
            "called": tt + ti + it + ii,
            "agree": tt + ii,
            "internal_precision": ii / (ii + it) if ii + it else None,
            "internal_recall": ii / (ii + ti + ui) if ii + ti + ui else None,
            "terminal_precision": tt / (tt + ti) if tt + ti else None,
            "terminal_recall": tt / (tt + it + ut) if tt + it + ut else None}
    return scores


def _flucut(counts: Counter, sites: list[dict], group: dict | None,
            trimmed_input: bool) -> dict:
    origin = _nested(counts, "origin")
    return {
        "reads_are_flucut_output": trimmed_input,
        "written": _nested(counts, "written"),
        "split": {f"{f},{t}": n
                  for (f, t), n in _select(counts, "split").items()},
        "origin": origin,
        "scores": _agreement(origin),
        "position": _nested(counts, "position"),
        "category": _nested(counts, "category"),
        "trim": _nested(counts, "trim"),
        "primer_bases_left": {str(k): v for k, v in sorted(
            _select(counts, "primer_bases_left").items())},
        "sites_only_flucut": _flucut_only_sites(counts, sites)[:50],
        "predicted_on_other_references": sum(
            _select(counts, "predicted_elsewhere").values()),
        "summary": None if group is None else {
            "reads": group.get("reads"), "ends": group.get("ends"),
            "short_complete": group.get("short_complete")}}


def _controls(counts: Counter) -> dict:
    """For each place of the search: the searches, the hits of each
    control and the middle one of them."""
    out = {}
    for context, values in _nested(counts, "control").items():
        found = sorted(v for k, v in values.items() if k != "searched")
        out[context] = {"searched": values.get("searched", 0),
                        "found_each": found,
                        "found": found[len(found) // 2] if found else 0}
    return out


def _top(counts: Counter, name: str, fields: tuple, n: int) -> list[dict]:
    items = sorted(_select(counts, name).items(), key=lambda kv: -kv[1])
    return [dict(zip(fields, key), reads=value) for key, value in items[:n]]


def build_summary(config, params, references, evidence, counts, measured,
                  predicted, clusters, irma, lifted, group, has_flucut,
                  trimmed_input, lift, shown) -> dict:
    """`predicted`: the coverage rebuilt from the tags of flucut (or
    None); `shown`: the coverage read from its report (or None), which is
    preferred; `lift`: the position maps of the other references."""
    sites = _sites(counts, references, params.oligos)
    ref_entries = []
    for ref in references:
        aligned = counts["segment_bases", ref.name, "aligned"]
        errors = counts["segment_bases", ref.name, "errors"]
        ref_entries.append({
            "name": ref.name, "description": ref.description,
            "segment": ref.segment, "gene": ref.gene, "length": ref.length,
            "terminus_12_bases": ref.t12, "terminus_13_bases": ref.t13,
            "adapted_bases": ref.adapted,
            "amplicons": counts["segment_units", ref.name],
            "identity": 1 - errors / aligned if aligned else None,
            "irma": {run.label: {k: v for k, v in lift[ref.name].items()
                                 if k != "depth"}
                     for run, lift in zip(irma, lifted)
                     if ref.name in lift}})
    classes = _nested(counts, "class")
    return {
        "format": 1, "program": "flucut-posthoc", "version": __version__,
        "updated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "sample": config.name, "command": config.command,
        "inputs": {"reads": config.reads, "irma": config.irma,
                   "flucut_reads": config.flucut_reads,
                   "flucut_summary": config.flucut_summary,
                   "flucut_report": config.flucut_report,
                   "references": config.references},
        "settings": {
            "primers": [{"name": o.name, "sequence": o.seq}
                        for o in params.oligos],
            "error_rate": params.error_rate,
            "terminal_window": params.terminal_window,
            "min_deletion": params.min_deletion,
            "min_aligned": params.min_aligned,
            "length_bin": LENGTH_BIN, "length_max": LENGTH_MAX,
            "short_fraction": SHORT_FRACTION},
        "references": ref_entries,
        "reference_choice": evidence,
        "reads": {"total": counts["reads"], "bases": counts["bases"],
                  "classes": classes,
                  "amplicons_per_read": _nested(counts, "units_per_read"),
                  "junctions": _nested(counts, "junction"),
                  "joined_not_influenza": counts["joined_not_influenza"],
                  "unmapped_lengths": {
                      kind: _histogram({k[1]: v for k, v in _select(
                          counts, "unmapped_length").items()
                          if k[0] == kind})
                      for kind in ("dimer", "offtarget", "one_primer",
                                   "no_primer")}},
        "amplicons": {"total": sum(_select(counts, "unit").values()),
                      "categories": _nested(counts, "unit"),
                      "primers": _nested(counts, "unit_primers"),
                      "strand": _nested(counts, "unit_strand")},
        "ends": _nested(counts, "end"),
        "primer_search": {
            "errors": _nested(counts, "primer_errors"),
            "truncated": _nested(counts, "primer_truncated"),
            "before_primer": {
                name: [[seq, n] for seq, n in sorted(
                    tails.items(), key=lambda kv: -kv[1])[:5]]
                + [["all", sum(tails.values())]]
                for name, tails in _nested(counts,
                                           "before_primer").items()},
            "control": _controls(counts)},
        "terminal_offsets": _nested(counts, "terminal_offset"),
        "internal_sites": sites[:MAX_SITES],
        "internal_sites_total": len(sites),
        "fragment_ends": _top(counts, "fragment_end",
                              ("reference", "side", "position"), 30),
        "deletions": _top(counts, "deletion",
                          ("reference", "from", "to"), 40),
        "complete_amplicons": _complete(counts, references, group, lift),
        "coverage": _coverage(references, measured, predicted, shown, irma,
                              lifted),
        "predicted_coverage_from": "report" if shown is not None else
        "tags" if predicted is not None else None,
        "offtarget": clusters,
        "irma": [{"label": run.label, "path": run.path,
                  "module": run.module, "read_counts": run.read_counts,
                  "reads_by_kind": _nested(counts, "irma").get(run.label,
                                                               {}),
                  "same_segment": _nested(counts, "irma_segment").get(
                      run.label, {})}
                 for run in irma],
        "flucut": _flucut(counts, sites, group, trimmed_input)
        if has_flucut else None,
    }
