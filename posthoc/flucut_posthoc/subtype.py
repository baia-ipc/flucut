"""The subtype of a sample: the one the mapping shows, the one flucut
assigned and the one IRMA assembled, and whether they are the same.

The mapping gives the H of the HA reference chosen for the sample and
the N of the NA reference chosen (references.choose_references()): the
subtype is read from the description of a reference, "(H5N1)". The
evidence is the number of sampled reads which map to the references of
each H and of each N. flucut is judged on HA and NA only; the references
of the other segments are compared as information."""

import re

from .references import GENES

MIN_READS = 20     # fewer sampled reads on HA or NA: the type is not told
PARTS = {"HA": "H", "NA": "N"}
_SUBTYPE = re.compile(r"\b(H\d+)(N\d+)\b")
_IRMA = re.compile(r"^4-A_(HA|NA)_([HN]\d+)$")
STATUS = {"same": "ok", "different": "bad", "unknown": "na"}


def subtype_of(description: str) -> str | None:
    """The subtype named in the description of a reference."""
    found = _SUBTYPE.search(description or "")
    return found.group(1) + found.group(2) if found else None


def part(subtype: str | None, letter: str) -> str | None:
    """The H or the N (`letter`) of a subtype, e.g. "H5" of "H5N1"."""
    found = re.search(letter + r"\d+", subtype or "")
    return found.group(0) if found else None


def _join(h: str | None, n: str | None) -> str | None:
    if h is None and n is None:
        return None
    return (h or "H?") + (n or "N?")


def irma_subtype(read_counts: list) -> str | None:
    """The subtype an IRMA run assembled: its HA and NA with most
    reads."""
    best: dict[str, tuple[int, str]] = {}
    for step, n in read_counts:
        found = _IRMA.match(step)
        if found and n > best.get(found.group(1), (0, ""))[0]:
            best[found.group(1)] = (n, found.group(2))
    return _join(best.get("HA", (0, None))[1], best.get("NA", (0, None))[1])


def _gene_evidence(gene: str, segment: int, evidence: list[dict],
                   subtypes: dict[str, str | None]) -> dict:
    """What the sampled reads say about the H (or the N): the reference
    chosen, the reads on the references of each type, and the type."""
    letter = PARTS[gene]
    entries = [e for e in evidence if e["segment"] == segment]
    chosen = next((e for e in entries if e["chosen"]), None)
    by_type: dict[str, int] = {}
    for e in entries:
        kind = part(subtypes.get(e["reference"]), letter) or "unknown"
        by_type[kind] = by_type.get(kind, 0) + e["sample_reads"]
    kind = part(subtypes.get(chosen["reference"]), letter) if chosen \
        else None
    reads = by_type.get(kind, 0) if kind else 0
    others = sorted(((n, k) for k, n in by_type.items() if k != kind),
                    reverse=True)
    return {"reference": chosen["reference"] if chosen else None,
            "type": kind if reads >= MIN_READS else None,
            "closest": kind, "reads": reads, "by_type": by_type,
            "other": {"type": others[0][1], "reads": others[0][0]}
            if others else None}


def _compare(measured: str | None, said: str | None) -> str:
    if measured is None or said is None:
        return "unknown"
    return "same" if measured == said else "different"


def _verdict(genes: dict, flucut: str | None) -> tuple[str, str]:
    """Whether the subtype of flucut is the one the mapping shows, and
    a line which says it."""
    measured = _join(genes["HA"]["type"], genes["NA"]["type"])
    states = [genes[gene]["status"] for gene in PARTS]
    if "different" in states:
        return "different", f"Subtype: flucut {flucut}, mapping {measured}"
    if states == ["same", "same"]:
        return "same", f"Subtype {flucut}: confirmed by the mapping"
    if flucut is None:
        why = "flucut assigned none"
    else:
        why = ", ".join(
            f"{gene} {genes[gene]['reads']} mapped reads" if
            genes[gene]["type"] is None else f"flucut gives no {letter}"
            for gene, letter in PARTS.items()
            if genes[gene]["status"] == "unknown")
    return "unknown", f"Subtype: cannot be told ({why})"


def _segments(references: list[dict], subtypes: dict,
              flucut_segments: dict) -> list[dict]:
    """For each segment, the reference of the mapping and that of flucut,
    with the subtype named in each."""
    rows = []
    for ref in references:
        said = (flucut_segments.get(ref.get("segment")) or {}) \
            .get("reference")
        rows.append({
            "gene": ref["gene"], "reference": ref["name"],
            "reference_subtype": subtypes.get(ref["name"]),
            "flucut_reference": said,
            "flucut_reference_subtype": subtypes.get(said),
            "same_reference": said == ref["name"] if said else None})
    return rows


def subtype_summary(evidence: list[dict], subtypes: dict[str, str | None],
                    references: list[dict], flucut_sample: dict | None,
                    irma: dict[str, list], sampled_reads: int) -> dict:
    """The subtype of a sample. `evidence`: of the choice of the
    references; `subtypes`: the subtype named by each reference;
    `references`: those chosen (name, gene, segment); `flucut_sample`:
    the member `sample` of the summary of flucut; `irma`: the read
    counts of each IRMA run by label."""
    flucut = (flucut_sample or {}).get("subtype")
    genes = {}
    for segment, gene in GENES.items():
        if gene in PARTS:
            entry = _gene_evidence(gene, segment, evidence, subtypes)
            entry["flucut"] = part(flucut, PARTS[gene])
            entry["status"] = _compare(entry["type"], entry["flucut"])
            genes[gene] = entry
    status, text = _verdict(genes, flucut)
    flucut_segments = {seg["segment"]: seg for seg in
                       (flucut_sample or {}).get("segments", [])}
    return {"status": status, "text": text,
            "measured": _join(genes["HA"]["type"], genes["NA"]["type"]),
            "flucut": flucut,
            "irma": {label: irma_subtype(counts)
                     for label, counts in irma.items()},
            "min_reads": MIN_READS, "sampled_reads": sampled_reads,
            "genes": genes,
            "segments": _segments(references, subtypes, flucut_segments)}
