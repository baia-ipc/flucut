"""What flucut said about each read, and its comparison with the mapping.

flucut appends tags to the titles of the reads it writes: the primers
found (`fl`, `fr`), the origin of each end (`ol`, `or`), where the bases
after a primer mapped (`ml`, `mr`), the category of the read (`ct`), the
bases removed (`tl`, `tr`) and the part of a split read (`sp`). They are
its predictions for each read; here they are compared with what the
mapping of the same read shows.
"""

from collections import Counter

from .amplicons import Unit, Unmapped
from .seqio import fastq_records, read_name

# the categories of the coverage plot of flucut
FLUCUT_CATEGORIES = ["terminal+", "terminal-", "internal+", "internal-",
                     "no_primer", "not_located", "deletion"]
ENDS = (("left", "12", "fl", "ol", "ml"), ("right", "13", "fr", "or", "mr"))
# the categories of the coverage plot (without the strand) of the values
# of the `ct` tag which name one
CT_CATEGORIES = {"full": "terminal", "internal": "internal"}


def load_flucut_reads(path: str, by_record: bool = False
                      ) -> dict[str, list[str]]:
    """The records of a FASTQ file written by flucut, by name of the
    original read: for each part, its kept length and its tags. With
    `by_record` each record is entered under its own name (for the
    analysis of the reads written by flucut themselves, where each part
    of a split read is a read)."""
    parts: dict[str, list[tuple[int, str]]] = {}
    for title, seq in fastq_records(path):
        name = read_name(title)
        tags = title.partition("\t")[2]
        k = 1
        if not by_record and "\tsp:Z:" in "\t" + tags:
            k = int(tags.rpartition("sp:Z:")[2].partition("/")[0])
            name = name.rpartition("_")[0]
        parts.setdefault(name, []).append((k, f"{len(seq)}\t{tags}"))
    return {name: [p for _, p in sorted(items)]
            for name, items in parts.items()}


def has_flucut_tags(path: str) -> bool:
    """Whether the first read of a FASTQ file has the tags of flucut."""
    for title, _ in fastq_records(path):
        return "\tfs:Z:" in title
    return False


def parse_part(text: str, lift: dict | None = None) -> dict:
    """A part as stored by load_flucut_reads(): the tags by name, `len`
    the kept length; `ml`/`mr` as (reference, 0-based position, reverse).
    A position on a reference of `lift` (see references.position_maps())
    is given on the reference chosen for the mapping instead."""
    fields = text.split("\t")
    part: dict = {"len": int(fields[0])}
    for field in fields[1:]:
        tag, _, value = field.partition(":")
        value = value[2:]
        if tag in ("ml", "mr"):
            ref, pos, strand = value.rsplit(":", 2)
            at = int(pos) - 1
            if lift and ref in lift:
                ref, positions = lift[ref]
                at = int(positions[min(max(at, 0), len(positions) - 1)])
            part[tag] = (ref, at, strand == "-")
        elif tag in ("tl", "tr", "ni"):
            part[tag] = int(value)
        else:
            part[tag] = value
    return part


def predicted_spans(part: dict, ref_lengths: dict[str, int]
                    ) -> tuple[str, str | None, list[tuple[int, int]]]:
    """The category and the span on its reference which flucut gives a
    read in its coverage plot: the category, the reference and the
    intervals. The span is worked out from the positions of the ends
    (`ml`, `mr`), as flucut does. The category is the one flucut wrote
    in the `ct` tag, with the strand; for reads without that tag it is
    worked out from the other tags too. A read which flucut placed in
    another way (by a probe from inside it, or by its length) has no
    tag for its position: it has no span, and without `ct` it is
    "not_located"."""
    category, ref, spans = _spans_from_tags(part, ref_lengths)
    said = part.get("ct")
    if said is None:
        return category, ref, spans
    if said in CT_CATEGORIES:
        strand = "-" if part.get("fs") == "-" else "+"
        return CT_CATEGORIES[said] + strand, ref, spans
    if said == "deletion":
        return "deletion", ref, spans
    return ("no_primer" if spans else "not_located"), ref, spans


def _spans_from_tags(part: dict, ref_lengths: dict[str, int]
                     ) -> tuple[str, str | None, list[tuple[int, int]]]:
    """The category and the span of a read from the tags of its ends
    (as read-cover.c does for a read with a located end)."""
    span = part["len"]
    minus = part.get("fs") == "-"
    strand = "-" if minus else "+"
    left, right = part.get("ml"), part.get("mr")
    both = False
    if left:
        ref, start, reverse = left
        both = right is not None and right[0] == ref
        if both:
            start, end = sorted((start, right[1]))
            end += 1
        elif reverse:
            start, end = start + 1 - span, start + 1
        else:
            end = start + span
        internal = part.get("ol") == "internal"
    elif right:
        ref, pos, reverse = right
        start, end = (pos + 1 - span, pos + 1) if reverse else \
            (pos, pos + span)
        if "fl" not in part:
            return "no_primer", ref, [(start, end)]
        if ref not in ref_lengths:
            return "not_located", None, []
        distance = start if reverse else ref_lengths[ref] - end
        internal = distance > 60 + span // 20
    else:
        return "not_located", None, []
    distance = end - start
    if both and distance - span > 100 and distance - span > distance / 10:
        half = span // 2
        return "deletion", ref, [(start, start + half),
                                 (end - (span - half), end)]
    return ("internal" if internal else "terminal") + strand, ref, \
        [(start, end)]


def _delta_bucket(delta: int) -> str:
    delta = abs(delta)
    for limit, label in ((0, "0"), (2, "1-2"), (5, "3-5"), (20, "6-20")):
        if delta <= limit:
            return label
    return ">20"


def _trim_bucket(origin: str, delta: int, terminus_len: int) -> str:
    """How the bases removed by flucut at an end compare with the primer
    found by the mapping; `delta`: removed minus non-template bases."""
    if abs(delta) <= 2:
        return "whole primer removed"
    if delta > 2:
        return "template removed"
    if origin == "terminal":
        if abs(delta + terminus_len) <= 2:
            return "terminus kept"
        if delta > -terminus_len:
            return "part of the terminus removed"
        return "tail or adapter bases left"
    return "primer bases left"


def _compare_unit(part: dict, category: str, unit: Unit, single: bool,
                  counts: Counter, trimmed: bool) -> None:
    for end, kind, found, origin, mapped in ENDS:
        said = part.get(origin, "unknown") if found in part else "no primer"
        truth_end = unit.end_of_kind(kind)
        if trimmed and found in part:
            # the primer is no longer in the read: the end where flucut
            # found it is judged by its position alone
            truth_end = unit.left if kind == "12" else unit.right
        truth = truth_end.origin if truth_end else "no primer"
        counts["origin", end, said, truth] += 1
        if truth_end and mapped in part:
            ref, pos, _ = part[mapped]
            bucket = _delta_bucket(pos - truth_end.pos) \
                if ref == unit.ref.name else "other reference"
            counts["position", end, bucket] += 1
    counts["category", category, unit.category] += 1
    if not single or trimmed:
        return
    start, stop = unit.read_span
    first, last = (unit.right, unit.left) if unit.reverse else \
        (unit.left, unit.right)
    for removed, expected, end in ((part.get("tl", 0), start, first),
                                   (part.get("tr", 0),
                                    unit.read_len - stop, last)):
        if end.hit is None:
            continue
        delta = removed - expected
        bucket = _trim_bucket(end.origin, delta, end.hit.oligo.terminus_len)
        counts["trim", end.origin, bucket] += 1
        if bucket == "primer bases left":
            counts["primer_bases_left", min(-delta, 30)] += 1


def compare_read(parts: list[str] | None, layout: list[Unit | None],
                 unmapped: Unmapped | None, ref_lengths: dict[str, int],
                 lift: dict, counts: Counter, trimmed: bool = False
                 ) -> list[tuple[str, str, int, int]]:
    """Count the agreement between the tags of a read and its mapping
    (`layout`: its amplicons in the order of the read, None for one which
    does not map); returns the coverage flucut predicts for it:
    (reference, category, from, to) intervals. `trimmed`: the reads
    mapped are those written by flucut, without their primers."""
    truth_class = "unmapped: " + unmapped.kind if unmapped else "influenza"
    if parts is None:
        counts["written", "not written", truth_class] += 1
        return []
    counts["written", "written", truth_class] += 1
    counts["split", min(len(parts), 4), min(len(layout), 4)] += 1
    predicted = []
    for i, text in enumerate(parts):
        part = parse_part(text, lift)
        category, ref, spans = predicted_spans(part, ref_lengths)
        for end, _, _, origin, mapped in ENDS:
            if part.get(origin) == "internal" and mapped in part:
                counts["flucut_site", part[mapped][0], end,
                       part[mapped][1] + 1] += 1
        predicted.extend((ref, category, a, b) for a, b in spans)
        if len(parts) == len(layout) and layout[i] is not None:
            _compare_unit(part, category, layout[i], len(parts) == 1,
                          counts, trimmed)
            continue
        if unmapped:
            other = truth_class
        elif len(parts) == len(layout):
            other = "joined amplicon, not influenza"
        else:
            other = "split differently"
        counts["category", category, other] += 1
        for end, _, found, origin, _ in ENDS:
            said = part.get(origin, "unknown") if found in part \
                else "no primer"
            counts["origin", end, said, other] += 1
    return predicted
