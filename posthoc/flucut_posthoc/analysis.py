"""Analysis of batches of mapped reads (run in worker processes).

Everything counted is a key of one Counter, so that the results of the
batches are merged by adding them. The keys are tuples whose first item
names the count (see process_batch()).
"""

from collections import Counter
from dataclasses import dataclass, field

from .amplicons import (Params, Unit, Unmapped, build_units,
                        classify_unmapped)
from .flucut_tags import compare_read
from .mapping import parse_record
from .references import Reference

LENGTH_BIN = 10
LENGTH_MAX = 4000
SHORT_FRACTION = 0.9

_state: dict = {}


@dataclass
class BatchResult:
    counts: Counter = field(default_factory=Counter)
    # coverage intervals: (reference index, category, from, to)
    coverage: list[tuple[int, str, int, int]] = field(default_factory=list)
    predicted: list[tuple[int, str, int, int]] = field(default_factory=list)
    rows: list[str] = field(default_factory=list)
    # reads not mapped which begin with a 12 primer: (first template
    # bases, template length, both primers, template)
    offtarget: list[tuple[str, int, bool, str]] = field(default_factory=list)


def init_worker(references: list[Reference], params: Params,
                ref_lengths: dict[str, int], lift: dict,
                trimmed: bool = False) -> None:
    _state["lift"] = lift
    _state["trimmed"] = trimmed
    _state["refs"] = {r.name: r for r in references}
    _state["index"] = {r.name: i for i, r in enumerate(references)}
    _state["params"] = params
    _state["ref_lengths"] = ref_lengths


def _end_text(end) -> str:
    if end.hit is None:
        return f"{end.pos + 1}\t.\t.\t{end.origin}"
    return f"{end.pos + 1}\t{end.hit.oligo.name}\t{end.hit.errors}" \
           f"\t{end.origin}"


def _count_end(counts: Counter, unit: Unit, side: str) -> None:
    end = unit.left if side == "left" else unit.right
    ref = unit.ref
    counts["end", side, end.kind, end.origin] += 1
    counts["control", "at an alignment boundary", "searched"] += 1
    for i, found in enumerate(end.control):
        counts["control", "at an alignment boundary", f"found {i}"] += found
    if end.hit is None:
        if end.origin == "internal":
            counts["fragment_end", ref.name, side, end.pos + 1] += 1
        return
    counts["primer_errors", end.hit.oligo.name, end.hit.errors] += 1
    if end.hit.missing:
        counts["primer_truncated", end.hit.oligo.name] += 1
    elif len(end.before) == _state["params"].tail_len:
        counts["before_primer", end.hit.oligo.name, end.before] += 1
    if end.origin == "internal":
        counts["site", ref.name, side, end.kind, end.pos + 1] += 1
    else:
        expected = ref.first_template if side == "left" else \
            ref.last_template
        counts["terminal_offset", side, end.pos - expected] += 1


def _count_unit(result: BatchResult, name: str, k: int, n: int,
                unit: Unit) -> None:
    counts, ref = result.counts, unit.ref
    category = unit.category
    counts["unit", category] += 1
    counts["segment_units", ref.name] += 1
    counts["segment_bases", ref.name, "aligned"] += unit.aligned
    counts["segment_bases", ref.name, "errors"] += unit.errors
    counts["unit_primers", unit.n_primers] += 1
    counts["unit_strand", "-" if unit.reverse else "+"] += 1
    _count_end(counts, unit, "left")
    _count_end(counts, unit, "right")
    length = unit.template_len
    if unit.left.kind == "12" and unit.right.kind == "13":
        counts["complete", ref.name] += 1
        counts["length", ref.name,
               min(length, LENGTH_MAX) // LENGTH_BIN * LENGTH_BIN] += 1
        if length < SHORT_FRACTION * ref.length:
            counts["short", ref.name] += 1
    for start, end in unit.deletions:
        counts["deletion", ref.name, start // 10 * 10, end // 10 * 10] += 1
    index = _state["index"][ref.name]
    for start, end in unit.blocks:
        result.coverage.append((index, category, start, end))
    deletions = ",".join(f"{a + 1}-{b}" for a, b in unit.deletions) or "."
    identity = 1 - unit.errors / unit.aligned if unit.aligned else 0
    result.rows.append(
        f"{name}\t{k}\t{n}\t{unit.read_len}\t{ref.name}"
        f"\t{'-' if unit.reverse else '+'}\t{_end_text(unit.left)}"
        f"\t{_end_text(unit.right)}\t{length}\t{deletions}\t{category}"
        f"\t{identity:.4f}\n")


def _count_unmapped(result: BatchResult, name: str, seq: str,
                    found: Unmapped) -> None:
    counts = result.counts
    counts["control", "at the end of a read not mapped", "searched"] += 2
    for i, n in enumerate(found.controls):
        counts["control", "at the end of a read not mapped",
               f"found {i}"] += n
    counts["unmapped_length", found.kind,
           min(found.template_len, LENGTH_MAX) // LENGTH_BIN
           * LENGTH_BIN] += 1
    if found.template and len(found.template) >= 24:
        result.offtarget.append((found.template[:24], len(found.template),
                                 found.kind == "offtarget",
                                 found.template))
    result.rows.append(f"{name}\t0\t0\t{len(seq)}\t.\t.\t.\t.\t.\t.\t."
                       f"\t.\t.\t.\t{found.template_len}\t.\t{found.kind}"
                       f"\t.\n")


def process_batch(batch: list[tuple[str, list[str], list[str] | None,
                                    tuple]]) -> BatchResult:
    """Analyse the reads of a batch: for each, its name, its SAM records,
    its parts in the output of flucut (None: not written, or no output of
    flucut given: see `flucut` below) and its use by each IRMA run."""
    refs, params = _state["refs"], _state["params"]
    result = BatchResult()
    counts = result.counts
    for name, lines, parts, irma in batch:
        alignments, seq = [], ""
        for line in lines:
            _, aln, seq = parse_record(line)
            if aln is not None:
                alignments.append(aln)
        units, junctions, layout = build_units(alignments, refs, params)
        counts["reads"] += 1
        counts["bases"] += len(seq)
        unmapped = None
        if units:
            counts["class", "influenza"] += 1
            counts["units_per_read", min(len(units), 4)] += 1
            if None in layout:
                counts["joined_not_influenza"] += 1
            for junction in junctions:
                counts["junction", junction] += 1
            for k, unit in enumerate(units, 1):
                _count_unit(result, name, k, len(units), unit)
            truth = units[0].category if len(units) == 1 else "several"
        else:
            unmapped = classify_unmapped(seq, params)
            counts["class", unmapped.kind] += 1
            _count_unmapped(result, name, seq, unmapped)
            truth = unmapped.kind
        for label, segment in irma:
            counts["irma", label, truth, segment is not None] += 1
            if segment is not None and units:
                same = any(u.ref.gene == segment for u in units)
                counts["irma_segment", label, same] += 1
        if parts is not False:
            index = _state["index"]
            for ref, category, start, end in compare_read(
                    parts, layout, unmapped, _state["ref_lengths"],
                    _state["lift"], counts, _state["trimmed"]):
                if ref in index:
                    result.predicted.append((index[ref], category, start,
                                             end))
                else:
                    counts["predicted_elsewhere", category] += 1
    return result
