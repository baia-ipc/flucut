"""From the alignments of a read to its amplicons.

A read is one amplicon, or several joined end to end, or none (a read
which does not map: a primer dimer, an amplicon of another template).
Each amplicon ("unit") is described by its two ends on the reference:
where its template begins and ends, whether a primer precedes it, and
whether that is at a terminus of the segment or inside it.

Left and right are those of the reference (its plus strand): the left
end of a unit is at the lower coordinate, whatever the orientation in
which the read was sequenced.
"""

from dataclasses import dataclass, field

import edlib

from .mapping import Alignment
from .primers import (MIN_PRIMER_BASES, N_CONTROLS, Hit, Oligo,
                      contains_primer, search_before, search_left,
                      search_read_start)
from .references import Reference
from .seqio import revcomp

_ANY = [("N", base) for base in "ACGT"]

CATEGORIES = ["terminal+", "terminal-", "internal+", "internal-",
              "internal3", "noprimer_terminal", "noprimer_fragment",
              "deletion"]


@dataclass
class Params:
    oligos: tuple[Oligo, ...]
    error_rate: float = 0.2
    terminal_window: int = 10    # distance from a terminus still terminal
    gap: int = 10                # clipped template between primer and
    edge_gap: int = 40           # alignment; at an incomplete reference end
    tail_len: int = 10           # bases before a primer which are reported
    short_flank: tuple = (6, 60)   # bases of a flank too short to map
    min_aligned: int = 30        # shorter alignments are ignored
    min_deletion: int = 100      # shorter deletions are indels of a read
    read_window: int = 100       # primer search in a read not mapped


@dataclass
class End:
    pos: int                 # reference position of the template base
    q: int                   # left: first template base of the oriented
                             # read; right: one past the last
    hit: Hit | None
    origin: str              # "terminal" or "internal"
    control: tuple = ()      # which controls were found here
    before: str = ""         # the bases of the read before the primer

    @property
    def kind(self) -> str:
        return self.hit.oligo.kind if self.hit else "none"


@dataclass
class Unit:
    ref: Reference
    reverse: bool
    blocks: list[tuple[int, int]]
    left: End
    right: End
    read_len: int
    deletions: list[tuple[int, int]] = field(default_factory=list)
    aligned: int = 0
    errors: int = 0

    @property
    def template_len(self) -> int:
        return self.right.q - self.left.q

    @property
    def read_span(self) -> tuple[int, int]:
        """The template in the read as sequenced: [from, to)."""
        if self.reverse:
            return self.read_len - self.right.q, self.read_len - self.left.q
        return self.left.q, self.right.q

    @property
    def read_extent(self) -> tuple[int, int]:
        """The template with its primers in the read as sequenced."""
        n = self.read_len
        low = self.left.hit.start if self.left.hit else self.left.q
        high = n - self.right.hit.start if self.right.hit else self.right.q
        return (n - high, n - low) if self.reverse else (low, high)

    @property
    def n_primers(self) -> int:
        return (self.left.hit is not None) + (self.right.hit is not None)

    def end_of_kind(self, kind: str) -> End | None:
        """The end with a primer of the terminus `kind`, the usual side
        first (12: left, 13: right)."""
        ends = (self.left, self.right) if kind == "12" else \
            (self.right, self.left)
        for end in ends:
            if end.kind == kind:
                return end
        return None

    @property
    def category(self) -> str:
        strand = "-" if self.reverse else "+"
        if self.deletions:
            return "deletion"
        if self.left.hit and self.left.origin == "internal":
            return "internal" + strand
        if self.right.hit and self.right.origin == "internal":
            return "internal3"
        if self.left.hit:
            return "terminal" + strand
        if self.left.origin == "terminal":
            return "noprimer_terminal"
        return "noprimer_fragment"


def _allowed_gap(distance: int, terminus: int, p: Params) -> int:
    """The bases of template allowed between a primer and the alignment
    of the read, for an alignment which begins `distance` bases inside
    the segment. Far from the terminus only a few (`gap`: bases which
    the aligner clipped). Near it, the sample can differ so much from
    the reference that the alignment begins some bases inside: the
    primer may then be anywhere out to the terminus; and a reference
    without its terminus (`terminus` 0) may lack more bases of the
    segment end, which the read has beyond the reference."""
    if distance > p.edge_gap:
        return p.gap
    return max(0, distance) + (p.gap if terminus else p.edge_gap)


def _searched(aln: Alignment, ref: Reference, p: Params
              ) -> tuple:
    """The primers at the left and right boundaries of an alignment (the
    right one in the coordinates of the reverse complement), whether the
    control was found at each, and the bases before each primer."""
    seq, n = aln.seq, len(aln.seq)
    lgap = _allowed_gap(aln.rs - ref.first_template, ref.t12, p)
    rgap = _allowed_gap(ref.last_template - (aln.re - 1), ref.t13, p)
    rc = revcomp(seq)
    left = search_left(seq, aln.qs, p.oligos, lgap, p.error_rate)
    right = search_left(rc, n - aln.qe, p.oligos, rgap, p.error_rate)
    lctl = tuple(search_left(seq, aln.qs, p.oligos, lgap, p.error_rate, i)
                 is not None for i in range(N_CONTROLS))
    rctl = tuple(search_left(rc, n - aln.qe, p.oligos, rgap, p.error_rate,
                             i) is not None for i in range(N_CONTROLS))
    # a flank too short to map: without a primer at a boundary, or with
    # one some bases away from it, which may be the bases of the segment
    # end beyond a deletion rather than template the aligner clipped
    lflank = rflank = None
    if left is None or aln.qs - left.end >= p.short_flank[0]:
        lflank = _short_flank(seq, aln.qs, ref.seq, ref.first_template,
                              aln.rs, p)
        left = None if lflank else left
    if right is None or n - aln.qe - right.end >= p.short_flank[0]:
        rflank = _short_flank(rc, n - aln.qe, revcomp(ref.seq), ref.t13,
                              ref.length - aln.re, p)
        right = None if rflank else right
    ltail = seq[max(0, left.start - p.tail_len):left.start] if left else ""
    rtail = rc[max(0, right.start - p.tail_len):right.start] if right \
        else ""
    return left, right, lctl, rctl, ltail, rtail, lflank, rflank


def _flank_errors(length: int) -> int:
    """Differences allowed between a short flank and the segment end: a
    fifth of its bases, but fewer for the shortest, which would otherwise
    fit the 20 or so positions searched by chance."""
    if length >= 12:
        return length // 5
    return 1 if length >= 9 else 0


def _short_flank(seq: str, boundary: int, ref_seq: str, first_template: int,
                 aligned_from: int, p: Params
                 ) -> tuple[Hit, int, int] | None:
    """A flank too short to be mapped, at the start of a read with an
    internal deletion: a primer some bases before the alignment (at
    `boundary` of `seq`, from `aligned_from` of the reference), and the
    bases between them at the terminus of the segment (`first_template`),
    at least `min_deletion` bases before the alignment. Returns the
    primer and the interval of the reference of those bases."""
    low, high = p.short_flank
    hit = search_before(seq, boundary, p.oligos, low, high, p.error_rate)
    if hit is None:
        return None
    flank = seq[hit.end:boundary]
    # a reference without its terminus can lack more bases of the segment
    # end: up to 8 bases of the flank may lie before its first base
    start = first_template - 8
    pad = max(0, -start)
    target = "N" * pad + ref_seq[max(0, start):
                                 first_template + len(flank) + 8]
    found = edlib.align(flank, target, mode="HW", task="locations",
                        k=_flank_errors(len(flank)),
                        additionalEqualities=_ANY)
    if found["editDistance"] < 0:
        return None
    a, b = found["locations"][0]
    if aligned_from - (start + b + 1) < p.min_deletion:
        return None
    return hit, max(0, start + a), start + b + 1


def _joined(first: Alignment, second: Alignment, first_right: Hit | None,
            second_left: Hit | None) -> bool:
    """Whether two alignments, in the order of the oriented read, are
    parts of one template: same reference and strand, in the same order
    in the read and on the reference, no primer between them. What lies
    between them is a deletion, or a stretch too different from the
    reference to align (see _make_unit()). The two can overlap in the
    read (bases at a junction which fit both flanks). A primer between
    them needs room: with fewer than MIN_PRIMER_BASES bases of the read
    between the alignments, a primer found there is template which
    resembles it."""
    if first.ref != second.ref or first.reverse != second.reverse:
        return False
    if second.qs - first.qe >= MIN_PRIMER_BASES and (
            first_right is not None or second_left is not None):
        return False
    return (second.rs >= first.re - 30 and second.qs > first.qs
            and second.qe > first.qe)


def _make_unit(chain: list, ref: Reference, p: Params) -> Unit:
    """A unit from alignments in the order of the oriented read, each
    with its searched primers."""
    first, last = chain[0][0], chain[-1][0]
    n = len(first.seq)
    lhit, lctl = chain[0][1], chain[0][3]
    rhit, rctl = chain[-1][2], chain[-1][4]
    if lhit is not None:
        q = min(lhit.end, first.qe - 1)
        left_pos = first.ref_of(q, True)
    else:
        q, left_pos = first.qs, first.rs
    origin = "terminal" if left_pos <= ref.first_template + \
        p.terminal_window else "internal"
    left = End(left_pos, q, lhit, origin, lctl, chain[0][5])
    if rhit is not None:
        q = max(n - rhit.end, last.qs + 1)
        right_pos = last.ref_of(q - 1, False)
    else:
        q, right_pos = last.qe, last.re - 1
    origin = "terminal" if right_pos >= ref.last_template - \
        p.terminal_window else "internal"
    right = End(right_pos, q, rhit, origin, rctl, chain[-1][6])
    unit = Unit(ref, first.reverse, [], left, right, n)
    # flanks too short to map (see _short_flank()): the end of the unit
    # is then at the terminus, with the primer found there
    lflank, rflank = chain[0][7], chain[-1][8]
    if lhit is None and lflank is not None:
        hit, start, end = lflank
        unit.left = End(start, hit.end, hit, "terminal", lctl)
        unit.blocks.append((start, end))
        unit.deletions.append((end, first.rs))
    if rhit is None and rflank is not None:
        hit, start, end = rflank          # on the minus strand
        unit.right = End(ref.length - 1 - start, n - hit.end, hit,
                         "terminal", rctl)
    # between two aligned pieces, the bases of the reference without a
    # counterpart in the read are a deletion; a stretch present in both
    # which did not align is template, and counted as covered
    prev_q = prev_r = None
    for aln, *_ in chain:
        for qs, qe, rs, re in aln.pieces(p.min_deletion):
            if prev_r is not None and qs < prev_q and \
                    rs - prev_r >= p.min_deletion:
                # bases of the read which fit both flanks of a deletion
                # were aligned to both: they are counted once, with the
                # second flank (as flucut does at a junction)
                start, end = unit.blocks[-1]
                prev_r = max(start, end - (prev_q - qs))
                unit.blocks[-1] = (start, prev_r)
            if prev_r is None:
                unit.blocks.append((rs, re))
            elif (rs - prev_r) - max(0, qs - prev_q) >= p.min_deletion:
                unit.deletions.append((prev_r, rs))
                unit.blocks.append((rs, re))
                if aln.bq[0] != qs:      # a deletion inside the alignment
                    unit.errors -= rs - prev_r
            else:
                unit.blocks[-1] = (unit.blocks[-1][0],
                                   max(re, unit.blocks[-1][1]))
            prev_q, prev_r = qe, max(re, prev_r or re)
        unit.aligned += sum(aln.bl)
        unit.errors += aln.mismatches
    if rhit is None and rflank is not None:
        _, start, end = rflank
        unit.deletions.append((last.re, ref.length - end))
        unit.blocks.append((ref.length - end, ref.length - start))
    unit.errors = max(0, unit.errors)
    return unit


def _parts(units: list[Unit], seq: str, p: Params) -> list[Unit | None]:
    """The amplicons of a read in its order: its units, and None for each
    stretch of the read outside them which contains a primer and at least
    `min_aligned` other bases: an amplicon which does not map (not
    influenza) joined to the others. `seq`: the read as sequenced."""
    parts: list[Unit | None] = []
    at = 0
    for unit in sorted(units, key=lambda u: u.read_extent[0]):
        low, high = unit.read_extent
        gap = seq[at:low]
        if len(gap) >= p.min_aligned + 12 and \
                contains_primer(gap, p.oligos):
            parts.append(None)
        parts.append(unit)
        at = max(at, high)
    gap = seq[at:]
    if len(gap) >= p.min_aligned + 12 and contains_primer(gap, p.oligos):
        parts.append(None)
    return parts


def build_units(alignments: list[Alignment], refs: dict[str, Reference],
                p: Params
                ) -> tuple[list[Unit], list[str], list[Unit | None]]:
    """The units of a read, in the order of the read as sequenced, the
    kind of each junction between consecutive units, and all the
    amplicons of the read (see _parts())."""
    alns = [a for a in alignments if a.qe - a.qs >= p.min_aligned]
    alns.sort(key=lambda a: a.read_start)
    searched = [(a, *_searched(a, refs[a.ref], p)) for a in alns]
    chains: list[list] = []
    for item in searched:
        if chains:
            prev = chains[-1][-1] if not item[0].reverse else chains[-1][0]
            if item[0].reverse:
                ok = _joined(item[0], prev[0], item[2], prev[1])
            else:
                ok = _joined(prev[0], item[0], prev[2], item[1])
            if ok:
                if item[0].reverse:
                    chains[-1].insert(0, item)
                else:
                    chains[-1].append(item)
                continue
        chains.append([item])
    units = [_make_unit(c, refs[c[0][0].ref], p) for c in chains]
    junctions = []
    for a, b in zip(units, units[1:]):
        facing = (a.left if a.reverse else a.right,
                  b.right if b.reverse else b.left)
        if facing[0].hit or facing[1].hit:
            junctions.append("amplicons joined end to end")
        elif a.ref is b.ref:
            junctions.append("same segment, rearranged")
        else:
            junctions.append("different segments, no primer")
    parts: list[Unit | None] = []
    if units:
        first = alns[0]
        seq = revcomp(first.seq) if first.reverse else first.seq
        parts = _parts(units, seq, p)
    return units, junctions, parts


@dataclass
class Unmapped:
    """A read without an alignment: the primers at its two ends."""
    kind: str                # dimer, offtarget, one_primer, no_primer
    template_len: int
    start: Hit | None
    end: Hit | None
    controls: tuple          # hits of each control at the two ends
    template: str            # after a 12 primer, oriented from it ("")


def classify_unmapped(seq: str, p: Params) -> Unmapped:
    rc = revcomp(seq)
    start = search_read_start(seq, p.oligos, p.read_window, p.error_rate)
    end = search_read_start(rc, p.oligos, p.read_window, p.error_rate)
    controls = tuple(sum(
        search_read_start(s, p.oligos, p.read_window, p.error_rate, i)
        is not None for s in (seq, rc)) for i in range(N_CONTROLS))
    n = len(seq)
    if start and end:
        length = n - start.end - end.end
        kind = "dimer" if length < p.min_aligned else "offtarget"
    else:
        length = n - (start.end if start else 0) - (end.end if end else 0)
        kind = "one_primer" if start or end else "no_primer"
    template = ""
    if start and start.oligo.kind == "12":
        template = seq[start.end:n - (end.end if end else 0)]
    elif end and end.oligo.kind == "12":
        template = rc[end.end:n - (start.end if start else 0)]
    return Unmapped(kind, max(0, length), start, end, controls, template)
