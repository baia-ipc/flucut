"""The universal primers and their search in a read.

What is searched is the part of a primer which every scheme shares: the
conserved terminus of the segments (12 bases for Uni12, 13 for Uni13).
The tail which a scheme adds before it differs between laboratories, so
it is not assumed: the bases found before the terminus are reported
instead (see analysis.py), and show which tail, barcode or adapter the
reads really carry.

A primer is searched next to the boundary of an alignment of the read:
the mapping says where the template of the read begins, so the primer
must end there. This is what makes the search more specific than one in
a window at the read end: a primer-like sequence elsewhere in the read is
not taken for the primer which started the amplicon.

Each search is repeated with controls: the bases of the primer in
another order, sequences of the same composition which no amplicon
should contain (each at least 5 edits from the primers, their usual
tails and their reverse complements). The hits of a control measure how
often a primer is found by chance. A control can still resemble an
adapter or a barcode of a library, and then be found in most reads:
three are used for each primer, and the middle one of their rates is
reported.
"""

from dataclasses import dataclass

import edlib

from .seqio import revcomp

# degenerate bases of the primers, and N (padding, see search_left())
_EQUALITIES = [("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("N", "A"),
               ("N", "C"), ("N", "G"), ("N", "T"), ("N", "R")]

# a read can begin inside a primer: at least this many of its bases must
# be present
MIN_PRIMER_BASES = 8


@dataclass(frozen=True)
class Oligo:
    """The conserved terminus of a primer, 5' to 3'; `kind` is the end of
    the segment it binds ("12": the start of the plus strand, "13": its
    end)."""
    name: str
    seq: str
    kind: str
    controls: tuple[str, ...]

    @property
    def terminus_len(self) -> int:
        return len(self.seq)


# Uni12 with its degenerate fourth base (A or G), and Uni13, each with
# its control
N_CONTROLS = 3
OLIGOS = (Oligo("Uni12", "AGCRAAAGCAGG", "12",
                ("RGGAAAGGAACC", "RGACACGGGAAA", "RCAAGGACGGAA")),
          Oligo("Uni13", "AGTAGAAACAAGG", "13",
                ("TAACAAGGGAAGA", "GGCAAAGAAATGA", "GGAAATACGGAAA")))


@dataclass(frozen=True)
class Hit:
    """A primer found in a read: its bases are seq[start:end]; `missing`
    bases of its 5' end are absent (the read begins inside it)."""
    oligo: Oligo
    start: int
    end: int
    errors: int
    missing: int


def search_left(seq: str, boundary: int, oligos: tuple[Oligo, ...],
                max_gap: int, error_rate: float,
                control: int | None = None) -> Hit | None:
    """Find a primer ending near `boundary`, the index of `seq` at which
    an alignment of the read begins.

    The primer may end up to `max_gap` bases before the boundary (template
    bases which the aligner clipped) or lie inside the alignment (primer
    bases which match the reference). With `control`, the controls of
    that number are searched instead of the primers.
    """
    best: Hit | None = None
    best_score: tuple[float, int] | None = None
    for oligo in oligos:
        query = oligo.seq if control is None else oligo.controls[control]
        size = len(query)
        first = max(0, boundary - size - max_gap)
        last = min(len(seq), boundary + size + 3)
        pad = size - MIN_PRIMER_BASES if first == 0 else 0
        target = "N" * pad + seq[first:last]
        if len(target) < size:
            continue
        result = edlib.align(query, target, mode="HW", task="locations",
                             k=int(error_rate * size),
                             additionalEqualities=_EQUALITIES)
        errors = result["editDistance"]
        if errors < 0:
            continue
        for start, end in result["locations"]:
            missing = max(0, pad - start)
            if errors > int(error_rate * (size - missing)):
                continue
            hit_end = first + end + 1 - pad
            hit_start = max(0, first + start - pad)
            if hit_end - hit_start < MIN_PRIMER_BASES - 2:
                continue
            if boundary - hit_end > max_gap:
                continue
            score = (errors / (size - missing), abs(hit_end - boundary))
            if best_score is None or score < best_score:
                best_score = score
                best = Hit(oligo, hit_start, hit_end, errors, missing)
    return best


def search_read_start(seq: str, oligos: tuple[Oligo, ...], window: int,
                      error_rate: float, control: int | None = None
                      ) -> Hit | None:
    """Find a primer within `window` bases of the start of a read which
    has no alignment (so no boundary to search at): the best match, and
    among equal ones the first."""
    best: Hit | None = None
    best_score: tuple[float, int] | None = None
    for oligo in oligos:
        query = oligo.seq if control is None else oligo.controls[control]
        size = len(query)
        pad = size - MIN_PRIMER_BASES
        target = "N" * pad + seq[:window + size]
        if len(target) < size:
            continue
        result = edlib.align(query, target, mode="HW", task="locations",
                             k=int(error_rate * size),
                             additionalEqualities=_EQUALITIES)
        errors = result["editDistance"]
        if errors < 0:
            continue
        for start, end in result["locations"]:
            missing = max(0, pad - start)
            if errors > int(error_rate * (size - missing)):
                continue
            hit_end = end + 1 - pad
            hit_start = max(0, start - pad)
            if hit_end - hit_start < MIN_PRIMER_BASES - 2:
                continue
            score = (errors / (size - missing), hit_start)
            if best_score is None or score < best_score:
                best_score = score
                best = Hit(oligo, hit_start, hit_end, errors, missing)
    return best


def contains_primer(seq: str, oligos: tuple[Oligo, ...],
                    max_errors: int = 1) -> bool:
    """Whether a primer, on either strand, lies anywhere in `seq`, with at
    most `max_errors` errors (strict: the whole sequence is searched)."""
    for oligo in oligos:
        for query in (oligo.seq, revcomp(oligo.seq)):
            if edlib.align(query, seq, mode="HW", task="distance",
                           k=max_errors,
                           additionalEqualities=_EQUALITIES
                           )["editDistance"] >= 0:
                return True
    return False


def search_before(seq: str, boundary: int, oligos: tuple[Oligo, ...],
                  min_gap: int, max_gap: int, error_rate: float
                  ) -> Hit | None:
    """Find a whole primer which ends between `min_gap` and `max_gap`
    bases before `boundary`; of several, the one nearest to it."""
    best: Hit | None = None
    first = max(0, boundary - max_gap - max(len(o.seq) for o in oligos))
    target = seq[first:max(first, boundary - min_gap)]
    for oligo in oligos:
        if len(target) < len(oligo.seq):
            continue
        result = edlib.align(oligo.seq, target, mode="HW", task="locations",
                             k=int(error_rate * len(oligo.seq)),
                             additionalEqualities=_EQUALITIES)
        if result["editDistance"] < 0:
            continue
        start, end = max(result["locations"], key=lambda loc: loc[1])
        hit = Hit(oligo, first + start, first + end + 1,
                  result["editDistance"], 0)
        if boundary - hit.end <= max_gap and (
                best is None or hit.end > best.end):
            best = hit
    return best
