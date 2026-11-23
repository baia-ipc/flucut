"""The segment references: loading, conserved termini, choice per sample.

The positions of the analysis are those of complete segment references
(with the untranslated regions), not of the IRMA consensus: the consensus
of IRMA starts and ends inside the segment, and can be the sequence of a
deleted genome. From a set of references of several subtypes, the one
closest to the sample is chosen for each segment by mapping a sample of
the reads to all of them. The chosen reference can still differ from the
sample by several percent, so it is then given the bases of the sample
(adapt_references()): the reads align better, and the positions stay
those of the reference.
"""

import os
import re
import subprocess
import tempfile
from dataclasses import dataclass

import edlib
import numpy as np

from .mapping import SENSITIVE, parse_record
from .seqio import fastq_records, read_fasta, write_fasta

GENES = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "M",
         8: "NS"}
TERMINUS_12 = "AGCRAAAGCAGG"
TERMINUS_13_RC = "CCTTGTTTCTACT"
_EQUALITIES = [("R", "A"), ("R", "G")]

CHOICE_READS = 20000    # reads sampled for the choice of the references
# the references which flucut downloads on first use
FLUCUT_CACHE = os.path.expanduser(
    "~/.cache/flucut/flu-a-references-v1.fasta")


@dataclass
class Reference:
    """A segment reference. `t12` and `t13` are the numbers of bases of
    the conserved termini present at its start and at its end (0: the
    record begins after the terminus, as many database records do)."""
    name: str
    description: str
    seq: str
    segment: int | None
    t12: int = 0
    t13: int = 0
    adapted: int = 0     # bases replaced by those of the sample

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def gene(self) -> str:
        return GENES.get(self.segment, self.name)

    @property
    def first_template(self) -> int:
        """0-based position of the first base after the 12 nt terminus."""
        return self.t12

    @property
    def last_template(self) -> int:
        """0-based position of the last base before the 13 nt terminus."""
        return self.length - 1 - self.t13


def _terminus_bases(seq: str) -> tuple[int, int]:
    """The bases of the conserved termini at the two ends of a segment."""
    t12 = t13 = 0
    hit = edlib.align(TERMINUS_12, seq[:20], mode="HW", task="locations",
                      k=2, additionalEqualities=_EQUALITIES)
    if hit["editDistance"] >= 0:
        t12 = hit["locations"][0][1] + 1
    tail = seq[-21:]
    hit = edlib.align(TERMINUS_13_RC, tail, mode="HW", task="locations", k=2)
    if hit["editDistance"] >= 0:
        t13 = len(tail) - hit["locations"][-1][0]
    return t12, t13


def load_references(path: str) -> list[Reference]:
    """The references of a FASTA file; the segment number is read from
    the description ("segment 4"), as flucut does."""
    references = []
    for name, description, seq in read_fasta(path):
        match = re.search(r"segment (\d)", description)
        segment = int(match.group(1)) if match else None
        t12, t13 = _terminus_bases(seq)
        references.append(Reference(name, description, seq, segment, t12,
                                    t13))
    if not references:
        raise ValueError(f"no sequences in {path}")
    return references


def write_read_sample(reads: list[str], path: str, n_reads: int) -> int:
    """Write the first `n_reads` reads of the FASTQ files to a FASTA file
    (for the choice of the references); returns the number written."""
    count = 0
    with open(path, "w") as out:
        for fastq in reads:
            for title, seq in fastq_records(fastq):
                out.write(f">r{count}\n{seq}\n")
                count += 1
                if count >= n_reads:
                    return count
    return count


def choose_references(references: list[Reference], reads: list[str],
                      minimap2: str, threads: int,
                      n_reads: int = CHOICE_READS
                      ) -> tuple[list[Reference], list[dict]]:
    """For each segment, the reference to which a sample of the reads
    maps best (most matching bases); references without a segment number
    are all kept. Returns the chosen references, ordered by segment, and
    the evidence: for each reference with reads, its matching bases."""
    if all(r.segment is None for r in references):
        return references, []
    matches = {r.name: 0 for r in references}
    n_mapped = {r.name: 0 for r in references}
    with tempfile.TemporaryDirectory() as tmp:
        ref_path = os.path.join(tmp, "references.fasta")
        sample_path = os.path.join(tmp, "sample.fasta")
        write_fasta(ref_path, [(r.name, r.seq) for r in references])
        write_read_sample(reads, sample_path, n_reads)
        proc = subprocess.run(
            [minimap2, "-x", "map-ont", "--secondary=no", "-t",
             str(threads), ref_path, sample_path],
            capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        matches[fields[5]] += int(fields[9])
        n_mapped[fields[5]] += 1
    chosen, evidence = [], []
    by_segment: dict[int | None, list[Reference]] = {}
    for ref in references:
        by_segment.setdefault(ref.segment, []).append(ref)
    for segment in sorted(by_segment, key=lambda s: (s is None, s)):
        candidates = by_segment[segment]
        if segment is None:
            chosen.extend(candidates)
            continue
        best = max(candidates, key=lambda r: matches[r.name])
        for ref in candidates:
            if matches[ref.name] > 0:
                evidence.append({"reference": ref.name, "segment": segment,
                                 "sample_reads": n_mapped[ref.name],
                                 "matching_bases": matches[ref.name],
                                 "chosen": ref is best})
        if matches[best.name] > 0:
            chosen.append(best)
    if not chosen:
        raise ValueError("no read of the sample maps to the references")
    return chosen, evidence


_CIGAR = re.compile(r"(\d+)([MIDN=X])")


def position_maps(chosen: list[Reference], references: list[Reference],
                  minimap2: str) -> dict[str, tuple[str, np.ndarray]]:
    """For each reference not chosen, of a segment with a chosen
    reference and similar enough to align to it: the name of the chosen
    reference and, for each of its positions, the corresponding position
    there (0-based). Positions which do not align (an insertion, the ends
    beyond the alignment) continue from the nearest aligned ones."""
    by_segment = {r.segment: r for r in chosen if r.segment is not None}
    others = [r for r in references
              if r.segment in by_segment and r not in chosen]
    if not others:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        target = os.path.join(tmp, "chosen.fasta")
        query = os.path.join(tmp, "others.fasta")
        write_fasta(target, [(r.name, r.seq) for r in chosen])
        write_fasta(query, [(r.name, r.seq) for r in others])
        proc = subprocess.run(
            [minimap2, "-c", *SENSITIVE, target, query],
            capture_output=True, text=True, check=True)
    segment_of = {r.name: r.segment for r in others}
    maps: dict[str, tuple[str, np.ndarray]] = {}
    best: dict[str, int] = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        name, to = f[0], by_segment[segment_of[f[0]]]
        if f[5] != to.name or f[4] != "+" or \
                int(f[9]) <= best.get(name, 0):
            continue
        best[name] = int(f[9])
        size, q, r = int(f[1]), int(f[2]), int(f[7])
        positions = np.full(size, -1, dtype=np.int64)
        cigar = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        for count, op in _CIGAR.findall(cigar):
            n = int(count)
            if op in "M=X":
                positions[q:q + n] = np.arange(r, r + n)
                q += n
                r += n
            elif op == "I":
                q += n
            else:
                r += n
        aligned = np.flatnonzero(positions >= 0)
        first, last = aligned[0], aligned[-1]
        for i in first + np.flatnonzero(positions[first:last + 1] < 0):
            positions[i] = positions[i - 1]
        # beyond the alignment: also before the first base of the chosen
        # reference (negative positions)
        positions[:first] = positions[first] - np.arange(first, 0, -1)
        positions[last + 1:] = positions[last] + \
            np.arange(1, size - last)
        maps[name] = (to.name, positions)
    return maps


# aligned bases at each end of a read left out of the adaptation
ADAPT_SKIP = 20

_BASE_INDEX = np.full(256, 4, dtype=np.int64)
for _i, _base in enumerate("ACGT"):
    _BASE_INDEX[ord(_base)] = _i


def adapt_references(chosen: list[Reference], reads: list[str],
                     minimap2: str, threads: int, n_reads: int = 200000,
                     rounds: int = 4) -> list[Reference]:
    """Give each reference the bases of the sample (see _adapt_once()),
    repeatedly: once the bases around a divergent stretch are those of
    the sample, the reads align into it, and it can be adapted in the
    next round. Stops when a round changes fewer than 3 bases."""
    with tempfile.TemporaryDirectory() as tmp:
        sample_path = os.path.join(tmp, "sample.fasta")
        write_read_sample(reads, sample_path, n_reads)
        for _ in range(rounds):
            before = sum(r.adapted for r in chosen)
            chosen = _adapt_once(chosen, sample_path, minimap2, threads)
            if sum(r.adapted for r in chosen) - before < 3:
                break
    return chosen


def _adapt_once(chosen: list[Reference], sample_path: str, minimap2: str,
                threads: int, min_depth: int = 5) -> list[Reference]:
    """Map a sample of the reads; each position covered by at least
    `min_depth` reads of which more than half agree on a base takes that
    base. Insertions and deletions are not applied, so that the positions
    remain those of the reference. The first and last ADAPT_SKIP aligned
    bases of each read are not counted: the bases under a primer are
    those of the primer, at a terminus and at an internal priming site
    alike, and where most reads start at such a site they would write
    the primer into the reference. The template there is taken from the
    reads which run through the site. The conserved termini are left as
    they are."""
    counts = {r.name: np.zeros((r.length, 5), dtype=np.int64)
              for r in chosen}
    with tempfile.TemporaryDirectory() as tmp:
        ref_path = os.path.join(tmp, "references.fasta")
        write_fasta(ref_path, [(r.name, r.seq) for r in chosen])
        proc = subprocess.Popen(
            [minimap2, "-a", *SENSITIVE, "-t", str(threads), ref_path,
             sample_path], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True)
        for line in proc.stdout:
            if line[0] == "@":
                continue
            _, aln, seq = parse_record(line)
            if aln is None:
                continue
            codes = _BASE_INDEX[np.frombuffer(seq.encode(), dtype=np.uint8)]
            table = counts[aln.ref]
            first, last = aln.qs + ADAPT_SKIP, aln.qe - ADAPT_SKIP
            for q, r, n in zip(aln.bq, aln.br, aln.bl):
                skip = max(0, first - q)
                n = min(n, last - q)
                if skip < n:
                    np.add.at(table, (np.arange(r + skip, r + n),
                                      codes[q + skip:q + n]), 1)
        proc.wait()
    adapted = []
    for ref in chosen:
        table = counts[ref.name][:, :4]
        depth = table.sum(axis=1)
        major = table.argmax(axis=1)
        sure = (depth >= min_depth) & (table.max(axis=1) * 2 > depth)
        sure[:ref.t12] = False
        if ref.t13:
            sure[-ref.t13:] = False
        bases = np.frombuffer(ref.seq.encode(), dtype=np.uint8).copy()
        new = np.frombuffer(b"ACGT", dtype=np.uint8)[major]
        changed = sure & (new != bases)
        bases[changed] = new[changed]
        adapted.append(Reference(ref.name, ref.description,
                                 bases.tobytes().decode(), ref.segment,
                                 ref.t12, ref.t13,
                                 ref.adapted + int(changed.sum())))
    return adapted
