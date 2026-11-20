"""Mapping of the reads with minimap2 and reading of its alignments.

All coordinates of an alignment are in the orientation of the reference
(as in a SAM record): `seq` is the read as it aligns to the plus strand,
`qs` and `qe` delimit its aligned part, `rs` and `re` the reference
interval, all 0-based and half-open.
"""

import re
import subprocess
from bisect import bisect_right
from dataclasses import dataclass, field
from typing import IO, Iterator

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

# Small k-mers and windows: the references are a few kilobases, and short
# pieces of template (parts of chimeric reads, short amplicons) must map.
# Chains are not joined over gaps of more than 100 bases: joined over a
# long gap, minimap2 can report one flank of a deleted genome and drop
# the other; kept apart, both flanks are reported, and joined here (see
# amplicons.py).
SENSITIVE = ["-x", "map-ont", "-k", "12", "-w", "6", "-g", "100", "-r",
             "100,100", "--secondary=no", "-s", "30", "-m", "25"]
MINIMAP2_OPTIONS = ["-a", "-Y", *SENSITIVE]


@dataclass
class Alignment:
    ref: str
    reverse: bool          # the read is the minus strand of the reference
    rs: int
    re: int
    qs: int
    qe: int
    seq: str
    mismatches: int        # edit distance (NM), large deletions included
    # the aligned blocks, split at every insertion or deletion: parallel
    # lists of the start in the read, the start in the reference, length
    bq: list[int] = field(default_factory=list)
    br: list[int] = field(default_factory=list)
    bl: list[int] = field(default_factory=list)

    @property
    def read_start(self) -> int:
        """Start of the aligned part in the read as sequenced."""
        return len(self.seq) - self.qe if self.reverse else self.qs

    def ref_of(self, q: int, towards_end: bool) -> int:
        """The reference position aligned to the read position `q`; a
        position inside an insertion or outside the alignment takes that
        of the next aligned base (`towards_end`) or of the previous one,
        extended base by base outside the alignment."""
        if q < self.qs:
            return self.rs - (self.qs - q)
        if q >= self.qe:
            return self.re - 1 + (q - self.qe + 1)
        i = bisect_right(self.bq, q) - 1
        if q < self.bq[i] + self.bl[i]:
            return self.br[i] + (q - self.bq[i])
        if towards_end:
            return self.br[i + 1]
        return self.br[i] + self.bl[i] - 1

    def pieces(self, min_deletion: int) -> list[tuple[int, int, int, int]]:
        """The alignment as pieces (qs, qe, rs, re) separated by the
        deletions of at least `min_deletion` bases."""
        out = []
        qs, rs = self.bq[0], self.br[0]
        for i in range(1, len(self.bq)):
            prev_r = self.br[i - 1] + self.bl[i - 1]
            if self.br[i] - prev_r >= min_deletion:
                out.append((qs, self.bq[i - 1] + self.bl[i - 1], rs, prev_r))
                qs, rs = self.bq[i], self.br[i]
        out.append((qs, self.qe, rs, self.re))
        return out


def parse_record(line: str) -> tuple[str, Alignment | None, str]:
    """A SAM record: the read name, its alignment (None if unmapped) and
    the sequence of the record."""
    fields = line.split("\t", 11)
    flag = int(fields[1])
    seq = fields[9]
    if flag & 4:
        return fields[0], None, seq
    q, r = 0, int(fields[3]) - 1
    aln = Alignment(fields[2], bool(flag & 16), r, r, 0, 0, seq, 0)
    first = True
    for count, op in _CIGAR.findall(fields[5]):
        n = int(count)
        if op in "M=X":
            if first:
                aln.qs = q
                first = False
            aln.bq.append(q)
            aln.br.append(r)
            aln.bl.append(n)
            q += n
            r += n
            aln.qe = q
        elif op in "IS":
            q += n
        elif op in "DN":
            r += n
    aln.re = aln.br[-1] + aln.bl[-1]
    match = re.search(r"\tNM:i:(\d+)", "\t" + fields[11]) \
        if len(fields) > 11 else None
    aln.mismatches = int(match.group(1)) if match else 0
    return fields[0], aln, seq


def start_minimap2(minimap2: str, reference_fasta: str, reads: list[str],
                   threads: int) -> subprocess.Popen:
    """Map the reads; the SAM records are read from the standard output
    of the process returned."""
    command = [minimap2, *MINIMAP2_OPTIONS, "-t", str(threads),
               reference_fasta, *reads]
    return subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True,
                            bufsize=1 << 20)


def read_groups(sam: IO[str]) -> Iterator[tuple[str, list[str]]]:
    """The records of each read (minimap2 writes them together)."""
    name, lines = None, []
    for line in sam:
        if line[0] == "@":
            continue
        current = line[:line.index("\t")]
        if current != name:
            if lines:
                yield name, lines
            name, lines = current, []
        lines.append(line)
    if lines:
        yield name, lines
