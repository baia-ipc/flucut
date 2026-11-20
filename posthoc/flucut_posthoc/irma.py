"""Reading of the results of an IRMA run.

Used from an IRMA output directory: the read counts of its steps, the
consensus of each segment, the reads it assembled (the names in its BAM
files) and the depth it reports. The depth is given along the consensus
of IRMA; it is transferred to the positions of the segment reference by
aligning the consensus to it.
"""

import glob
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field

import edlib
import numpy as np

from .references import Reference
from .mapping import SENSITIVE
from .primers import OLIGOS
from .seqio import read_fasta, revcomp, write_fasta

_CIGAR = re.compile(r"(\d+)([MIDN=X])")
_PART = re.compile(r"_\d+$")


def gene_of(irma_name: str) -> str:
    """The gene of an IRMA assembly name (A_HA_H5: HA; A_MP: M)."""
    gene = irma_name.split("_")[1] if "_" in irma_name else irma_name
    return "M" if gene == "MP" else gene


@dataclass
class IrmaRun:
    label: str
    path: str
    read_counts: list[tuple[str, int]] = field(default_factory=list)
    consensus: dict[str, str] = field(default_factory=dict)   # by gene
    depth: dict[str, np.ndarray] = field(default_factory=dict)
    used: dict[str, str] = field(default_factory=dict)  # read: gene
    module: str = ""


def _read_depth(path: str) -> np.ndarray:
    values = []
    with open(path) as handle:
        next(handle)
        for line in handle:
            values.append(int(line.split("\t")[2]))
    return np.array(values, dtype=np.int64)


def load_irma(path: str, samtools: str, label: str | None = None
              ) -> IrmaRun:
    path = path.rstrip("/")
    run = IrmaRun(label or os.path.basename(path), path)
    counts = os.path.join(path, "tables", "READ_COUNTS.txt")
    if os.path.exists(counts):
        with open(counts) as handle:
            next(handle)
            for line in handle:
                fields = line.split("\t")
                run.read_counts.append((fields[0], int(fields[1])))
    info = os.path.join(path, "logs", "run_info.txt")
    if os.path.exists(info):
        for line in open(info):
            if line.startswith("module_param_call"):
                run.module = line.rstrip("\n").split("\t")[2]
    for fasta in sorted(glob.glob(os.path.join(path, "*.fasta"))):
        name = os.path.basename(fasta)[:-len(".fasta")]
        gene = gene_of(name)
        run.consensus[gene] = read_fasta(fasta)[0][2]
        table = os.path.join(path, "tables", name + "-coverage.txt")
        if os.path.exists(table):
            run.depth[gene] = _read_depth(table)
        bam = os.path.join(path, name + ".bam")
        if os.path.exists(bam):
            _read_names(bam, samtools, gene, run.used)
    if not run.consensus:
        raise ValueError(f"no consensus (*.fasta) in {path}: not the "
                         "output directory of an IRMA run?")
    return run


def _read_names(bam: str, samtools: str, gene: str,
                used: dict[str, str]) -> None:
    """The reads of a BAM file of IRMA. The parts of a read which flucut
    split (<read>_1, <read>_2) are also entered under the name of the
    read."""
    proc = subprocess.Popen([samtools, "view", bam], stdout=subprocess.PIPE,
                            text=True)
    for line in proc.stdout:
        name = line[:line.index("\t")]
        used[name] = gene
        if _PART.search(name):
            used.setdefault(_PART.sub("", name), gene)
    proc.wait()


def _deletions(pieces: list[tuple[int, int, int, int]]) -> list[list[int]]:
    """The deletions of a consensus from its aligned pieces (consensus
    from, to, reference from, to): between consecutive pieces, at least
    100 bases of the reference more than of the consensus; 1-based,
    inclusive."""
    out = []
    pieces = sorted(pieces)
    for (_, qe, _, re), (qs, _, rs, _) in zip(pieces, pieces[1:]):
        if (rs - re) - max(0, qs - qe) >= 100:
            out.append([re + 1, rs])
    return out


def _primers_inside(seq: str) -> int:
    """The universal primers found more than 30 bases from both ends of
    a consensus, where no primer of a single amplicon can be: the
    consensus was built from amplicons joined end to end."""
    count = 0
    for oligo in OLIGOS:
        for query in (oligo.seq, revcomp(oligo.seq)):
            hits = edlib.align(query, seq, mode="HW", task="locations", k=1,
                               additionalEqualities=[("R", "A"), ("R", "G"),
                                                     ("Y", "C"), ("Y", "T")])
            ends = {end for _, end in hits["locations"] or []}
            count += sum(1 for end in ends
                         if 30 < end < len(seq) - 30 and end - 1 not in ends)
    return count


def lift_to_references(run: IrmaRun, references: list[Reference],
                       minimap2: str) -> dict[str, dict]:
    """Align each consensus of IRMA to the reference of its segment. For
    each reference name: the length of the consensus, the bases of the
    reference it covers, its deletions of 100 bases or more with respect
    to the reference (1-based, inclusive), and the depth reported by IRMA
    at each position of the reference (-1 where the consensus has no
    base), and the number of primers inside the consensus (see
    _primers_inside())."""
    by_gene = {r.gene: r for r in references}
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        for gene, seq in run.consensus.items():
            ref = by_gene.get(gene)
            if ref is None:
                continue
            ref_path = os.path.join(tmp, "ref.fasta")
            query_path = os.path.join(tmp, "consensus.fasta")
            write_fasta(ref_path, [(ref.name, ref.seq)])
            write_fasta(query_path, [(gene, seq)])
            proc = subprocess.run(
                [minimap2, "-c", *SENSITIVE, ref_path, query_path],
                capture_output=True, text=True, check=True)
            depth = np.full(ref.length, -1, dtype=np.int64)
            reported = run.depth.get(gene)
            pieces = []
            for line in proc.stdout.splitlines():
                fields = line.split("\t")
                if fields[4] != "+":
                    continue
                q, r = int(fields[2]), int(fields[7])
                pieces.append((q, int(fields[3]), r, int(fields[8])))
                cigar = next(f[5:] for f in fields[12:]
                             if f.startswith("cg:Z:"))
                for count, op in _CIGAR.findall(cigar):
                    n = int(count)
                    if op in "M=X":
                        if reported is not None and q + n <= len(reported):
                            depth[r:r + n] = reported[q:q + n]
                        else:
                            depth[r:r + n] = 0
                        q += n
                        r += n
                    elif op == "I":
                        q += n
                    else:
                        r += n
            covered = np.flatnonzero(depth >= 0)
            out[ref.name] = {
                "consensus_length": len(seq),
                "reference_covered": int(len(covered)),
                "deletions": _deletions(pieces),
                "primers_inside": _primers_inside(seq),
                "depth": depth.tolist() if reported is not None else None}
    return out
