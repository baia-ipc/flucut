"""Reading of FASTA and FASTQ files, and small sequence helpers."""

import gzip
import subprocess
from typing import IO, Iterator

_COMPLEMENT = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn",
                            "TGCAYRMKVHDBNtgcayrmkvhdbn")


def revcomp(seq: str) -> str:
    """The reverse complement of a sequence (IUPAC codes included)."""
    return seq.translate(_COMPLEMENT)[::-1]


def open_text(path: str) -> IO[str]:
    """Open a text file for reading, gzipped or not."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path: str) -> list[tuple[str, str, str]]:
    """The records of a FASTA file: (name, description, sequence)."""
    records: list[tuple[str, str, list[str]]] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name, _, description = line[1:].partition(" ")
                records.append((name, description, []))
            elif records:
                records[-1][2].append(line.upper())
    return [(n, d, "".join(s)) for n, d, s in records]


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    """Write (name, sequence) records to a FASTA file."""
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def fastq_lines(path: str) -> tuple[Iterator[str], subprocess.Popen | None]:
    """The lines of a FASTQ file, decompressed by gzip if needed.

    Returns the iterator and the gzip process (None for a plain file),
    which the caller waits for after reading.
    """
    if str(path).endswith(".gz"):
        proc = subprocess.Popen(["gzip", "-cd", str(path)],
                                stdout=subprocess.PIPE, text=True)
        return iter(proc.stdout), proc
    return iter(open(path)), None


def fastq_records(path: str) -> Iterator[tuple[str, str]]:
    """The (title, sequence) pairs of a FASTQ file; the title is the
    header line without its leading '@'."""
    lines, proc = fastq_lines(path)
    for i, line in enumerate(lines):
        if i % 4 == 0:
            title = line[1:].rstrip("\n")
        elif i % 4 == 1:
            yield title, line.rstrip("\n")
    if proc is not None:
        proc.wait()


def read_name(title: str) -> str:
    """The name of a read: the first word of its title."""
    for i, char in enumerate(title):
        if char in " \t":
            return title[:i]
    return title
