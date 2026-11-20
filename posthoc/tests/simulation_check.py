#!/usr/bin/env python3
"""
Compare the results of flucut-posthoc on simulated reads with their truth.

The reads come from benchmark/simulate.py, which writes for each read
where each of its bases comes from (flank, primer tail, primer at a
terminus, primer at the internal site, template). The table of reads
written by flucut-posthoc (reads.tsv.gz) is checked against it: which
reads are primer dimers, what each read end is (a primer at a terminus,
a primer at the internal site, no primer) and the length of the
template.

Usage:
    simulation_check.py [options] <truth.tsv.gz> <reads.tsv.gz>

Arguments:
    <truth.tsv.gz>  Truth written by benchmark/simulate.py.
    <reads.tsv.gz>  Table of reads written by flucut-posthoc.

Options:
    -h --help     Show this help message and exit.
    -V --version  Show the script version and exit.
"""

import gzip
import sys
from collections import Counter

MIN_PRIMER_BASES = 8


def truth_ends(labels: str) -> tuple[str, str, int]:
    """What the start and the end of a simulated read are ("terminal",
    "internal" or "no primer"), and the length of its template."""
    first, last = labels.find("T"), labels.rfind("T")
    if first < 0:
        return "no template", "no template", 0

    def state(part: str) -> str:
        if part.count("i") >= MIN_PRIMER_BASES:
            return "internal"
        if part.count("p") >= MIN_PRIMER_BASES:
            return "terminal"
        return "no primer"

    return state(labels[:first]), state(labels[last + 1:]), last - first + 1


def load_truth(path: str) -> dict[str, tuple[str, str, str, int]]:
    truth = {}
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            name, _, kind, _, labels = line.rstrip("\n").split("\t")[:5]
            truth[name] = (kind, *truth_ends(labels))
    return truth


def measured_ends(fields: list[str]) -> tuple[str, str, int]:
    """The same for a row of reads.tsv.gz of a read with one amplicon."""
    def state(primer: str, origin: str) -> str:
        return origin if primer != "." else "no primer"

    left = state(fields[7], fields[9])
    right = state(fields[11], fields[13])
    if fields[5] == "-":
        left, right = right, left
    return left, right, int(fields[14])


def simulation_check(truth_path: str, reads_path: str) -> dict:
    """Agreement of the reads table with the truth; returns the counts
    and prints them."""
    truth = load_truth(truth_path)
    ends: Counter = Counter()
    lengths: Counter = Counter()
    classes: Counter = Counter()
    with gzip.open(reads_path, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            kind, start, end, length = truth[fields[0]]
            mapped = fields[4] != "."
            classes[kind, "mapped" if mapped else fields[16]] += 1
            if not mapped or fields[2] != "1":
                continue
            got_start, got_end, got_length = measured_ends(fields)
            ends[start, got_start] += 1
            ends[end, got_end] += 1
            if "no primer" not in (start, end, got_start, got_end):
                delta = abs(got_length - length)
                lengths["exact" if delta == 0 else
                        "within 2" if delta <= 2 else "other"] += 1
    states = ["terminal", "internal", "no primer"]
    print("read kind (truth)\tmeasured\treads")
    for (kind, got), n in sorted(classes.items()):
        print(f"{kind}\t{got}\t{n}")
    print("\nend (truth)\t" + "\t".join(states))
    for state in states:
        print(state + "\t" + "\t".join(str(ends[state, s]) for s in states))
    print("\ntemplate length of the reads with both primers found")
    for key in ("exact", "within 2", "other"):
        print(f"{key}\t{lengths[key]}")
    correct = sum(ends[s, s] for s in states)
    return {"ends": sum(ends.values()), "ends_correct": correct,
            "lengths": dict(lengths), "classes": dict(classes)}


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    if sys.argv[1] in ("-V", "--version"):
        sys.exit("simulation_check.py 0.1")
    simulation_check(sys.argv[1], sys.argv[2])
