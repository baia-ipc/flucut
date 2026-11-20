"""Shared helpers of the tests."""

import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from flucut_posthoc.amplicons import Params  # noqa: E402
from flucut_posthoc.mapping import Alignment  # noqa: E402
from flucut_posthoc.primers import OLIGOS  # noqa: E402
from flucut_posthoc.references import Reference  # noqa: E402
from flucut_posthoc.seqio import revcomp  # noqa: E402

UNI12 = "AGCAAAAGCAGG"
UNI13_RC = "CCTTGTTTCTACT"


def random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def params() -> Params:
    return Params(OLIGOS)


@pytest.fixture
def reference() -> Reference:
    """A segment of 600 bases with its conserved termini."""
    seq = UNI12 + random_seq(600 - 25, 1) + UNI13_RC
    return Reference("ref", "test segment 4", seq, 4, 12, 13)


def exact_alignment(ref: Reference, read: str, qs: int, rs: int, n: int,
                    reverse: bool = False) -> Alignment:
    """The alignment of read[qs:qs+n] to ref[rs:rs+n], without errors;
    for a reverse alignment `read` is given as it aligns (the reverse
    complement of the read sequenced)."""
    assert read[qs:qs + n] == ref.seq[rs:rs + n]
    return Alignment(ref.name, reverse, rs, rs + n, qs, qs + n, read, 0,
                     [qs], [rs], [n])


__all__ = ["UNI12", "UNI13_RC", "random_seq", "exact_alignment", "revcomp"]
