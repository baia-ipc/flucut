"""The whole analysis on simulated reads with a known truth."""

import os
import shutil
import subprocess
import sys

import pytest

from conftest import UNI12, UNI13_RC, random_seq
from flucut_posthoc.pipeline import Config, run
from simulation_check import simulation_check

SIMULATE = os.path.join(os.path.dirname(__file__), "..", "..", "benchmark",
                        "simulate.py")


@pytest.mark.skipif(shutil.which("minimap2") is None,
                    reason="minimap2 not in the PATH")
def test_simulated_reads(tmp_path):
    genome = tmp_path / "genome.fasta"
    with open(genome, "w") as out:
        for segment, length in ((1, 2300), (4, 1750), (8, 880)):
            seq = UNI12 + random_seq(length - 25, segment) + UNI13_RC
            out.write(f">seg{segment} test virus segment {segment}\n"
                      f"{seq}\n")
    sim = tmp_path / "sim"
    subprocess.run([sys.executable, "-I", SIMULATE, str(genome), str(sim),
                    "--reads-per-segment", "300", "--internal-share", "0.4",
                    "--seed", "3", "--error-rate", "0.02"], check=True)
    out = tmp_path / "posthoc"
    summary = run(Config(reads=[str(sim / "reads.fastq.gz")],
                         outdir=str(out), name="sim",
                         references=str(genome), threads=2))
    result = simulation_check(str(sim / "truth.tsv.gz"),
                              str(out / "reads.tsv.gz"))
    assert result["ends_correct"] / result["ends"] > 0.97
    lengths = result["lengths"]
    assert lengths.get("other", 0) / sum(lengths.values()) < 0.01
    for (kind, got), n in result["classes"].items():
        if kind != "dimer":
            assert got == "mapped"
    # the internal priming site of the simulation is the largest site
    with open(sim / "site.tsv") as handle:
        start = int(handle.readlines()[1].split("\t")[1])
    site = summary["internal_sites"][0]
    assert site["reference"] == "seg4" and site["side"] == "left"
    assert abs(site["position"] - (start + 12 + 1)) <= 1
    assert os.path.exists(out / "posthoc_report.html")
