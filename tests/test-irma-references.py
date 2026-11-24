# (c) 2026 Giorgio Gonnella
"""Tests of scripts/irma_references.py.

Builds references from seeds and donors made up here (no network, no
IRMA installation) and checks that

- each reference begins and ends with the conserved termini of the
  segments, and holds the coding sequence of its seed unchanged;
- the untranslated regions come from a donor of the same subtype, so
  that the coding sequence starts at the right distance from the
  terminus (the length of the untranslated regions of the HA and the NA
  differs between subtypes);
- a seed of influenza B, and one of a segment without a donor, are left
  out and reported;
- the provenance table has a row per reference;
- the configuration file is written with the paths and read back.

Usage (from the repository root):

  python3 -I tests/test-irma-references.py
"""

import importlib.util
import os
import sys
import tempfile

TERMINUS_5 = "AGCAAAAGCAGG"
TERMINUS_3 = "CCTTGTTTCTACT"


def load(path):
  spec = importlib.util.spec_from_file_location("irma_references", path)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


def check(ok, what):
  if not ok:
    sys.exit("FAILED: %s" % what)


def bases(n, letter="A"):
  return letter * n


def donor(accession, segment, subtype, utr5, cds, utr3):
  """A complete record: the 5' terminus, the untranslated region, the
  coding sequence, the other untranslated region and the 3' terminus."""
  head = TERMINUS_5 + bases(utr5 - len(TERMINUS_5), "C")
  tail = bases(utr3 - len(TERMINUS_3), "G") + TERMINUS_3
  return (">%s Influenza A virus (A/duck/Nowhere/1/2020(%s)) segment %d "
          "gene, complete cds\n%s\n" %
          (accession, subtype, segment, head + cds + tail))


def seed_cds(letter, length):
  """A coding sequence: ATG, a distinctive stretch, a stop codon."""
  return "ATG" + (letter * (length - 6)) + "TAA"


def write(path, text):
  with open(path, "w") as handle:
    handle.write(text)


def read_fasta(module, path):
  return {name: (description, sequence)
          for name, description, sequence in module.read_fasta(path)}


def test_build(module, workdir):
  """The references built from seeds and donors."""
  ha_h5 = seed_cds("T", 300)
  ha_h2 = seed_cds("G", 300)
  na_n1 = seed_cds("C", 240)
  pb2 = seed_cds("A", 600)
  seeds = os.path.join(workdir, "seeds.fasta")
  write(seeds,
        ">A_HA_H5\n%s\n>A_HA_H2\n%s\n>A_NA_N1\n%s\n>A_PB2\n%s\n"
        ">B_HA\n%s\n>A_NP\n%s\n" %
        (ha_h5, ha_h2, na_n1, pb2, seed_cds("T", 300), seed_cds("A", 300)))
  # the HA donors differ in the length of their untranslated regions, as
  # the subtypes do; there is no donor of segment 5 (the NP seed)
  donors = os.path.join(workdir, "donors.fasta")
  write(donors,
        donor("D1.1", 4, "H2N9", 43, ha_h2, 41) +
        donor("D2.1", 4, "H5N1", 28, ha_h5, 44) +
        donor("D3.1", 6, "H5N1", 20, na_n1, 28) +
        donor("D4.1", 1, "H1N1", 27, pb2, 34))
  output = os.path.join(workdir, "refs.fasta")
  provenance = os.path.join(workdir, "refs.tsv")
  written, skipped = module.build(module.read_fasta(seeds),
                                 module.read_fasta(donors), output,
                                 provenance)
  check(written == 4, "four references are built, not %d" % written)
  left_out = dict(skipped)
  check("B_HA" in left_out, "the seed of influenza B is left out")
  check("A_NP" in left_out and "segment 5" in left_out["A_NP"],
        "the seed without a donor is left out: %r" % left_out)
  records = read_fasta(module, output)
  check(sorted(records) == ["A_HA_H2", "A_HA_H5", "A_NA_N1", "A_PB2"],
        "the references built: %r" % sorted(records))
  for name, (description, sequence) in records.items():
    check(sequence.startswith(TERMINUS_5),
          "%s begins with the conserved terminus" % name)
    check(sequence.endswith(TERMINUS_3),
          "%s ends with the conserved terminus" % name)
    check("segment" in description,
          "%s names its segment: %r" % (name, description))
  # the coding sequence of each seed is kept, and the untranslated
  # regions are those of the donor of its own subtype
  check(records["A_HA_H5"][1] == TERMINUS_5 + bases(28 - len(TERMINUS_5),
        "C") + ha_h5 + bases(44 - len(TERMINUS_3), "G") + TERMINUS_3,
        "the H5 reference is the H5 donor's untranslated regions around "
        "the H5 seed")
  check(records["A_HA_H2"][1].index(ha_h2) == 43,
        "the H2 coding sequence starts 43 bases inside, as in its donor")
  check("D2.1" in records["A_HA_H5"][0] and "D1.1" in records["A_HA_H2"][0],
        "each description names the donor of its untranslated regions")
  rows = open(provenance).read().strip().split("\n")
  check(len(rows) == 5, "the provenance table has a row per reference and "
        "a header, not %d lines" % len(rows))
  check(rows[0].startswith("reference\tsegment"),
        "the provenance table has a header: %r" % rows[0])


def test_config(module, workdir):
  """The paths written to the configuration file and read back."""
  os.environ["XDG_CONFIG_HOME"] = workdir
  module.write_config({"irma": "/opt/IRMA", "irma_module": "FLU",
                       "references": "/tmp/refs.fasta"})
  config = module.read_config()
  check(config.get("irma") == "/opt/IRMA" and
        config.get("references") == "/tmp/refs.fasta",
        "the paths are read back: %r" % config)
  module.write_config({"references": "/tmp/other.fasta"})
  config = module.read_config()
  check(config.get("references") == "/tmp/other.fasta" and
        config.get("irma") == "/opt/IRMA",
        "a new value replaces one key and keeps the others: %r" % config)


def main():
  here = os.path.dirname(os.path.abspath(__file__))
  module = load(os.path.join(here, "..", "scripts", "irma_references.py"))
  module.log = lambda message: None
  with tempfile.TemporaryDirectory(prefix="fc-irma-refs-") as workdir:
    test_build(module, workdir)
    test_config(module, workdir)
  print("[OK] irma_references")


if __name__ == "__main__":
  main()
