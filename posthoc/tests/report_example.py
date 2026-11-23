"""A small summary of a sample, as flucut-posthoc writes it, for the
tests of the report."""

import copy

from flucut_posthoc.subtype import subtype_summary

LENGTH = 400


def _coverage(scale: float) -> dict:
    measured = {"terminal+": [100] * LENGTH, "terminal-": [10] * LENGTH,
                "internal+": [0] * 200 + [40] * 200,
                "deletion": [30] * 100 + [0] * 200 + [30] * 100}
    predicted = {"terminal+": [100 * scale] * LENGTH,
                 "terminal-": [10 * scale] * LENGTH,
                 "internal+": [0] * 200 + [40 * scale] * 200,
                 "deletion": [30 * scale] * 100 + [0] * 200
                 + [30 * scale] * 100, "not_located": [1.0] * LENGTH}
    total = [sum(v) for v in zip(*measured.values())]
    said = [sum(v) for v in zip(*predicted.values())]
    return {"measured": measured, "predicted": predicted, "irma": {},
            "agreement": {
                "measured_mean_depth": sum(total) / LENGTH,
                "predicted_mean_depth": sum(said) / LENGTH,
                "correlation": 0.99,
                "mean_absolute_difference": sum(
                    abs(a - b) for a, b in zip(total, said)) / LENGTH}}


def example_subtype(flucut: str | None = "H5N1", ha_reads: int = 500,
                    other_h: int = 0) -> dict:
    """The subtype of a sample whose reads map to an H5 and to an N1
    reference; `other_h`: reads on an H7 reference."""
    evidence = [
        {"reference": "ha5", "segment": 4, "sample_reads": ha_reads,
         "matching_bases": 1000 * ha_reads, "chosen": True},
        {"reference": "ha7", "segment": 4, "sample_reads": other_h,
         "matching_bases": 900 * other_h, "chosen": False},
        {"reference": "na1", "segment": 6, "sample_reads": 300,
         "matching_bases": 300000, "chosen": True},
        {"reference": "pb2", "segment": 1, "sample_reads": 900,
         "matching_bases": 900000, "chosen": True}]
    subtypes = {"ha5": "H5N1", "ha7": "H7N9", "na1": "H5N1", "pb2": "H7N9",
                "pb2b": "H5N1"}
    references = [{"name": "pb2", "gene": "PB2", "segment": 1},
                  {"name": "ha5", "gene": "HA", "segment": 4},
                  {"name": "na1", "gene": "NA", "segment": 6}]
    sample = {"subtype": flucut, "segments": [
        {"segment": 1, "reference": "pb2b"},
        {"segment": 4, "reference": "ha5"},
        {"segment": 6, "reference": "na1"}]}
    return subtype_summary(evidence, subtypes, references, sample,
                           {"IRMA": [["4-A_HA_H5", 400], ["4-A_NA_N1", 250],
                                     ["5-A_HA_H2", 3]]}, 20000)


def example_summary(scale: float = 1.0, name: str = "one") -> dict:
    """A sample with one segment (HA); `scale`: the depth of flucut as a
    share of the depth measured."""
    origin = {"terminal": {"terminal": 900, "internal": 2},
              "internal": {"internal": 95},
              "unknown": {"terminal": 30, "internal": 5}}
    return {
        "version": "0.1", "updated": "2026-01-01T00:00:00+0000",
        "sample": name,
        "inputs": {"reads": ["reads.fastq.gz"], "irma": [],
                   "flucut_reads": "trimmed.fastq.gz",
                   "flucut_summary": "flucut_summary.json",
                   "references": "refs.fasta"},
        "settings": {"primers": [{"name": "Uni12", "sequence": "AGC"}],
                     "error_rate": 0.2, "terminal_window": 10,
                     "min_deletion": 100, "length_bin": 10,
                     "short_fraction": 0.9},
        "references": [{
            "name": "r", "gene": "HA", "segment": 4, "length": LENGTH,
            "terminus_12_bases": 12, "terminus_13_bases": 13,
            "adapted_bases": 3, "amplicons": 1100, "identity": 0.98,
            "irma": {}}],
        "reads": {"total": 1200, "bases": 480000,
                  "classes": {"influenza": 1100, "offtarget": 80,
                              "dimer": 20},
                  "amplicons_per_read": {"1": 1090, "2": 10},
                  "junctions": {"amplicons joined end to end": 10},
                  "joined_not_influenza": 4},
        "amplicons": {"total": 1110,
                      "categories": {"terminal+": 800, "terminal-": 80,
                                     "internal+": 100, "deletion": 130},
                      "primers": {"2": 1000, "1": 100, "0": 10},
                      "strand": {"+": 1030, "-": 80}},
        "ends": {"left": {"12": {"terminal": 960, "internal": 100},
                          "none": {"internal": 50}},
                 "right": {"13": {"terminal": 1000, "internal": 2}}},
        "primer_search": {"errors": {"Uni12": {"0": 900, "1": 100}},
                          "truncated": {"Uni12": 3},
                          "before_primer": {"Uni12": [["ACGT", 700],
                                                      ["all", 1000]]},
                          "control": {"at an alignment boundary": {
                              "searched": 1000, "found_each": [1, 2, 4],
                              "found": 2}}},
        "terminal_offsets": {"left": {"0": 950, "5": 10}, "right": {}},
        "internal_sites": [{
            "gene": "HA", "reference": "r", "side": "left", "primer": "12",
            "position": 201, "from": 200, "to": 202, "reads": 95,
            "template": "AGCTATAGCAGG", "mismatches": 2,
            "flucut_reads": 90}],
        "internal_sites_total": 1,
        "fragment_ends": [{"reference": "r", "side": "left",
                           "position": 77, "reads": 20}],
        "deletions": [{"reference": "r", "from": 100, "to": 300,
                       "reads": 120}],
        "complete_amplicons": [{
            "reference": "r", "gene": "HA", "reference_length": LENGTH,
            "amplicons": 1000, "main_length": 390, "short": 150,
            "short_main_length": 190,
            "histogram": [0] * 19 + [150] + [0] * 19 + [850],
            "flucut": {"amplicons": 990, "main_length": 390, "short": 148,
                       "short_main_length": 190,
                       "histogram": [0] * 19 + [148] + [0] * 19 + [842]}}],
        "subtype": example_subtype(),
        "coverage": {"r": _coverage(scale)},
        "predicted_coverage_from": "summary",
        "offtarget": [{"share_of_sample": 0.6, "template_length": 150,
                       "representative": "ACGT" * 30,
                       "sampled_reads": 6}],
        "irma": [],
        "flucut": {
            "reads_are_flucut_output": False,
            "written": {"written": {"influenza": 1100},
                        "not written": {"unmapped: dimer": 20}},
            "split": {"1,1": 1080, "2,2": 10, "1,0": 80, "1,2": 2},
            "origin": {"left": origin, "right": copy.deepcopy(origin)},
            "scores": {},
            "position": {"left": {"0": 990, "1-2": 5, ">20": 1}},
            "category": {"terminal+": {"terminal+": 790, "deletion": 20},
                         "internal+": {"internal+": 95},
                         "deletion": {"deletion": 105},
                         "not_located": {"terminal+": 10,
                                         "unmapped: offtarget": 80}},
            "trim": {"internal": {"whole primer removed": 99,
                                  "primer bases left": 1},
                     "terminal": {"terminus kept": 900}},
            "primer_bases_left": {"12": 1},
            "sites_only_flucut": [],
            "predicted_on_other_references": 0,
            "summary": {
                "reads": {"in": 1200, "both_primers": 1080,
                          "one_primer": 100, "no_primer": 30, "dimers": 18,
                          "split": 12, "too_short": 25},
                "ends": {"left": {"terminal": 930, "internal": 96}},
                "short_complete": 148,
                "written_by_category": {"full": 880, "internal": 96,
                                        "deletion": 125, "other": 90},
                "read_lengths": {"input": {"median": 400, "n50": 410}},
                "read_quality": {"input": {"median": 15.5}},
                "sequencing": {"flow_cell": "FAX00001"},
                "min_depth": 50, "complete_fraction": 0.9,
                "sample": {"segments": [{
                    "segment": 4, "gene": "HA", "state": "complete",
                    "median_depth": 141 * scale, "mean_depth": 146 * scale,
                    "covered": 1.0, "well_covered": 1.0, "internal": 0.14,
                    "deletion": 0.1, "minus_strand": 0.07}]}}},
    }
