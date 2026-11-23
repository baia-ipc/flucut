"""What flucut reports for a sample in its JSON summary
(flucut_summary.json): the section of the sample, and the coverage by
read origin which flucut estimates without mapping.

The coverage of a reference is given at up to 200 positions, as the
depth of each class of reads; these are what flucut predicts, as the
reader of its report sees it.
"""

import json

import numpy as np

# the classes of the coverage of flucut, by their names in its summary
COVERAGE_CLASSES = {"terminal_plus": "terminal+",
                    "terminal_minus": "terminal-",
                    "internal_plus": "internal+",
                    "internal_minus": "internal-",
                    "no_primer": "no_primer", "not_located": "not_located",
                    "deletion": "deletion"}
# what the summary says about the sample and its reads, shown next to
# the values measured
SAMPLE_MEMBERS = ("reads", "ends", "short_complete", "written_by_category",
                  "sample", "read_lengths", "read_quality", "sequencing")


def load_group(path: str, group: str | None, name: str) -> dict:
    """The section of a sample in a JSON summary of flucut: the group
    named, or the one whose name the sample name contains, or all reads
    if the run had one group. The depth and the share of the positions
    which make a segment complete are added to it (`min_depth`,
    `complete_fraction`)."""
    with open(path) as handle:
        data = json.load(handle)
    section = _find_group(data, group, name, path)
    for key in ("min_depth", "complete_fraction"):
        if key in data:
            section[key] = data[key]
    return section


def _find_group(data: dict, group: str | None, name: str, path: str) -> dict:
    groups = [g for g in data.get("groups", []) if g.get("processed")]
    for g in groups:
        if g["name"] == (group or name):
            return g
    if group is None:
        matching = [g for g in groups if g["name"] in name
                    or name in g["name"]]
        if len(matching) == 1:
            return matching[0]
        if len(groups) <= 1:
            return data.get("all")
    raise ValueError(f"group {group or name!r} not found in {path}; "
                     f"groups: {', '.join(g['name'] for g in groups)}")


def sample_members(section: dict) -> dict:
    """The members of a section which describe the sample (None for
    those the summary does not have)."""
    out = {key: section.get(key) for key in SAMPLE_MEMBERS}
    out["min_depth"] = section.get("min_depth", 50)
    out["complete_fraction"] = section.get("complete_fraction", 0.9)
    return out


def read_coverage(section: dict) -> dict[str, dict[str, np.ndarray]]:
    """The coverage of a section: for each reference name, the positions
    given (0-based, key "positions") and the depth of each class at
    them."""
    out = {}
    for entry in section.get("coverage") or []:
        positions = np.array(entry["positions"], dtype=float)
        if len(positions) < 2:
            continue
        profiles = {COVERAGE_CLASSES[key]: np.array(depth, dtype=float)
                    for key, depth in entry["depth"].items()
                    if key in COVERAGE_CLASSES}
        profiles["positions"] = positions
        out[entry["reference"]] = profiles
    return out


def on_references(coverage: dict[str, dict[str, np.ndarray]],
                  references: list, lift: dict
                  ) -> dict[str, dict[str, list[float]]]:
    """The profiles of read_coverage() at each position of the chosen
    references: a profile given on another reference of the same segment
    is transferred by the position map `lift` (see
    references.position_maps()), and interpolated between the positions
    given. The profiles of several references of flucut which land on
    one chosen reference are added (flucut can locate the reads of a
    segment in more than one of its references); a reference without
    depth tells nothing and is left out."""
    lengths = {r.name: r.length for r in references}
    sums: dict[str, dict[str, np.ndarray]] = {}
    for name, profiles in coverage.items():
        categories = [c for c in profiles if c != "positions"]
        if not categories:
            continue
        positions = profiles["positions"]
        if name in lift:
            name, mapped = lift[name]
            positions = mapped[np.minimum(positions.astype(int),
                                          len(mapped) - 1)]
        if name not in lengths:
            continue
        positions = np.maximum.accumulate(positions.astype(float))
        positions += np.arange(len(positions)) * 1e-6
        at = np.arange(lengths[name])
        on_reference = sums.setdefault(name, {})
        for category in categories:
            depth = np.interp(at, positions, profiles[category], left=0,
                              right=0)
            on_reference[category] = on_reference.get(category, 0) + depth
    return {name: {category: depth.round(1).tolist()
                   for category, depth in profiles.items()}
            for name, profiles in sums.items()}
