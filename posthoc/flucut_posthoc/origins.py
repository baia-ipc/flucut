"""The origin of the reads, in the four classes of the report of flucut,
and what the depth by origin says about a segment."""

import numpy as np

# (key, label), in the order of stacking of the coverage plots
ORIGINS = [("oth", "other"), ("full", "full template"),
           ("int", "internally primed"), ("del", "internal deletion")]
ORIGIN_LABELS = dict(ORIGINS)
# the categories of the amplicons measured, by origin
MEASURED_ORIGIN = {"terminal+": "full", "terminal-": "full",
                   "internal+": "int", "internal-": "int",
                   "internal3": "int", "deletion": "del",
                   "noprimer_terminal": "oth", "noprimer_fragment": "oth"}
# the categories of the reads of flucut, by origin
PREDICTED_ORIGIN = {"terminal+": "full", "terminal-": "full",
                    "internal+": "int", "internal-": "int",
                    "deletion": "del", "no_primer": "oth",
                    "not_located": "oth", "unmapped_probe": "oth"}
MINUS = ("terminal-", "internal-")
STRANDED = ("terminal+", "terminal-", "internal+", "internal-")
STATE_LABELS = {"complete": "complete", "partial": "partial",
                "low": "low depth", "missing": "not found"}
STATE_STATUS = {"complete": "ok", "partial": "warn", "low": "bad",
                "missing": "bad"}
LOW_FRACTION = 0.5    # below this share of positions deep enough: low


def by_origin(profiles: dict[str, list], classes: dict[str, str]
              ) -> dict[str, np.ndarray]:
    """The depth of each origin along a reference, from the depth of
    each category (`classes`: the origin of each category)."""
    out: dict[str, np.ndarray] = {}
    for category, values in profiles.items():
        origin = classes.get(category, "oth")
        depth = np.asarray(values, dtype=float)
        out[origin] = out[origin] + depth if origin in out else depth
    return out


def depth_stats(profiles: dict[str, list], classes: dict[str, str],
                min_depth: int = 50, complete_fraction: float = 0.9
                ) -> dict:
    """What the depth by category says about a segment, as flucut
    summarises its own estimate: the median and mean depth, the shares
    of the positions with a read and with at least `min_depth` reads,
    the shares of the depth from internally primed reads, from reads
    with an internal deletion and from the minus strand, and the state
    of the segment."""
    depth = by_origin(profiles, classes)
    total = sum(depth.values())
    bases = float(total.sum())
    if bases == 0:
        return {"state": "missing", "median_depth": 0, "mean_depth": 0.0,
                "covered": 0.0, "well_covered": 0.0, "internal": None,
                "deletion": None, "minus_strand": None}
    well = float((total >= min_depth).mean())
    stranded = sum(float(np.sum(profiles.get(c, 0))) for c in STRANDED)
    minus = sum(float(np.sum(profiles.get(c, 0))) for c in MINUS)
    return {
        "state": "complete" if well >= complete_fraction else
        "low" if well < LOW_FRACTION else "partial",
        "median_depth": float(np.sort(total)[len(total) // 2]),
        "mean_depth": float(total.mean()),
        "covered": float((total > 0).mean()), "well_covered": well,
        "internal": float(depth["int"].sum()) / bases
        if "int" in depth else 0.0,
        "deletion": float(depth["del"].sum()) / bases
        if "del" in depth else 0.0,
        "minus_strand": minus / stranded if stranded else None}


def low_intervals(profiles: dict[str, list], min_depth: int
                  ) -> list[tuple[int, int]]:
    """The intervals of a reference (0-based, end excluded) with fewer
    than `min_depth` reads."""
    total = sum(np.asarray(v, dtype=float) for v in profiles.values())
    low = np.concatenate(([0], (total < min_depth).astype(int), [0]))
    edges = np.flatnonzero(np.diff(low))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2])]


def origin_counts(counts: dict[str, int], classes: dict[str, str]
                  ) -> dict[str, int]:
    """Counts by category as counts by origin."""
    out = {key: 0 for key, _ in ORIGINS}
    for category, n in counts.items():
        out[classes.get(category, "oth")] += n
    return out
