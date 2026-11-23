"""How well the predictions of flucut agree with the mapping: one check
per measure, each with its share, what it is a share of, and a state
(ok, warn, bad; na when there are too few cases to judge)."""

from dataclasses import dataclass

from .amplicons import CATEGORIES
from .html import pct
from .origins import (MEASURED_ORIGIN, ORIGINS, PREDICTED_ORIGIN,
                      depth_stats)
from .subtype import STATUS

MIN_CASES = 20   # fewer cases are not judged
# measure: (ok from, warn from); for the measures of LOWER_IS_BETTER, up to
LIMITS = {"same": (0.99, 0.95), "confirmed": (0.99, 0.95),
          "found": (0.90, 0.70), "position": (0.99, 0.95),
          "origin": (0.90, 0.75), "split": (0.95, 0.85),
          "primer_left": (0.02, 0.10), "depth_ratio": (0.10, 0.25),
          "depth_difference": (0.10, 0.25)}
LOWER_IS_BETTER = ("primer_left", "depth_ratio", "depth_difference")
END_LABELS = {"left": "5′ ends (Uni12)", "right": "3′ ends (Uni13)"}
POSITION_NEAR = ("0", "1-2", "3-5", "6-20")
ORDER = {"ok": 0, "na": 0, "warn": 1, "bad": 2}


@dataclass
class Check:
    kind: str            # the measure, a key of LIMITS
    label: str
    share: float | None  # None: nothing to compare
    cases: int           # what `share` is a share of ...
    of: str              # ... and its name
    end: str | None = None   # left, right: a measure of the read ends
    text: str = ""       # the value as written (default: the share)

    @property
    def status(self) -> str:
        return grade(self.kind, self.share, self.cases)

    @property
    def value(self) -> str:
        if self.share is None:
            return "–"
        return self.text or pct(self.share, 1, 2 if self.share > 0.9 else 1)


def grade(kind: str, share: float | None, cases: int = MIN_CASES) -> str:
    """The state of a measure from its limits."""
    if share is None or cases < MIN_CASES:
        return "na"
    good, fair = LIMITS[kind]
    if kind in LOWER_IS_BETTER:
        value = abs(share - 1) if kind == "depth_ratio" else share
        return "ok" if value <= good else "warn" if value <= fair else "bad"
    return "ok" if share >= good else "warn" if share >= fair else "bad"


def limit_text(kind: str) -> str:
    """The limits of a measure, as written in the report."""
    good, fair = LIMITS[kind]
    if kind == "depth_ratio":
        return f"within {good:.0%}, ! within {fair:.0%}"
    if kind in LOWER_IS_BETTER:
        return f"up to {good:.0%}, ! up to {fair:.0%}"
    return f"from {good:.0%}, ! from {fair:.0%}"


def worst(statuses: list[str]) -> str:
    return max(statuses, key=lambda s: ORDER[s], default="na")


def _ratio(part: float, total: float) -> float | None:
    return part / total if total else None


def origin_scores(origin: dict) -> dict:
    """Scores of the origin calls of flucut against the mapping, for
    each end, from the counts origin[end][flucut][mapping].

    `called` and `agree` are over the ends which both assign to a
    terminus or to the interior. Precision is over all the ends flucut
    calls so (`*_called`), also where the mapping finds no primer;
    recall over all the ends the mapping calls so (`*_measured`), also
    where flucut leaves the origin unknown or finds no primer."""
    scores = {}
    for end, table in origin.items():
        get = lambda said, truth: table.get(said, {}).get(truth, 0)
        tt, ti = get("terminal", "terminal"), get("terminal", "internal")
        it, ii = get("internal", "terminal"), get("internal", "internal")
        # flucut calls an origin where the mapping finds no primer
        tn, i_n = get("terminal", "no primer"), get("internal", "no primer")
        # the mapping calls an origin which flucut does not
        missed_t = get("unknown", "terminal") + get("no primer", "terminal")
        missed_i = get("unknown", "internal") + get("no primer", "internal")
        scores[end] = {
            "called": tt + ti + it + ii, "agree": tt + ii,
            "internal_called": ii + it + i_n,
            "internal_measured": ii + ti + missed_i,
            "terminal_called": tt + ti + tn,
            "terminal_measured": tt + it + missed_t,
            "internal_precision": _ratio(ii, ii + it + i_n),
            "internal_recall": _ratio(ii, ii + ti + missed_i),
            "terminal_precision": _ratio(tt, tt + ti + tn),
            "terminal_recall": _ratio(tt, tt + it + missed_t)}
    return scores


def origin_table(category: dict) -> dict[str, dict[str, int]]:
    """The amplicons compared read by read, by the origin flucut gives
    them and the origin measured: table[flucut][measured]."""
    out = {a: {b: 0 for b, _ in ORIGINS} for a, _ in ORIGINS}
    for said, row in category.items():
        for truth, n in row.items():
            if truth in CATEGORIES:
                out[PREDICTED_ORIGIN.get(said, "oth")][
                    MEASURED_ORIGIN[truth]] += n
    return out


def depth_agreement(s: dict) -> dict | None:
    """Predicted and measured depth over all segments: their ratio and
    the mean absolute difference per position, as a share of the
    measured depth."""
    pairs = [(r["length"], s["coverage"][r["name"]].get("agreement"))
             for r in s["references"]]
    pairs = [(n, a) for n, a in pairs if a]
    measured = sum(n * a["measured_mean_depth"] for n, a in pairs)
    if not measured:
        return None
    return {"ratio": sum(n * a["predicted_mean_depth"]
                         for n, a in pairs) / measured,
            "difference": sum(n * a["mean_absolute_difference"]
                              for n, a in pairs) / measured,
            "positions": sum(n for n, _ in pairs)}


def _end_checks(f: dict) -> list[Check]:
    out = []
    for end, sc in origin_scores(f.get("origin", {})).items():
        if end not in END_LABELS or not sc["called"]:
            continue
        out += [
            Check("same", "same origin", _ratio(sc["agree"], sc["called"]),
                  sc["called"], "ends both call terminal or internal", end),
            Check("confirmed", "internal calls confirmed",
                  sc["internal_precision"], sc["internal_called"],
                  "ends flucut calls internal", end),
            Check("found", "internally primed ends found",
                  sc["internal_recall"], sc["internal_measured"],
                  "ends the mapping places inside a segment", end),
            Check("confirmed", "terminal calls confirmed",
                  sc["terminal_precision"], sc["terminal_called"],
                  "ends flucut calls terminal", end),
            Check("found", "terminal ends found", sc["terminal_recall"],
                  sc["terminal_measured"],
                  "ends the mapping places at a terminus", end)]
        counts = f.get("position", {}).get(end, {})
        total = sum(counts.values())
        if total:
            near = sum(counts.get(b, 0) for b in POSITION_NEAR)
            out.append(Check("position", "position within 20 bases",
                             near / total, total, "ends located by both",
                             end))
    return out


def sample_checks(s: dict) -> list[Check]:
    """The checks of a sample (none without the reads of flucut or the
    coverage of its summary)."""
    f = s.get("flucut")
    if not f:
        return []
    out = _end_checks(f)
    trim = f.get("trim", {}).get("internal", {})
    if trim:
        total = sum(trim.values())
        out.append(Check("primer_left", "primer bases left",
                         trim.get("primer bases left", 0) / total, total,
                         "internally primed ends trimmed by flucut"))
    table = origin_table(f.get("category", {}))
    compared = sum(sum(row.values()) for row in table.values())
    if compared:
        same = sum(table[key][key] for key, _ in ORIGINS)
        out.append(Check("origin", "same read origin", same / compared,
                         compared, "amplicons compared read by read"))
    split = {tuple(k.split(",")): n for k, n in f.get("split", {}).items()}
    mapped = sum(n for (_, truth), n in split.items() if truth != "0")
    if mapped:
        same = sum(n for (said, truth), n in split.items()
                   if said == truth)
        out.append(Check("split", "same number of amplicons in the read",
                         same / mapped, mapped,
                         "reads written by flucut which map"))
    depth = depth_agreement(s)
    if depth:
        out += [Check("depth_ratio", "depth, flucut / measured",
                      depth["ratio"], depth["positions"],
                      "the depth measured, all segments",
                      text=pct(depth["ratio"], 1)),
                Check("depth_difference", "depth difference per position",
                      depth["difference"], depth["positions"],
                      "the depth measured, all segments",
                      text=pct(depth["difference"], 1))]
    return out


def _flucut_segments(s: dict) -> dict:
    """What the summary of flucut says about each segment, by segment
    number."""
    summary = (s.get("flucut") or {}).get("summary") or {}
    sample = summary.get("sample") or {}
    return {seg["segment"]: seg for seg in sample.get("segments", [])}


def depth_limits(s: dict) -> tuple[int, float]:
    """The depth and the share of the positions at that depth which
    make a segment complete (those of flucut, if its summary is
    given)."""
    summary = (s.get("flucut") or {}).get("summary") or {}
    return summary.get("min_depth", 50), summary.get("complete_fraction",
                                                     0.9)


def segment_rows(s: dict) -> list[dict]:
    """For each reference: what the measured depth says about the
    segment, what flucut says (from its summary, else from its coverage)
    and the state of their agreement."""
    min_depth, fraction = depth_limits(s)
    said = _flucut_segments(s)
    rows = []
    for ref in s["references"]:
        cover = s["coverage"][ref["name"]]
        measured = depth_stats(cover["measured"], MEASURED_ORIGIN,
                               min_depth, fraction)
        flucut = said.get(ref.get("segment"))
        if flucut is None and cover.get("predicted"):
            flucut = depth_stats(cover["predicted"], PREDICTED_ORIGIN,
                                 min_depth, fraction)
        agreement = cover.get("agreement")
        row = {"gene": ref["gene"], "reference": ref["name"],
               "measured": measured, "flucut": flucut, "ratio": None,
               "difference": None, "correlation": None, "status": "na",
               "same_state": None}
        if agreement and agreement["measured_mean_depth"]:
            mean = agreement["measured_mean_depth"]
            row["ratio"] = agreement["predicted_mean_depth"] / mean
            row["difference"] = agreement["mean_absolute_difference"] / mean
            row["correlation"] = agreement["correlation"]
        if flucut is not None:
            row["same_state"] = flucut["state"] == measured["state"]
            row["status"] = worst(
                [grade("depth_ratio", row["ratio"]),
                 grade("depth_difference", row["difference"]),
                 "ok" if row["same_state"] else "warn"])
        rows.append(row)
    return rows


def verdict(checks: list[Check], segments: list[dict],
            subtype: dict | None = None) -> tuple[str, str]:
    """The state of the agreement of a sample, and a line which says
    it; the subtype (see subtype.subtype_summary()) is one of the
    checks."""
    states = [c.status for c in checks] + \
        [row["status"] for row in segments if row["flucut"] is not None]
    if subtype:
        states.append(STATUS[subtype["status"]])
    judged = [state for state in states if state != "na"]
    if not judged:
        return "na", "nothing to compare"
    outside = sum(state != "ok" for state in judged)
    if not outside:
        return "ok", f"flucut agrees with the mapping in all " \
                     f"{len(judged)} checks"
    return worst(judged), f"{outside} of {len(judged)} checks outside " \
                          f"their limits"
