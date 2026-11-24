"""The first part of the report: whether the predictions of flucut agree
with the mapping, by measure and by segment."""

from html import escape

from .agreement import (END_LABELS, Check, limit_text, segment_rows,
                        verdict)
from .html import (Markup, badge, card, chip, chips, note, pct, section,
                   table, tile, tiles)
from .origins import STATE_LABELS, STATE_STATUS
from .report_subtype import subtype_badge, subtype_card, subtype_chips

MAX_CHIPS = 6


def _cell(check: Check | None) -> Markup:
    """The value of a check with its state, and what it is a share of."""
    if check is None or check.share is None:
        return Markup("–")
    return Markup(f"{badge(check.status, check.value)} "
                  f"<small>of {check.cases:,}</small>")


def _label(check: Check, ends: bool = False) -> Markup:
    of = ("of the " if ends else "of ") + check.of
    return Markup(f"{escape(check.label)} <small>{escape(of)}</small>")


def _verdict(checks: list[Check], rows: list[dict],
             subtype: dict | None) -> str:
    status, text = verdict(checks, rows, subtype)
    outside = [chip("subtype differs")] \
        if subtype and subtype["status"] == "different" else []
    outside += [chip(f"{END_LABELS[c.end].split(' (')[0] + ': ' if c.end else ''}"
                    f"{c.label} {c.value}")
               for c in checks if c.status in ("warn", "bad")]
    outside += [chip(f"{row['gene']}: "
                     + ("state differs" if not row["same_state"]
                        else "depth differs"))
                for row in rows if row["status"] in ("warn", "bad")]
    if len(outside) > MAX_CHIPS:
        outside = outside[:MAX_CHIPS] + [chip(
            f"{len(outside) - MAX_CHIPS} more below")]
    return chips(badge(status, text), *outside)


def _pick(checks: list[Check], kind: str, label: str = "",
          end: str | None = None) -> Check | None:
    return next((c for c in checks if c.kind == kind and c.end == end
                 and (not label or c.label == label)), None)


def _tiles(checks: list[Check], rows: list[dict]) -> str:
    out = []
    same = [c for c in checks if c.kind == "same"]
    if same:
        cases = sum(c.cases for c in same)
        pooled = Check("same", "", sum(c.share * c.cases for c in same)
                       / cases, cases, "")
        out.append(tile(pooled.value, "end origins the same",
                        f"of {cases:,} ends both call terminal or internal",
                        pooled.status))
    for check, label in (
            (_pick(checks, "found", "internally primed ends found", "left"),
             "internally primed 5′ ends found"),
            (_pick(checks, "primer_left"), "primer bases left"),
            (_pick(checks, "origin"), "same read origin")):
        if check and check.share is not None:
            out.append(tile(check.value, label,
                            f"of {check.cases:,} {check.of}", check.status))
    ratio = _pick(checks, "depth_ratio")
    difference = _pick(checks, "depth_difference")
    if ratio:
        out.append(tile(ratio.value, "depth, flucut / measured",
                        f"{difference.value} difference per position",
                        ratio.status))
    compared = [row for row in rows if row["same_state"] is not None]
    if compared:
        same_state = sum(row["same_state"] for row in compared)
        out.append(tile(f"{same_state}/{len(compared)}",
                        "segments in the same state",
                        "complete, partial, low depth or not found",
                        "ok" if same_state == len(compared) else "warn"))
    return tiles(*out)


def _measures(checks: list[Check]) -> str:
    """The checks: those of the read ends with a column for each end,
    then the others."""
    parts = []
    ends = [end for end in END_LABELS if any(c.end == end for c in checks)]
    if ends:
        by_label: dict[str, dict[str, Check]] = {}
        for c in checks:
            if c.end:
                by_label.setdefault(c.label, {})[c.end] = c
        lines = []
        for found in by_label.values():
            first = next(iter(found.values()))
            lines.append([_label(first, ends=True)]
                         + [_cell(found.get(end)) for end in ends]
                         + ["✓ " + limit_text(first.kind)])
        parts.append(table(["read ends"] + [END_LABELS[e] for e in ends]
                           + ["limits"], lines))
    others = [c for c in checks if not c.end]
    if others:
        parts.append(table(
            ["reads and depth", "", "limits"],
            [[_label(c), _cell(c), "✓ " + limit_text(c.kind)]
             for c in others]))
    return card("Agreement by measure", *parts)


def _state(stats: dict) -> Markup:
    return Markup(badge(STATE_STATUS[stats["state"]],
                        STATE_LABELS[stats["state"]]))


def _segments(rows: list[dict]) -> str:
    lines = []
    for row in rows:
        m, f = row["measured"], row["flucut"]
        if f is None:
            continue
        r = row["correlation"]
        lines.append([
            Markup(f'<a href="#seg-{escape(row["gene"])}">'
                   f'{escape(row["gene"])}</a>'),
            Markup(badge(row["status"], "agrees" if row["status"] == "ok"
                         else "not compared" if row["status"] == "na"
                         else "differs")),
            _state(m), _state(f),
            f"{m['median_depth']:,.0f}×", f"{f['median_depth']:,.0f}×"
            if f.get("median_depth") is not None else None,
            f"{m['mean_depth']:,.0f}×", f"{f['mean_depth']:,.0f}×"
            if f.get("mean_depth") is not None else None,
            pct(row["ratio"], 1) if row["ratio"] is not None else None,
            pct(row["difference"], 1) if row["difference"] is not None
            else None,
            f"{r:.3f}" if r is not None else None])
    if not lines:
        return ""
    return card("Agreement by segment", table(
        ["segment", "", "measured", "flucut", "measured", "flucut",
         "measured", "flucut", "flucut / measured",
         "difference per position", "correlation"], lines, left=2,
        groups=[("", 2), ("state", 2), ("median depth", 2),
                ("mean depth", 2), ("depth along the segment", 3)]))


def agreement_section(s: dict, checks: list[Check]) -> str:
    """The section of the report on the agreement of flucut with the
    mapping (empty without anything of flucut to compare)."""
    rows = segment_rows(s)
    if not checks and all(row["flucut"] is None for row in rows):
        return ""
    subtype = s.get("subtype")
    parts = [chips(subtype_badge(subtype), *subtype_chips(subtype)),
             _verdict(checks, rows, subtype), _tiles(checks, rows)]
    if s["flucut"].get("reads_are_flucut_output"):
        parts.append(note(
            "The reads mapped are those written by flucut, without their "
            "primers: an end is judged by its position alone, and the "
            "bases removed cannot be checked."))
    elif not s["flucut"].get("origin"):
        parts.append(note(
            "Without the reads written by flucut, only its summary is "
            "compared: the depth and the amplicon lengths."))
    return section("agreement", "flucut and the mapping", *parts,
                   subtype_card(s, verdict=False), _measures(checks),
                   _segments(rows))
