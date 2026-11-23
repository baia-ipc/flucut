"""The subtype in the report: whether the one flucut assigned is the
one the mapping shows, with the evidence for HA and NA, the subtype
IRMA assembled, and the references of all segments."""

from .html import Markup, badge, card, chip, note, table
from .subtype import PARTS, STATUS, part

MARKS = {"same": "same", "different": "different",
         "unknown": "cannot be told"}


def subtype_badge(subtype: dict | None) -> str:
    """The verdict on the subtype, as a state with its line."""
    if not subtype:
        return ""
    return badge(STATUS[subtype["status"]], subtype["text"])


def subtype_chips(subtype: dict | None) -> list[str]:
    """The subtype of the mapping and of each IRMA run."""
    if not subtype:
        return []
    out = [chip(f"mapping: {subtype['measured'] or 'no subtype'}")]
    out += [chip(f"{label}: {found or 'no subtype'}")
            for label, found in subtype["irma"].items()]
    return out


def _reference(name: str | None, named: str | None) -> str | None:
    if name is None:
        return None
    return f"{name} ({named})" if named else name


def _typed_line(row: dict, entry: dict, subtype: dict) -> list:
    """The line of HA or NA: its type by the mapping, by flucut and by
    IRMA, with the reads behind the mapping."""
    letter = PARTS[row["gene"]]
    other = entry["other"]
    return [
        (row["gene"], "hi"),
        Markup(badge(STATUS[entry["status"]], MARKS[entry["status"]])),
        (entry["type"] or f"{entry['closest'] or '–'}?", "hi"),
        _reference(row["reference"], row["reference_subtype"]),
        entry["reads"],
        f"{other['type']}: {other['reads']:,}" if other else "none",
        (entry["flucut"], "hi"),
        _reference(row["flucut_reference"],
                   row["flucut_reference_subtype"]),
    ] + [(part(found, letter), "hi") for found in subtype["irma"].values()]


def _other_line(row: dict, subtype: dict) -> list:
    """The line of another segment: the two references, as
    information."""
    same = row["same_reference"]
    mark = "" if same is None else Markup(
        badge("ok", "same reference") if same
        else badge("na", "other reference"))
    return [row["gene"], mark, None,
            _reference(row["reference"], row["reference_subtype"]), None,
            None, None,
            _reference(row["flucut_reference"],
                       row["flucut_reference_subtype"])] \
        + [None] * len(subtype["irma"])


def subtype_card(s: dict, verdict: bool = True) -> str:
    """The subtype by segment: HA and NA with their evidence, then the
    references of the other segments; `verdict`: with the line of the
    verdict above the table."""
    subtype = s.get("subtype")
    if not subtype:
        return ""
    typed = [row for row in subtype["segments"] if row["gene"] in PARTS]
    others = [row for row in subtype["segments"] if row["gene"] not in PARTS]
    lines = [_typed_line(row, subtype["genes"][row["gene"]], subtype)
             for row in typed] + [_other_line(row, subtype) for row in others]
    labels = list(subtype["irma"])
    groups = [("", 2), ("mapping", 4), ("flucut", 2)]
    if labels:
        groups.append(("IRMA", len(labels)))
    return card(
        "Subtype",
        f'<p class="chips">{subtype_badge(subtype)}</p>' if verdict else "",
        table(["segment", "", "type", "reference", "reads", "other type",
               "type", "reference"] + labels, lines, left=2, groups=groups),
        note(f"Reads: of the first {subtype['sampled_reads']:,} reads, "
             f"those which map to the references of that type; other "
             f"type: the next one. Under {subtype['min_reads']} reads on "
             f"HA or NA the subtype cannot be told. The verdict is on HA "
             f"and NA; the references of the other segments (with the "
             f"subtype of their strain) are information."))
