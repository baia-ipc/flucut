"""The coverage which flucut shows in its HTML report.

The summary of flucut (flucut_summary.json) has no coverage, so the
profiles are read back from the plots of the report: each plot is a
stack of polygons, one per category, drawn from the depth at up to 704
positions of the reference, with the greatest depth written beside the
axis. They are what flucut predicts, as the reader of its report sees it.
"""

import re

import numpy as np

# the fill colours of the categories in the report of flucut
_COLOURS = {"#0b9e73": "terminal+", "#7fc7a9": "terminal-",
            "#e8a013": "internal+", "#2e7fe0": "internal-",
            "#9e9e9e": "no_primer", "#cccccc": "unmapped_probe",
            "#8e44ad": "deletion"}
_SECTION = re.compile(r'<h2 id="s\d+"[^>]*>(.*?)</h2>', re.S)
_PLOT = re.compile(r'<div class="coverplot">\s*<b>(.*?)</b>\s*'
                   r'(<svg.*?</svg>)', re.S)
_POLYGON = re.compile(r'<polygon points="([^"]*)" fill="(#[0-9a-f]{6})"')
_TEXT = re.compile(r"<text[^>]*>\s*([\d,]+)\s*</text>")
_PLOT_HEIGHT = 150.0


def _section(html: str, group: str | None) -> str:
    """The part of the report about a group (all of it, if the report
    has no sections of groups)."""
    sections = list(_SECTION.finditer(html))
    if not sections:
        return html
    names = [re.sub(r"<[^>]*>", "", m.group(1)).strip() for m in sections]
    matching = [i for i, name in enumerate(names)
                if group and (name == group or group in name
                              or name in group)]
    if len(matching) != 1:
        if len(sections) == 1:
            matching = [0]
        else:
            raise ValueError(f"group {group!r} not found in the report of "
                             f"flucut; groups: {', '.join(names)}")
    i = matching[0]
    end = sections[i + 1].start() if i + 1 < len(sections) else len(html)
    return html[sections[i].start():end]


def read_coverage(path: str, group: str | None
                  ) -> dict[str, dict[str, np.ndarray]]:
    """The coverage plots of a group: for each reference name, the
    positions drawn (0-based, key "positions") and the depth of each
    category at them."""
    with open(path) as handle:
        html = _section(handle.read(), group)
    out = {}
    for label, svg in _PLOT.findall(html):
        name = label.strip().rsplit(" ", 1)[-1]
        numbers = [int(t.replace(",", "")) for t in _TEXT.findall(svg)]
        if len(numbers) < 2:
            continue
        ymax, length = numbers[0], numbers[1]
        profiles: dict[str, np.ndarray] = {}
        columns = 0
        for points, colour in _POLYGON.findall(svg):
            y = np.array([float(p.split(",")[1]) for p in points.split()])
            columns = len(y) // 2
            top, base = y[:columns], y[columns:][::-1]
            profiles[_COLOURS[colour]] = (base - top) / _PLOT_HEIGHT * ymax
        if columns < 2:
            continue
        profiles["positions"] = np.floor(
            np.arange(columns) * (length - 1) / (columns - 1))
        out[name] = profiles
    return out


def on_references(coverage: dict[str, dict[str, np.ndarray]],
                  references: list, lift: dict
                  ) -> dict[str, dict[str, list[float]]]:
    """The profiles of read_coverage() at each position of the chosen
    references: a profile drawn on another reference of the same segment
    is transferred by the position map `lift` (see
    references.position_maps()), and interpolated between the positions
    drawn."""
    lengths = {r.name: r.length for r in references}
    out: dict[str, dict[str, list[float]]] = {}
    for name, profiles in coverage.items():
        positions = profiles["positions"]
        if name in lift:
            name, mapped = lift[name]
            positions = mapped[np.minimum(positions.astype(int),
                                          len(mapped) - 1)]
        if name not in lengths or name in out:
            continue
        positions = np.maximum.accumulate(positions.astype(float))
        positions += np.arange(len(positions)) * 1e-6
        at = np.arange(lengths[name])
        out[name] = {
            category: np.interp(at, positions, values, left=0,
                                right=0).round(1).tolist()
            for category, values in profiles.items()
            if category != "positions"}
    return out
