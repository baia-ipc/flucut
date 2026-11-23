"""The plots of the report, as inline SVG: the coverage of a segment by
read origin and the lengths of its complete amplicons, drawn as in the
report of flucut."""

from html import escape

import numpy as np

from .origins import ORIGIN_LABELS, ORIGINS

COLUMNS = 250              # positions drawn along a reference
LENGTH_BARS = 62           # bars of a length plot
TOP, BOTTOM = 10.0, 104.0  # of the plot area
COVER = (560, 52.0, 552.0)   # width, left and right of the plot area
LENGTHS = (300, 44.0, 292.0)


def _columns(values, columns: int = COLUMNS) -> np.ndarray:
    """`columns` values of a profile, evenly spaced (all, if fewer)."""
    values = np.asarray(values, dtype=float)
    if len(values) <= columns:
        return values
    at = (np.arange(columns) * (len(values) - 1) / (columns - 1)).astype(int)
    return values[at]


def _y(value: float, ymax: float) -> float:
    return BOTTOM - (BOTTOM - TOP) * min(value, ymax) / ymax


def _svg(frame: tuple, body: list[str], top: str, first: str,
         last: str = "") -> str:
    width, left, right = frame
    labels = (f'<line class="axis" x1="{left}" y1="{BOTTOM}" x2="{right}" '
              f'y2="{BOTTOM}"/><text class="muted" x="{left - 4}" y="16" '
              f'text-anchor="end">{top}</text><text class="muted" '
              f'x="{left}" y="117">{first}</text>')
    if last:
        labels += (f'<text class="muted" x="{right}" y="117" '
                   f'text-anchor="end">{last}</text>')
    return (f'<svg viewBox="0 0 {width} 122" role="img">{"".join(body)}'
            f"{labels}</svg>")


def _areas(depth: dict[str, np.ndarray], ymax: float) -> list[str]:
    """The depth of each origin as stacked areas."""
    _, left, right = COVER
    body, base = [], None
    for key, _ in ORIGINS:
        if key not in depth or not np.any(depth[key]):
            continue
        values = _columns(depth[key])
        n = len(values)
        xs = left + np.arange(n) * (right - left) / max(1, n - 1)
        if base is None:
            base = np.zeros(n)
        top = base + values
        points = [f"{x:.1f},{_y(v, ymax):.1f}" for x, v in zip(xs, top)]
        points += [f"{x:.1f},{_y(v, ymax):.1f}"
                   for x, v in zip(xs[::-1], base[::-1])]
        body.append(f'<polygon class="c-{key}" points="{" ".join(points)}">'
                    f"<title>{ORIGIN_LABELS[key]}</title></polygon>")
        base = top
    return body


def _depth_lines(lines: dict[str, list], ymax: float) -> list[str]:
    """Depth profiles as lines (a value below 0: no line there)."""
    _, left, right = COVER
    body = []
    for k, (label, values) in enumerate(lines.items()):
        sampled = _columns(values)
        n = len(sampled)
        path, pen = [], False
        for i, v in enumerate(sampled):
            if v < 0:
                pen = False
                continue
            x = left + i * (right - left) / max(1, n - 1)
            path.append(f"{'L' if pen else 'M'}{x:.1f},{_y(v, ymax):.1f}")
            pen = True
        body.append(f'<path class="depth d{k % 3}" d="{" ".join(path)}">'
                    f"<title>{escape(label)}</title></path>")
    return body


def coverage_plot(depth: dict[str, np.ndarray], length: int, ymax: float,
                  min_depth: int = 0, low: list[tuple[int, int]] = (),
                  lines: dict[str, list] | None = None,
                  site: int | None = None) -> str:
    """The depth along a reference, stacked by read origin (`depth`: the
    depth of each origin at each position). `min_depth` is marked, and
    the intervals `low` below it; `lines` are other depth profiles;
    `site` is the 1-based position of an internal priming site."""
    _, left, right = COVER
    scale = (right - left) / max(1, length - 1)
    body = _areas(depth, ymax)
    for start, end in low:
        body.append(f'<rect class="c-low" x="{left + start * scale:.1f}" '
                    f'y="105" width="{max((end - start) * scale, 1):.1f}" '
                    f'height="3"><title>{start + 1}-{end}: below '
                    f"{min_depth}×</title></rect>")
    if 0 < min_depth < ymax * 0.9:
        y = _y(min_depth, ymax)
        body.append(f'<line class="mark" x1="{left}" y1="{y:.1f}" '
                    f'x2="{right}" y2="{y:.1f}"/><text class="muted" '
                    f'x="{left - 4}" y="{y + 4:.1f}" text-anchor="end">'
                    f"{min_depth}×</text>")
    body += _depth_lines(lines or {}, ymax)
    if site is not None and length > 1:
        x = left + (site - 1) * scale
        anchor, dx = ("end", -3) if x > (left + right) / 2 else ("start", 3)
        body.append(f'<line class="site" x1="{x:.1f}" y1="{TOP}" '
                    f'x2="{x:.1f}" y2="{BOTTOM}"/><text class="sitel" '
                    f'x="{x + dx:.1f}" y="18" text-anchor="{anchor}">'
                    f"{site}</text>")
    return _svg(COVER, body, f"{ymax:,.0f}×", "1", f"{length:,}")


def _bars(histogram: list[int], step: int, bars: int) -> np.ndarray:
    """A histogram with `step` of its bins in each of `bars` bars."""
    values = np.zeros(step * bars)
    part = np.asarray(histogram[:step * bars], dtype=float)
    values[:len(part)] = part
    return values.reshape(bars, step).sum(axis=1)


def length_plot(measured: list[int], said: list[int] | None, bin_size: int,
                reference_length: int, short_fraction: float = 0.9) -> str:
    """The lengths of the complete amplicons of a segment. Bars: those
    measured (short ones in another colour); line: those of flucut; the
    length of the reference is marked."""
    _, left, right = LENGTHS
    bins = int(reference_length * 1.07) // bin_size + 1
    step = -(-bins // LENGTH_BARS)
    bars = -(-bins // step)
    values = _bars(measured, step, bars)
    line = _bars(said, step, bars) if said else None
    ymax = max(values.max(), line.max() if line is not None else 0, 1)
    width = (right - left) / bars
    span = step * bin_size
    body = []
    for i, v in enumerate(values):
        if v:
            y = _y(v, ymax)
            short = (i + 1) * span <= short_fraction * reference_length
            body.append(
                f'<rect class="c-{"short" if short else "full"}" '
                f'x="{left + i * width:.1f}" y="{y:.1f}" '
                f'width="{width * 0.85:.1f}" height="{BOTTOM - y:.1f}">'
                f"<title>{i * span}-{(i + 1) * span - 1} bases: {v:,.0f}"
                f"</title></rect>")
    if line is not None and line.any():
        points = " ".join(f"{left + (i + 0.5) * width:.1f},"
                          f"{_y(v, ymax):.1f}" for i, v in enumerate(line))
        body.append(f'<polyline class="said" points="{points}">'
                    f"<title>flucut</title></polyline>")
    x = left + reference_length / span * width
    body.append(f'<line class="mark" x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" '
                f'y2="{BOTTOM}"/><text class="muted" x="{x:.1f}" y="117" '
                f'text-anchor="middle">{reference_length:,}</text>')
    return _svg(LENGTHS, body, f"{ymax:,.0f}", "0")
