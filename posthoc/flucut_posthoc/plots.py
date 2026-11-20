"""The plots of the report, as inline SVG."""

from html import escape

from .html import pct

# (key, label, colour), in the order of stacking from the bottom: the
# categories which cover a segment evenly first, so that a step in the
# bands above them is that of the band itself
MEASURED_BANDS = [
    ("noprimer_terminal", "no primer at the 5′ end, starts at the terminus",
     "#9e9e9e"),
    ("terminal+", "terminal +", "#0b9e73"),
    ("terminal-", "terminal −", "#7fc7a9"),
    ("internal3", "internally primed at the 3′ end", "#d6589f"),
    ("internal+", "internally primed +", "#e8a013"),
    ("internal-", "internally primed −", "#2e7fe0"),
    # fragments start anywhere: not flat, so above the bands of the
    # amplicons with a primer
    ("noprimer_fragment", "no primer at the 5′ end, starts inside",
     "#cccccc"),
    ("deletion", "internal deletion", "#8e44ad"),
]
PREDICTED_BANDS = [
    ("unmapped_probe", "5′ probe not located", "#cccccc"),
    ("no_primer", "no 5′ primer", "#9e9e9e"),
    ("terminal+", "terminal +", "#0b9e73"),
    ("terminal-", "terminal −", "#7fc7a9"),
    ("internal+", "internal +", "#e8a013"),
    ("internal-", "internal −", "#2e7fe0"),
    ("deletion", "internal deletion", "#8e44ad"),
]
LINE_COLOURS = ["#111111", "#d7263d", "#0a72f5", "#7a3e9d"]

WIDTH, HEIGHT = 520, 150
LEFT, RIGHT, TOP, BOTTOM = 52, 8, 10, 20


def _x(i: float, n: int) -> float:
    return LEFT + i * (WIDTH - LEFT - RIGHT) / max(1, n - 1)


def _y(value: float, ymax: float) -> float:
    return TOP + HEIGHT * (1 - value / ymax)


def _sample(values: list, columns: int) -> list:
    """`columns` values of a profile, evenly spaced (all, if fewer)."""
    n = len(values)
    if n <= columns:
        return list(values)
    return [values[int(i * (n - 1) / (columns - 1))] for i in range(columns)]


def _frame(ymax: float, length: int, body: list[str]) -> str:
    bottom = TOP + HEIGHT
    axes = (
        f'<line x1="{LEFT}" y1="{TOP}" x2="{LEFT}" y2="{bottom}" '
        f'stroke="currentColor" opacity="0.4"/>'
        f'<line x1="{LEFT}" y1="{bottom}" x2="{WIDTH - RIGHT}" '
        f'y2="{bottom}" stroke="currentColor" opacity="0.4"/>'
        f'<text x="{LEFT - 4}" y="{TOP + 4}" text-anchor="end" '
        f'font-size="10" fill="currentColor">{ymax:,.0f}</text>'
        f'<text x="{LEFT - 4}" y="{bottom}" text-anchor="end" '
        f'font-size="10" fill="currentColor">0</text>'
        f'<text x="{LEFT}" y="{bottom + 14}" font-size="10" '
        f'fill="currentColor">1</text>'
        f'<text x="{WIDTH - RIGHT}" y="{bottom + 14}" text-anchor="end" '
        f'font-size="10" fill="currentColor">{length:,}</text>')
    return (f'<svg viewBox="0 0 {WIDTH} {TOP + HEIGHT + BOTTOM}" '
            f'width="100%" style="max-width:{WIDTH}px;height:auto" '
            f'role="img">{axes}{"".join(body)}</svg>')


def stacked_coverage(profiles: dict[str, list[int]], bands: list[tuple],
                     ymax: float, lines: dict[str, list[int]] | None = None,
                     marker: int | None = None) -> str:
    """A stacked-area plot of depth along a reference by category, with
    optional depth lines (values below 0: no line there) and a vertical
    marker at a 1-based position."""
    length = len(next(iter(profiles.values())))
    columns = WIDTH - LEFT - RIGHT
    base = None
    body = []
    for key, _, colour in bands:
        values = _sample(profiles.get(key, [0] * length), columns)
        n = len(values)
        if base is None:
            base = [0] * n
        top = [b + v for b, v in zip(base, values)]
        if any(values):
            points = [f"{_x(i, n):.1f},{_y(top[i], ymax):.1f}"
                      for i in range(n)]
            points += [f"{_x(i, n):.1f},{_y(base[i], ymax):.1f}"
                       for i in range(n - 1, -1, -1)]
            body.append(f'<polygon points="{" ".join(points)}" '
                        f'fill="{colour}"/>')
        base = top
    for k, (label, values) in enumerate((lines or {}).items()):
        sampled = _sample(values, columns)
        n = len(sampled)
        path, pen = [], False
        for i, v in enumerate(sampled):
            if v < 0:
                pen = False
                continue
            path.append(f"{'L' if pen else 'M'}{_x(i, n):.1f},"
                        f"{_y(min(v, ymax), ymax):.1f}")
            pen = True
        dash = ' stroke-dasharray="5 3"' if k % 2 else ""
        body.append(f'<path d="{" ".join(path)}" fill="none" '
                    f'stroke="{LINE_COLOURS[k % len(LINE_COLOURS)]}" '
                    f'stroke-width="1.3"{dash}><title>{escape(label)}'
                    f'</title></path>')
    if marker is not None and length > 1:
        x = _x(marker - 1, length)
        body.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" '
                    f'y2="{TOP + HEIGHT}" stroke="#d7263d" '
                    f'stroke-width="1" stroke-dasharray="2 2"/>')
    return _frame(ymax, length, body)


def legend(bands: list[tuple], lines: list[str] = ()) -> str:
    items = [f'<span><i style="background:{colour}"></i>{escape(label)}'
             f'</span>' for _, label, colour in reversed(bands)]
    for k, label in enumerate(lines):
        colour = LINE_COLOURS[k % len(LINE_COLOURS)]
        style = "dashed" if k % 2 else "solid"
        items.append(f'<span><i class="line" style="border-top:2px {style} '
                     f'{colour}"></i>{escape(label)}</span>')
    return f'<div class="legend">{"".join(items)}</div>'


def length_histogram(measured: list[int], said: list[int] | None,
                     bin_size: int, reference_length: int) -> str:
    """Bars: the lengths measured; line: those of flucut."""
    n = max(len(measured), len(said or []),
            reference_length // bin_size + 5)
    ymax = max(max(measured, default=0), max(said or [0]), 1)
    width = (WIDTH - LEFT - RIGHT) / n
    body = []
    for i, v in enumerate(measured):
        if v:
            y = _y(v, ymax)
            body.append(f'<rect x="{LEFT + i * width:.1f}" y="{y:.1f}" '
                        f'width="{max(width, 0.6):.1f}" '
                        f'height="{TOP + HEIGHT - y:.1f}" fill="#0a72f5"/>')
    if said:
        points = " ".join(
            f"{LEFT + (i + 0.5) * width:.1f},{_y(v, ymax):.1f}"
            for i, v in enumerate(said))
        body.append(f'<polyline points="{points}" fill="none" '
                    f'stroke="#f84a2a" stroke-width="1.2"/>')
    x = LEFT + reference_length / bin_size * width
    body.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" '
                f'y2="{TOP + HEIGHT}" stroke="currentColor" opacity="0.5" '
                f'stroke-dasharray="2 2"/>')
    return _frame(ymax, n * bin_size, body)


def bars(items: list[tuple[str, float, str]], total: float) -> str:
    """Horizontal bars: (label, value, colour) as shares of `total`."""
    rows = []
    for label, value, colour in items:
        share = 100 * value / total if total else 0
        rows.append(
            f'<div class="bar"><span>{escape(label)}</span>'
            f'<div><i style="width:{share:.2f}%;background:{colour}"></i>'
            f'</div><b>{value:,.0f}</b>'
            f'<em>{escape(pct(value, total) if total else "0.0%")}</em>'
            f'</div>')
    return f'<div class="bars">{"".join(rows)}</div>'
