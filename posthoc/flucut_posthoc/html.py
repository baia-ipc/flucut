"""Small helpers for writing the HTML report."""

import os
from html import escape

from .style import STYLE

SYMBOLS = {"ok": "✓", "warn": "!", "bad": "✕", "na": "·"}


def num(value) -> str:
    if value is None:
        return "–"
    if isinstance(value, float):
        return f"{value:,.1f}"
    return f"{value:,}"


def pct(part: float, total: float, digits: int = 1) -> str:
    """`part` as a percentage of `total`, with `digits` decimals. Zero and
    100 are written only for none and for a part equal to the total: with
    one decimal, a part which is not zero, but below 0.1%, gives "<0.1%",
    and one above 99.9%, but below the total, ">99.9%" (a part just above
    the total, e.g. of a ratio of two measures, gives ">100.0%")."""
    if not total:
        return "–"
    value = 100 * part / total
    step = 10 ** -digits
    text = f"{value:.{digits}f}%"
    if part != 0 and 0 < value < step:
        return f"<{step:.{digits}f}%"
    if part < total and value > 100 - step:
        return f">{100 - step:.{digits}f}%"
    if part > total and text == f"{100:.{digits}f}%":
        return ">" + text
    return text


class Markup(str):
    """Text which is already HTML."""


def table(headers: list[str], rows: list[list], left: int = 1,
          cls: str = "", groups: list[tuple[str, int]] | None = None
          ) -> str:
    """A table; a cell is a value or a (value, class) pair; the first
    `left` columns are left-aligned. `groups`: a line of titles above the
    headers, each over a number of columns."""
    out = [f'<div class="scroll"><table{f" class={cls}" if cls else ""}>']
    if groups:
        out.append("<tr>" + "".join(
            f'<th colspan="{span}" class="g">{escape(title)}</th>'
            for title, span in groups) + "</tr>")
    out.append("<tr>")
    for i, h in enumerate(headers):
        out.append(f'<th{" class=l" if i < left else ""}>{escape(h)}</th>')
    out.append("</tr>")
    for row in rows:
        out.append("<tr>")
        for i, cell in enumerate(row):
            value, kind = cell if isinstance(cell, tuple) else (cell, "")
            if i < left:
                kind = (kind + " l").strip()
            text = value if isinstance(value, Markup) else \
                num(value) if isinstance(value, (int, float)) or \
                value is None else escape(str(value))
            out.append(f'<td{f" class=\"{kind}\"" if kind else ""}>{text}</td>')
        out.append("</tr>")
    out.append("</table></div>")
    return "".join(out)


def cross_table(corner: str, counts: dict, rows: list[str], cols: list[str],
                agree: dict | None = None, rest: bool = True) -> str:
    """A table of counts[row][column] with row totals; `agree` maps a row
    to the column which agrees with it (shown in bold); `rest` adds a
    column with the counts of the columns not named."""
    body = []
    for r in rows:
        line = counts.get(r, {})
        total = sum(line.values())
        if total == 0:
            continue
        cells = [r]
        for c in cols:
            v = line.get(c, 0)
            cls = ""
            if agree and r in agree and c in agree.values() and v:
                cls = "hi" if agree[r] == c else ""
            cells.append((v, cls if v else "dim"))
        other = total - sum(line.get(c, 0) for c in cols)
        body.append(cells + ([(other, "dim" if not other else "")]
                             if rest else []) + [total])
    return table([corner] + cols + (["other"] if rest else []) + ["total"],
                 body)


def badge(status: str, text: str) -> str:
    """A state: its symbol and `text`, on the tint of the state
    (`status`: ok, warn, bad or na)."""
    return (f'<span class="st {status}">{SYMBOLS[status]} '
            f"{escape(text)}</span>")


def chips(*items: str) -> str:
    return f'<p class="chips">{"".join(i for i in items if i)}</p>'


def chip(text: str) -> str:
    return f'<span class="chip">{escape(text)}</span>'


def tile(value: str, label: str, sub: str = "", status: str = "") -> str:
    """A large value with its name and, in a smaller line, what it is a
    share of; `status` adds the symbol of a state to the value."""
    mark = f'<span class="st {status}">{SYMBOLS[status]}</span>' \
        if status else ""
    return (f'<div class="tile"><div class="v">{escape(value)}{mark}</div>'
            f'<div class="l">{escape(label)}</div>'
            + (f'<div class="u">{escape(sub)}</div>' if sub else "")
            + "</div>")


def tiles(*items: str) -> str:
    return f'<div class="tiles">{"".join(items)}</div>'


def card(title: str, *parts: str, open_: bool = True) -> str:
    """A subsection: a box which folds under its title."""
    return (f'<details class="card"{" open" if open_ else ""}>'
            f"<summary><h3>{escape(title)}</h3></summary>"
            f'{"".join(parts)}</details>')


def legend(items: list[tuple[str, str]]) -> str:
    """Swatches (class, label) with their names."""
    return '<div class="lg">' + "".join(
        f'<span><i class="{cls}"></i>{escape(label)}</span>'
        for cls, label in items) + "</div>"


def stacked_bar(title: str, of: str, items: list[tuple[str, str, float]]
                ) -> str:
    """A bar of the parts (class, label, value) of a whole, with the
    share of each in the legend; `of` names the whole in the title."""
    total = sum(v for _, _, v in items)
    if not total:
        return ""
    parts, names = [], []
    for cls, label, value in items:
        if value:
            parts.append(f'<i class="{cls}" style="flex:{value}" '
                         f'title="{escape(label)}: {value:,}"></i>')
        names.append(f'<span><i class="{cls}"></i>{escape(label)} '
                     f"<b>{escape(pct(value, total))}</b></span>")
    return (f'<div class="bar"><div class="bl">{escape(title)} '
            f'<span class="note">of {total:,} {escape(of)}</span></div>'
            f'<div class="sbar">{"".join(parts)}</div>'
            f'<div class="lg">{"".join(names)}</div></div>')


def key_values(items: list[tuple[str, str, str]]) -> str:
    """Small counts with their names: (name, value, smaller text)."""
    return '<dl class="kv">' + "".join(
        f"<div><dt>{escape(k)}</dt><dd>{escape(v)}"
        + (f" <small>{escape(small)}</small>" if small else "")
        + "</dd></div>" for k, v, small in items) + "</dl>"


def fold(title: str, *parts: str, hint: str = "", anchor: str = "") -> str:
    """Content shown on request, under a title."""
    body = "".join(parts)
    if not body:
        return ""
    return (f'<details class="fold"{f" id={anchor}" if anchor else ""}>'
            f"<summary>{escape(title)}"
            + (f"<small>{escape(hint)}</small>" if hint else "")
            + f"</summary>{body}</details>")


def more(title: str, *parts: str) -> str:
    """A few more lines shown on request, inside a card."""
    return f'<details><summary>{escape(title)}</summary>' \
           f'{"".join(parts)}</details>'


def terms(rows: list[tuple[str, str]]) -> str:
    """The terms of a page, explained in a table shown on request."""
    body = "".join(f"<tr><td>{escape(k)}</td><td>{escape(v)}</td></tr>"
                   for k, v in rows)
    return ('<details class="terms"><summary>Terms used in this page'
            f"</summary><table>{body}</table></details>")


def note(text: str) -> str:
    """One line of explanation."""
    return f'<p class="note">{escape(text)}</p>'


def section(anchor: str, title: str, *parts: str, open_: bool = True
            ) -> str:
    """A section of a page, which folds under its title."""
    return (f'<details class="sec" id="{anchor}"{" open" if open_ else ""}>'
            f"<summary><h2>{escape(title)}</h2></summary>"
            f'{"".join(parts)}</details>')


# a link inside the page opens the folded sections around its target
SCRIPT = """<script>
function openTarget() {
  var t = location.hash && document.getElementById(location.hash.slice(1));
  for (var e = t; e; e = e.parentElement)
    if (e.tagName === "DETAILS") e.open = true;
  if (t) t.scrollIntoView();
}
window.addEventListener("hashchange", openTarget);
openTarget();
</script>"""


def page(title: str, *parts: str, wide: bool = False) -> str:
    """A complete page with the stylesheet; `wide`: for a table with
    many columns."""
    style = STYLE + ("main { max-width:1800px; }" if wide else "")
    return ('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,'
            f'initial-scale=1"><title>{escape(title)}</title>'
            f"<style>{style}</style></head><body><main>"
            f'{"".join(parts)}<p class="note">Generated by flucut-posthoc.'
            f"</p></main>{SCRIPT}</body></html>")


def write_page(path: str, text: str) -> None:
    with open(path + ".tmp", "w") as handle:
        handle.write(text)
    os.replace(path + ".tmp", path)
