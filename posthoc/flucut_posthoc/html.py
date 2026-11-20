"""Small helpers for writing the HTML report."""

from html import escape

STYLE = """
:root { --primary:#0a72f5; --accent:#f84a2a; --brand:#03215a;
  --bg:#ffffff; --fg:#1c2330; --card:#f3f5f9; --line:#d9dee7;
  --good:#0b9e73; --bad:#d7263d; }
@media (prefers-color-scheme: dark) {
  :root { --primary:#4d9bff; --accent:#ff6b4f; --brand:#9cc2ff;
    --bg:#12161d; --fg:#e4e8ef; --card:#1c222c; --line:#333c4a; } }
body { margin:0; background:var(--bg); color:var(--fg);
  font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1120px; margin:0 auto; padding:20px 24px 60px; }
h1 { color:var(--brand); margin:0 0 2px; font-size:1.6em; }
h2 { color:var(--brand); border-top:2px solid var(--accent);
  padding-top:14px; margin:34px 0 8px; font-size:1.25em; }
h3 { color:var(--brand); margin:20px 0 6px; font-size:1.02em; }
p { margin:6px 0 10px; max-width:860px; }
.sub, .note { color:color-mix(in srgb,var(--fg) 70%,transparent); }
.note { font-size:0.9em; }
a { color:var(--primary); }
code { font:0.88em ui-monospace,Menlo,Consolas,monospace; }
.tiles { display:grid; gap:10px; margin:14px 0;
  grid-template-columns:repeat(auto-fill,minmax(170px,1fr)); }
.tile { background:var(--card); border-radius:8px; padding:10px 12px; }
.tile b { display:block; font-size:1.45em; color:var(--brand); }
.tile span { font-size:0.85em; }
table { border-collapse:collapse; margin:8px 0 14px; font-size:0.9em; }
th, td { padding:3px 10px; text-align:right;
  border-bottom:1px solid var(--line); white-space:nowrap; }
th { font-weight:600; background:var(--card); }
th:first-child, td:first-child, td.l, th.l { text-align:left; }
td.seq { font-family:ui-monospace,Menlo,Consolas,monospace; }
td.good { color:var(--good); font-weight:600; }
td.bad { color:var(--bad); font-weight:600; }
td.dim { opacity:0.55; }
.scroll { overflow-x:auto; }
.params td { white-space:normal; overflow-wrap:anywhere;
  text-align:left; }
.plots { display:grid; gap:6px 26px;
  grid-template-columns:repeat(auto-fill,minmax(330px,1fr)); }
.pair { display:grid; gap:4px 22px; grid-template-columns:1fr 1fr;
  margin:4px 0 14px; }
@media (max-width:760px) { .pair { grid-template-columns:1fr; } }
.plot b { font-size:0.9em; }
.plot svg { display:block; }
.legend { font-size:0.85em; margin:6px 0 10px; line-height:1.8; }
.legend span { white-space:nowrap; margin-right:16px; }
.legend i { width:11px; height:11px; border-radius:2px;
  display:inline-block; vertical-align:middle; margin-right:4px; }
.legend i.line { height:0; width:18px; border-radius:0; }
.bars { max-width:760px; margin:8px 0 14px; font-size:0.9em; }
.bar { display:grid; grid-template-columns:250px 1fr 90px 60px; gap:8px;
  align-items:center; margin:3px 0; }
.bar div { background:var(--card); border-radius:3px; height:13px; }
.bar i { display:block; height:13px; border-radius:3px; }
.bar b, .bar em { text-align:right; font-style:normal; }
nav { font-size:0.9em; margin:10px 0 0; }
nav a { margin-right:14px; white-space:nowrap; }
"""


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


def table(headers: list[str], rows: list[list], left: int = 1) -> str:
    """A table; a cell is a value or a (value, class) pair; the first
    `left` columns are left-aligned."""
    out = ['<div class="scroll"><table><tr>']
    for i, h in enumerate(headers):
        out.append(f'<th{" class=l" if i < left else ""}>{escape(h)}</th>')
    out.append("</tr>")
    for row in rows:
        out.append("<tr>")
        for i, cell in enumerate(row):
            value, cls = cell if isinstance(cell, tuple) else (cell, "")
            if i < left:
                cls = (cls + " l").strip()
            text = num(value) if isinstance(value, (int, float)) or \
                value is None else escape(str(value))
            out.append(f'<td{f" class=\"{cls}\"" if cls else ""}>{text}</td>')
        out.append("</tr>")
    out.append("</table></div>")
    return "".join(out)


def cross_table(corner: str, counts: dict, rows: list[str], cols: list[str],
                agree: dict | None = None) -> str:
    """A table of counts[row][column] with row totals; `agree` maps a row
    to the column which agrees with it (shown as good, the others of the
    columns named in it as bad)."""
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
                cls = "good" if agree[r] == c else "bad"
            cells.append((v, cls if v else "dim"))
        other = total - sum(line.get(c, 0) for c in cols)
        body.append(cells + [(other, "dim" if not other else ""), total])
    return table([corner] + cols + ["other", "total"], body)


def tile(label: str, value: str) -> str:
    return f'<div class="tile"><b>{escape(value)}</b>' \
           f'<span>{escape(label)}</span></div>'


def section(anchor: str, title: str, *parts: str) -> str:
    return f'<h2 id="{anchor}">{escape(title)}</h2>' + "".join(parts)


def para(text: str, cls: str = "") -> str:
    """A paragraph of trusted HTML text."""
    return f'<p{f" class={cls}" if cls else ""}>{text}</p>'
