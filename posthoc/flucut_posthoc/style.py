"""The stylesheet of the HTML pages: the visual language of the report
of flucut (colours, cards, tiles, badges, strips), light and dark.

The origins of the reads (full, int, del, oth) have their own hues,
never those of the states (ok, warn, bad), which always come with a
symbol or a word."""

STYLE = """
:root { --bg:#ffffff; --card:#f6f5f2; --ink:#0b0b0b; --ink2:#52514e;
  --grid:#e4e3df; --primary:#0a72f5; --accent:#f84a2a; --brand:#03215a;
  --full:#2a78d6; --int:#d55181; --del:#6250d6; --oth:#b9b8b2;
  --short:#8a8983; --p1:#86b6ef;
  --ok:#0ca30c; --warn:#c98500; --bad:#d03b3b;
  --okbg:#e3f4e3; --warnbg:#fdf0d0; --badbg:#f9e1e1; --nabg:#ecebe7; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#1b1b1a; --card:#262624; --ink:#f2f1ee; --ink2:#b5b3ad;
    --grid:#3a3936; --primary:#4d9bff; --accent:#ff6b4f; --brand:#9cc2ff;
    --full:#3987e5; --int:#d55181; --del:#9085e9; --oth:#6d6b66;
    --short:#a3a199; --p1:#1c5cab;
    --ok:#0ca30c; --warn:#fab219; --bad:#e66767;
    --okbg:#17361a; --warnbg:#3d3115; --badbg:#432020; --nabg:#33322f; } }
body { background:var(--bg); color:var(--ink); margin:0; line-height:1.45;
  font-family:system-ui,-apple-system,'Segoe UI',Roboto,sans-serif; }
main { max-width:1280px; margin:0 auto; padding:16px; }
h1 { font-size:1.6em; margin:0.4em 0 0.1em; color:var(--brand); }
h1 small { font-size:0.5em; font-weight:700; letter-spacing:0.35em;
  text-transform:uppercase; color:var(--ink2); display:block; }
h2 { font-size:1.3em; margin:1.6em 0 0.4em; padding-top:0.5em;
  border-top:2px solid var(--accent); color:var(--brand); }
h3 { margin:0 0 8px; font-size:1em; color:var(--brand); }
details.sec { margin-top:1.6em; padding-top:0.5em;
  border-top:2px solid var(--accent); }
details.sec > summary h2 { display:inline; margin:0; padding:0;
  border:0; }
details.card > summary h3 { display:inline; margin:0; }
details.sec > summary, details.card > summary { list-style:none;
  font-size:1em; color:inherit; user-select:none; margin:0 0 8px; }
details.sec > summary::-webkit-details-marker,
details.card > summary::-webkit-details-marker { display:none; }
details.sec > summary::before, details.card > summary::before {
  content:"\\25b8"; color:var(--accent); display:inline-block;
  width:1.1em; }
details.sec[open] > summary::before,
details.card[open] > summary::before { content:"\\25be"; }
details.card:not([open]) > summary { margin:0; }
details.sec > summary:hover h2,
details.card > summary:hover h3 { text-decoration:underline; }
h4 { margin:12px 0 2px; font-size:0.9em; }
a { color:var(--primary); }
.sub { color:var(--ink2); margin:0; font-size:0.9em; }
.note { color:var(--ink2); font-size:0.85em; margin:4px 0; }
code { font-family:ui-monospace,monospace; font-size:0.9em;
  word-break:break-all; }
table { border-collapse:collapse; font-size:0.9em; margin:0.4em 0; }
th, td { padding:3px 10px; text-align:right;
  border-bottom:1px solid var(--grid);
  font-variant-numeric:tabular-nums; }
th { color:var(--ink2); font-weight:600; vertical-align:bottom; }
td.l, th.l { text-align:left; }
th.g { text-align:center; border-bottom:0; padding-bottom:0; }
td small { color:var(--ink2); }
td.seq { font-family:ui-monospace,monospace; text-align:left; }
td.dim, .dim { color:var(--ink2); }
td.hi { font-weight:700; }
td.sc { background:var(--nabg); border-left:2px solid var(--bg); }
td.sc.ok { background:var(--okbg); }
td.sc.warn { background:var(--warnbg); }
td.sc.bad { background:var(--badbg); }
.scroll { overflow-x:auto; }
details > summary { cursor:pointer; color:var(--ink2); font-size:0.85em;
  margin:4px 0; }
details.fold { border-top:1px solid var(--grid); padding:4px 0; }
details.fold > summary { font-size:1em; font-weight:600;
  color:var(--brand); padding:4px 0; }
details.fold > summary small { font-weight:400; color:var(--ink2);
  margin-left:8px; }
.params td { text-align:left; overflow-wrap:anywhere; }
.terms { margin-top:28px; }
.terms td { text-align:left; vertical-align:top; }
.terms td:first-child { white-space:nowrap; font-weight:600; }
.chips { display:flex; flex-wrap:wrap; align-items:center; gap:6px 8px;
  margin:6px 0 10px; font-size:0.9em; }
.st { display:inline-block; padding:1px 9px; border-radius:10px;
  font-weight:600; white-space:nowrap; background:var(--nabg); }
.st.ok { background:var(--okbg); box-shadow:inset 3px 0 var(--ok); }
.st.warn { background:var(--warnbg); box-shadow:inset 3px 0 var(--warn); }
.st.bad { background:var(--badbg); box-shadow:inset 3px 0 var(--bad); }
td .st, .sname .st { font-size:0.85em; }
.chip { padding:1px 9px; border-radius:10px; border:1px solid var(--grid);
  color:var(--ink2); }
.tiles { display:flex; flex-wrap:wrap; gap:10px; margin:0.6em 0; }
.tile { background:var(--card); border-radius:8px; padding:8px 14px;
  flex:1 1 150px; min-width:150px; box-sizing:border-box; }
.tile .v { font-size:1.5em; font-weight:600; }
.tile .v .st { font-size:0.55em; vertical-align:4px; margin-left:4px; }
.tile .l { color:var(--ink2); font-size:0.85em; }
.tile .u { color:var(--ink2); font-size:0.75em; }
.cols { display:grid; grid-template-columns:3fr 2fr; gap:12px;
  align-items:start; }
.card { background:var(--card); border-radius:10px; padding:10px 14px;
  margin:12px 0 0; min-width:0; }
.card table { width:100%; }
.bar { margin:0 0 12px; }
.bl { font-size:0.9em; font-weight:600; margin-bottom:3px; }
.bl .note { font-weight:400; }
.sbar { display:flex; gap:2px; height:16px; }
.sbar i { display:block; min-width:2px; border-radius:2px; }
.lg { font-size:0.8em; margin-top:3px; line-height:1.7; }
.lg span { white-space:nowrap; margin-right:14px; }
.lg i { width:10px; height:10px; border-radius:2px; display:inline-block;
  vertical-align:-1px; margin-right:4px; }
i.c-full { background:var(--full); }
i.c-int { background:var(--int); }
i.c-del { background:var(--del); }
i.c-oth { background:var(--oth); }
i.c-p1 { background:var(--p1); }
i.c-short { background:var(--short); }
i.c-low { background:var(--bad); }
i.c-nw { background:repeating-linear-gradient(45deg,var(--oth) 0 2px,
  transparent 2px 5px); outline:1px solid var(--oth); outline-offset:-1px; }
i.c-line, i.c-dash { height:0; width:18px; border-radius:0;
  vertical-align:3px; border-top:2px solid var(--ink); }
i.c-dash { border-top-style:dashed; }
i.c-site { width:0; height:11px; border-left:2px solid var(--int);
  border-radius:0; }
i.t-ok { background:var(--okbg); outline:1px solid var(--ok); }
i.t-warn { background:var(--warnbg); outline:1px solid var(--warn); }
i.t-bad { background:var(--badbg); outline:1px solid var(--bad); }
i.t-na { background:var(--nabg); outline:1px solid var(--oth); }
.kv { display:grid; grid-template-columns:repeat(auto-fill,
  minmax(190px,1fr)); gap:0 18px; margin:4px 0 0; font-size:0.85em; }
.kv > div { display:flex; justify-content:space-between; gap:8px;
  border-bottom:1px solid var(--grid); padding:2px 0; }
.kv dt { color:var(--ink2); }
.kv dd { margin:0; font-variant-numeric:tabular-nums; white-space:nowrap; }
.kv small { color:var(--ink2); }
.seghead, .seg { display:grid; gap:4px 18px; align-items:center;
  grid-template-columns:320px minmax(0,1.9fr) minmax(0,1fr); }
.seghead { font-size:0.9em; align-items:end; }
.seg { border-top:1px solid var(--grid); padding:8px 0; }
.sname { display:flex; align-items:center; gap:8px; flex-wrap:wrap; }
.sname b { font-size:1.15em; }
.sinfo .note, .slen .note, .scov .note { margin:0; font-size:0.75em; }
.sinfo table { width:100%; font-size:0.8em; }
.sinfo th, .sinfo td { padding:2px 0 2px 8px; white-space:nowrap; }
.sinfo td.l, .sinfo th.l { padding-left:0; color:var(--ink2); }
.sinfo td.kh { text-align:left; padding:6px 0 0; border:0;
  color:var(--ink2); font-size:0.9em; }
.sinfo td .st { font-size:0.9em; padding:0 7px; }
.seg svg, .card svg { width:100%; }
svg { max-width:100%; height:auto; display:block; }
svg text { fill:var(--ink); font-size:11px; font-family:inherit; }
svg .muted { fill:var(--ink2); }
svg .axis { stroke:var(--ink2); }
svg .c-full { fill:var(--full); }
svg .c-int { fill:var(--int); }
svg .c-del { fill:var(--del); }
svg .c-oth { fill:var(--oth); }
svg .c-short { fill:var(--short); }
svg .c-low { fill:var(--bad); }
svg .mark { stroke:var(--ink2); stroke-dasharray:4 3; }
svg .site { stroke:var(--int); stroke-width:1; }
svg .sitel { font-weight:600; }
svg .depth { fill:none; stroke:var(--ink); stroke-width:1.3;
  stroke-linejoin:round; }
svg .depth.d1 { stroke-dasharray:5 3; }
svg .depth.d2 { stroke-dasharray:1 3; }
svg .said { fill:none; stroke:var(--ink); stroke-width:1.2;
  stroke-linejoin:round; }
.samples { width:100%; }
.samples td { white-space:nowrap; }
.samples td.sc { text-align:center; }
@media (max-width:900px) {
  .cols, .seghead, .seg { grid-template-columns:minmax(0,1fr); }
  .seghead { display:none; } }
@media print { .card, .seg { break-inside:avoid; }
  details > summary::before { display:none; } }
"""
