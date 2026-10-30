#!/usr/bin/env python3
"""
Terminal screenshots of the libcbaia manual: run commands in a
pseudo-terminal (so that the logger uses colors, as on a real terminal),
and write their output, with the ANSI colors, as an HTML page with a
terminal look; screenshot.js turns it into a PNG.

Usage: terminal.py OUT.html [-w COLUMNS] -- CMD [-- CMD ...]
  each CMD is a shell command, shown with a "$ " prompt before its output;
  with the suffix " @N", only the first N lines of its output are shown
"""
import html
import os
import pty
import re
import subprocess
import sys

# xterm-like palette on a dark background
PALETTE = {
  30: "#3b3b3b", 31: "#e06c75", 32: "#98c379", 33: "#d19a66",
  34: "#61afef", 35: "#c678dd", 36: "#56b6c2", 37: "#dcdfe4",
  90: "#7f848e", 91: "#ff7b85", 92: "#a8e890", 93: "#f0d070",
  94: "#7cc4ff", 95: "#e09cf5", 96: "#6fd9e6", 97: "#ffffff",
}
FG = "#dcdfe4"
BG = "#1e2127"

def run(cmd, cols):
  """output of cmd, run in a pseudo-terminal"""
  master, slave = pty.openpty()
  env = dict(os.environ, COLUMNS=str(cols), TERM="xterm")
  p = subprocess.Popen(cmd, shell=True, stdin=slave, stdout=slave,
                       stderr=slave, env=env, close_fds=True)
  os.close(slave)
  out = b""
  while True:
    try:
      data = os.read(master, 65536)
    except OSError:
      break
    if not data:
      break
    out += data
  p.wait()
  os.close(master)
  return out.decode("utf8", "replace").replace("\r\n", "\n")

def apply_cr(text):
  """carriage returns: the text after the last one of each line"""
  return "\n".join(l.split("\r")[-1] if "\r" in l.rstrip("\r") else
                   l.rstrip("\r") for l in text.split("\n"))

def to_html(text):
  out, color = [], None
  pos = 0
  for m in re.finditer(r"\x1b\[([0-9;]*)m", text):
    out.append(span(text[pos:m.start()], color))
    codes = [int(c) for c in m.group(1).split(";") if c]
    for c in codes or [0]:
      color = PALETTE.get(c) if c != 0 else None
    pos = m.end()
  out.append(span(text[pos:], color))
  return "".join(out)

def span(s, color):
  s = html.escape(re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", s))
  return f'<span style="color:{color}">{s}</span>' if color and s else s

def main():
  args = sys.argv[1:]
  outpath = args.pop(0)
  cols = 100
  if args[0] == "-w":
    cols = int(args[1])
    args = args[2:]
  cmds = " ".join(args).split(" -- ")
  cmds[0] = cmds[0].lstrip("- ").strip()
  body = []
  for k, cmd in enumerate(cmds):
    # "CMD @N": only the first N lines of the output, then "..."
    max_lines = None
    m = re.search(r"\s@(\d+)$", cmd)
    if m:
      max_lines = int(m.group(1))
      cmd = cmd[:m.start()]
    if k > 0:
      body.append("")
    body.append(f'<span style="color:{PALETTE[92]}">$</span> '
                + html.escape(cmd))
    out = apply_cr(run(cmd, cols)).split("\n")
    visible = lambda l: re.sub(r"\x1b\[[0-9;]*m", "", l).strip()
    while out and not visible(out[-1]):   # e.g. a final color reset
      out[-2:] = [out[-2] + out[-1]] if len(out) > 1 else []
    if max_lines is not None and len(out) > max_lines:
      out = out[:max_lines] + ["\x1b[0m..."]
    body.append(to_html("\n".join(out)))
  page = f"""<!doctype html><html><head><meta charset="utf-8"><style>
body{{margin:0;background:#fff}}
h1{{font-size:1px;margin:0;color:#fff}}
pre{{display:inline-block;margin:0;padding:16px 20px;background:{BG};
color:{FG};font:14px/1.4 "DejaVu Sans Mono",monospace;border-radius:6px;
white-space:pre}}
</style></head><body><h1>t</h1><pre>{chr(10).join(body)}</pre></body></html>
"""
  open(outpath, "w").write(page)

if __name__ == "__main__":
  main()
