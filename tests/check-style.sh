#!/bin/bash
#
# Check the formatting rules of the developer guide (libcbaia
# docs/dev_manual.md, "Source files") on the C sources of flucut:
# lines of at most 80 characters, no tabs, no trailing whitespace, a
# copyright line at the start of each file, include guards named after
# the header file. The same rules hold for the HTML templates
# (src/*.html), with a copyright comment, except that only their comment
# lines must be at most 80 characters long: the other lines are written
# to the report as they are (wrapping a paragraph would change it). Exits
# with status 1 and lists the problems, if any.
#
# The Python and shell scripts (*.py and *.sh outside ext/) are
# checked too: no tabs or trailing whitespace, lines of at
# most 80 characters (except lines with a URL and those listed in
# LONG_LINES_OK) and, for Python, 2 spaces per indentation level (each
# indented block 2 spaces deeper than the enclosing one; continuation
# lines, comments and strings are not checked).
#
# Usage: tests/check-style.sh   (from the repository root)
#

set -uo pipefail
export LC_ALL=C
FILES=$(ls src/*.[ch] tests/fuzz/*.[ch] 2>/dev/null)
TEMPLATES=$(ls src/*.html 2>/dev/null)
# the scripts tracked by git (without git: those found outside ext/)
SCRIPTS=$(git ls-files -- '*.py' '*.sh' ':!ext' 2>/dev/null ||
  find . -path ./ext -prune -o -path './.*' -prune -o \
    \( -name '*.py' -o -name '*.sh' \) -print | sed 's|^\./||' | sort)
PYTHON=$(for f in $SCRIPTS; do case $f in *.sh) ;; *) echo $f ;; esac; done)
# lines of the scripts allowed to be longer than 80 characters, one per
# line, as <file>:<extended regular expression> (none at the moment; lines
# with a URL are always allowed)
LONG_LINES_OK=""
status=0

problem() { echo "$1: $2"; status=1; }

# tabs and trailing whitespace in the file <f>
check_whitespace() {
  grep -n $'\t' $1 | sed "s|^|$1:|; s|$| (tab)|" | grep . && status=1
  grep -nE ' +$' $1 | sed "s|^|$1:|; s|$| (trailing whitespace)|" \
    | grep . && status=1
}

for f in $FILES; do
  # (the sources are ASCII: bytes and characters are the same)
  awk -v f=$f 'length($0) > 80 { printf "%s:%d: line longer than 80 '\
'characters (%d)\n", f, FNR, length($0) }' $f | grep . && status=1
  check_whitespace $f
  head -1 $f | grep -qE '^/\* \(c\) [0-9]{4}(-[0-9]{4})? [A-Z][a-z]+ ' \
    || problem $f "missing copyright line, e.g. /* (c) 2026 Name */"
  case $f in
    *.h)
      guard=$(basename $f .h | tr 'a-z-' 'A-Z_')_H
      grep -q "^#ifndef $guard$" $f && grep -q "^#define $guard$" $f \
        || problem $f "missing include guard $guard" ;;
  esac
done
for f in $TEMPLATES; do
  awk -v f=$f '/^[[:space:]]*<!--.*-->[[:space:]]*$/ && length($0) > 80 {
    printf "%s:%d: comment line longer than 80 characters (%d)\n", f, FNR,
    length($0) }' $f | grep . && status=1
  check_whitespace $f
  head -1 $f | grep -qE '^<!-- \(c\) [0-9]{4}(-[0-9]{4})? [A-Z][a-z]+ ' \
    || problem $f "missing copyright line, e.g. <!-- (c) 2026 Name -->"
done
for f in $SCRIPTS; do
  ok=$(echo "$LONG_LINES_OK" | sed -n "s|^$f:||p")
  awk -v f=$f -v ok="$ok" 'length($0) > 80 && !/https?:\/\// &&
    (ok == "" || $0 !~ ok) { printf "%s:%d: line longer than 80 '\
'characters (%d)\n", f, FNR, length($0) }' $f | grep . && status=1
  check_whitespace $f
done
# indentation of the Python code, from its tokens: each indented block is
# 2 spaces deeper than the enclosing one
python3 - $PYTHON <<'EOF' || status=1
import sys, tokenize
bad = False
for f in sys.argv[1:]:
  with open(f, "rb") as fh:
    levels = [0]
    for t in tokenize.tokenize(fh.readline):
      if t.type == tokenize.INDENT:
        if len(t.string) != levels[-1] + 2:
          print(f"{f}:{t.start[0]}: indentation of {len(t.string)} "
                f"spaces, expected {levels[-1] + 2}")
          bad = True
        levels.append(len(t.string))
      elif t.type == tokenize.DEDENT:
        levels.pop()
sys.exit(1 if bad else 0)
EOF
[ $status -eq 0 ] && echo "[OK] coding style ($(echo $FILES $TEMPLATES \
  $SCRIPTS | wc -w) files)"
exit $status
