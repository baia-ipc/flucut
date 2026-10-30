#!/bin/sh
#
# Write a data file (e.g. a CSS style or an HTML template) as a C string
# literal, to be included by a source; the lines made only of a comment
# (/* ... */ or <!-- ... -->, e.g. the copyright line) are removed, so
# that the data files can be documented without bloating the program.
# A CSS file (.css) is also compacted: the empty lines and the indentation
# are removed and each rule (up to its closing brace) is written on one
# line, so that it can be formatted for reading.
#
# Usage: dev/embed.sh <file> > <file.inc>
#   then, in a source: const char *const x =
#                      #include "file.inc"
#                        ;
#

case "$1" in
  *.css) compact=1 ;;
  *) compact=0 ;;
esac
sed -e '/^[[:space:]]*\/\*.*\*\/[[:space:]]*$/d' \
    -e '/^[[:space:]]*<!--.*-->[[:space:]]*$/d' "$1" |
  if [ $compact = 1 ]; then
    awk '{ sub(/^[[:space:]]+/, ""); sub(/[[:space:]]+$/, "") }
         $0 == "" { next }
         /\{$/ { sub(/ \{$/, "{") }
         /;$/ { sub(/: /, ":") }
         { line = line $0 }
         /\}$/ { print line; line = "" }
         END { if (line != "") print line }'
  else
    cat
  fi |
  sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' -e 's/^/"/' -e 's/$/\\n"/'
