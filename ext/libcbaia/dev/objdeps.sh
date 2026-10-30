#!/usr/bin/env bash
#
# List the object files to link into a program: for each header included
# (directly or indirectly, also through the .c files of the modules) by the
# program source, the object of the corresponding .c file, if it exists.
#
# Usage: dev/objdeps.sh <srcdir> <objdir> <program.c>
#
# (compatible with bash 3, e.g. on macOS)
#

SRCDIR=$1
OBJDIR=$2
PROG=$3

seen=" "
queue=("$PROG")
objs=()

while [ ${#queue[@]} -gt 0 ]; do
  file=${queue[0]}
  queue=("${queue[@]:1}")
  for header in $(sed -n 's/^#include "\(.*\)"/\1/p' "$file"); do
    case "$seen" in
      *" $header "*) continue ;;
    esac
    seen="$seen$header "
    [ -e "$SRCDIR/$header" ] && queue+=("$SRCDIR/$header")
    module=${header%.h}
    if [ "$module" != "$header" ] && [ -e "$SRCDIR/$module.c" ]; then
      objs+=("$OBJDIR/$module.o")
      queue+=("$SRCDIR/$module.c")
    fi
  done
done

echo "${objs[@]+"${objs[@]}"}"
