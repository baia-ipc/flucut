#!/usr/bin/env bash
#
# List the object files to link into a program: for each header included
# (directly or indirectly, also through the .c files of the modules) by the
# program source, the object of the corresponding .c file, if it exists.
# The sources are in several directories, searched in the given order (the
# headers are included without the directory, and are unique over all of
# them); the objects are all in the same directory.
#
# Usage: scripts/objdeps.sh "<srcdir> ..." <objdir> <program.c>
#
# (compatible with bash 3, e.g. on macOS)
#

SRCDIRS=$1
OBJDIR=$2
PROG=$3

# print the path of the file <name> in the first source directory with it
find_source() {
  for dir in $SRCDIRS; do
    [ -e "$dir/$1" ] && { echo "$dir/$1"; return 0; }
  done
  return 1
}

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
    path=$(find_source "$header") && queue+=("$path")
    module=${header%.h}
    if [ "$module" != "$header" ] && path=$(find_source "$module.c"); then
      objs+=("$OBJDIR/$module.o")
      queue+=("$path")
    fi
  done
done

echo "${objs[@]+"${objs[@]}"}"
