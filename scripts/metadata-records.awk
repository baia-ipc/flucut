# Read a metadata file of the metadata directory.
#
# The files are records of one kind or another, each a line "- kind:
# name" followed by indented "key: value" lines, with blank lines and
# "#" comments ignored:
#
#   - primer: uni12
#     virus: A
#     end: left
#
# This script reads them and calls flush() at the end of each record,
# with <kind>, <name> and the fields in <f>; the script which writes the
# table defines flush() and an END which writes what it collected. Both
# are given to awk, this one first:
#
#   awk -f scripts/metadata-records.awk -f scripts/<x>-table.awk <file>
#
# A value may be quoted, for one which would otherwise be read as
# something else; ": " and "#" are not allowed in a value.

function fail(msg) {
  print FILENAME ": line " FNR ": " msg > "/dev/stderr"
  exit 1
}

# a C string literal of <s>, which may be empty
function cstr(s,   out) {
  out = s
  gsub(/\\/, "\\\\", out)
  gsub(/"/, "\\\"", out)
  return "\"" out "\""
}

BEGIN { kind = ""; name = "" }

/^[ \t]*#/ { next }
/^[ \t]*$/ { next }

{
  line = $0
  sub(/^- /, "  ", line)          # a record begins with "- "
  if (line !~ /^  [a-z0-9]+: /)
    fail("expected \"  key: value\"")
  key = line
  sub(/^  /, "", key)
  sub(/:.*/, "", key)
  value = line
  sub(/^  [a-z0-9]+: /, "", value)
  sub(/[ \t]+$/, "", value)
  if (value ~ /^".*"$/)
    value = substr(value, 2, length(value) - 2)
  if ($0 ~ /^- /) {
    flush()
    kind = key
    name = value
  } else {
    if (kind == "")
      fail("a field before the first record")
    f[key] = value
  }
}

END { flush() }
