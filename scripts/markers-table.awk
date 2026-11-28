# Turn metadata/markers.yaml into the table of the program.
#
# Writes a header with the names of the markers, their number and, for a
# translation unit which defines MARKERS_TABLE_DEFINE, the table itself
# (see markers.h). Run by the Makefile; it needs nothing but awk, so
# that building flucut still needs only a C compiler and make.
#
# Usage: awk -f scripts/markers-table.awk metadata/markers.yaml > out.h

function fail(msg) {
  print "markers-table.awk: " FILENAME " line " FNR ": " msg > "/dev/stderr"
  exit 1
}

# a C string literal of <s>, which may be empty
function cstr(s,   out) {
  out = s
  gsub(/\\/, "\\\\", out)
  gsub(/"/, "\\\"", out)
  return "\"" out "\""
}

# the identifier of the marker <s>: upper case, - and . becoming _
function ident(s,   out) {
  out = toupper(s)
  gsub(/[^A-Z0-9]/, "_", out)
  return "MARKER_" out
}

function yesno(s) {
  if (s == "yes") return "true"
  if (s == "no" || s == "") return "false"
  fail("expected yes or no, found \"" s "\"")
}

function kind_of(s) {
  if (s == "residues") return "MARKER_RESIDUES"
  if (s == "to_stop") return "MARKER_TO_STOP"
  if (s == "frame_length") return "MARKER_FRAME_LENGTH"
  if (s == "span") return "MARKER_SPAN"
  fail("unknown kind \"" s "\"")
}

# the record read so far, written out
function flush(   i) {
  if (name == "") return
  names[n] = name
  rows[n] = \
    "  [" ident(name) "] = {\n" \
    "    .name = " cstr(name) ", .protein = " cstr(f["protein"]) \
      ", .segment = " (f["segment"] + 0) ",\n" \
    "    .peptide = " cstr(f["anchor"]) ",\n" \
    "    .peptide2 = " (f["anchor2"] == "" ? "NULL" : cstr(f["anchor2"])) \
      ",\n" \
    "    .offset = " (f["offset"] + 0) ", .window = " (f["window"] + 0) \
      ", .kind = " kind_of(f["kind"]) ",\n" \
    "    .numbering = " cstr(f["numbering"]) ",\n" \
    "    .positions = " cstr(f["positions"]) ",\n" \
    "    .usual = " cstr(f["usual"]) ", .attention = " cstr(f["attention"]) \
      ",\n" \
    "    .attention_below = " (f["below"] + 0) ", .flag = " yesno(f["flag"]) \
      ",\n" \
    "    .meaning = " cstr(f["meaning"]) ",\n" \
    "  },"
  n++
  name = ""
  split("", f)
}

BEGIN { n = 0; name = "" }

/^[ \t]*#/ { next }
/^[ \t]*$/ { next }

{
  line = $0
  sub(/^- /, "  ", line)          # a record begins with "- "
  if ($0 ~ /^- /) flush()
  if (line !~ /^  [a-z0-9]+: /) fail("expected \"  key: value\"")
  key = line
  sub(/^  /, "", key)
  sub(/:.*/, "", key)
  value = line
  sub(/^  [a-z0-9]+: /, "", value)
  sub(/[ \t]+$/, "", value)
  # a value may be quoted, for one which YAML would read otherwise
  if (value ~ /^".*"$/) {
    value = substr(value, 2, length(value) - 2)
  }
  if (key == "name") {
    name = value
  } else {
    if (name == "") fail("a field before the first name")
    f[key] = value
  }
}

END {
  flush()
  if (n == 0) fail("no marker in the file")
  print "/* Generated from metadata/markers.yaml by scripts/markers-table.awk."
  print "   Do not edit: edit the metadata file. */"
  print ""
  # the names are included once; the table is included by the one
  # translation unit which asks for it, which may be after them
  print "#ifndef MARKERS_TABLE_H"
  print "#define MARKERS_TABLE_H"
  print ""
  print "/* the markers, in the order they are reported */"
  print "enum {"
  for (i = 0; i < n; i++)
    printf "  %s,\n", ident(names[i])
  print "  MARKER_N_MARKERS"
  print "};"
  print ""
  print "#endif"
  print ""
  print "#if defined(MARKERS_TABLE_DEFINE) && !defined(MARKERS_TABLE_DEFINED)"
  print "#define MARKERS_TABLE_DEFINED"
  print "static const marker_def_t MARKERS_TABLE[MARKER_N_MARKERS] = {"
  for (i = 0; i < n; i++)
    print rows[i]
  print "};"
  print "#endif"
}
