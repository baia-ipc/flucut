# Turn metadata/markers.yaml into the table of the program.
#
# Writes a header with the names of the markers, their number and, for a
# translation unit which defines MARKERS_TABLE_DEFINE, the table itself
# (see markers.h). Run by the Makefile; it needs nothing but awk, so
# that building flucut still needs only a C compiler and make.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/markers-table.awk \
#          metadata/markers.yaml > out.h

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
# a site: each of its fields becomes a constant named after the site
# and the field, a number as it is written and anything else as a string
function site_macros(   key, out) {
  out = ""
  for (key in f)
    out = out sprintf("#define %s_%s %s\n", toupper(name), toupper(key),
                      f[key] ~ /^-?[0-9]+$/ ? f[key] : cstr(f[key]))
  return out
}

function flush(   i) {
  if (kind == "")
    return
  if (kind == "site") {
    sites[n_sites] = site_macros()
    n_sites++
    kind = ""
    name = ""
    split("", f)
    return
  }
  if (kind != "marker")
    fail("unknown record \"" kind "\"")
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
  kind = ""
  name = ""
  split("", f)
}

BEGIN { n = 0; n_sites = 0 }

END {
  if (n == 0) fail("no marker in the file")
  print "/* Generated from metadata/markers.yaml by scripts/markers-table.awk."
  print "   Do not edit: edit the metadata file. */"
  print ""
  # the names are included once; the table is included by the one
  # translation unit which asks for it, which may be after them
  print "#ifndef MARKERS_TABLE_H"
  print "#define MARKERS_TABLE_H"
  print ""
  for (i = 0; i < n_sites; i++) {
    print "/* a site read beside a conserved peptide, reported in its own"
    print "   way (metadata/markers.yaml) */"
    printf "%s\n", sites[i]
  }
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
