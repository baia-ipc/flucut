# Turn metadata/examples.yaml into the constants of the program.
#
# Writes a header with one constant per example sequence, used by the
# unit tests (see the metadata file). Run by the Makefile.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/examples-table.awk \
#          metadata/examples.yaml > out.h

function flush(   macro) {
  if (kind == "example") {
    if (f["bases"] !~ /^[ACGTURYSWKMBDHVN]+$/)
      fail("the example " name " has no bases")
    macro = toupper(name)
    gsub(/[^A-Z0-9]/, "_", macro)
    rows[n] = sprintf("#define EXAMPLE_%s %s", macro, cstr(f["bases"]))
    n++
  } else if (kind != "") {
    fail("unknown record \"" kind "\"")
  }
  kind = ""
  name = ""
  split("", f)
}

BEGIN { n = 0 }

END {
  if (n == 0)
    fail("no example in the file")
  print "/* Generated from metadata/examples.yaml by"
  print "   scripts/examples-table.awk. Do not edit: edit the metadata. */"
  print ""
  print "#ifndef EXAMPLES_TABLE_H"
  print "#define EXAMPLES_TABLE_H"
  print ""
  for (i = 0; i < n; i++)
    print rows[i]
  print ""
  print "#endif"
}
