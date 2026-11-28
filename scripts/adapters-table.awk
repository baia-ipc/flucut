# Turn metadata/adapters.yaml into the table of the program.
#
# Writes a header with the names of the kits, for the option --adapters,
# and, for the translation unit which defines ADAPTERS_TABLE_DEFINE, the
# table itself (see adapters.h). Run by the Makefile.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/adapters-table.awk \
#          metadata/adapters.yaml > out.h

# the bases of <list>, a space-separated list of named sequences and
# literal bases; "none" gives NULL
function bases_of(list,   parts, i, k, out, p) {
  if (list == "none" || list == "")
    return "NULL"
  out = ""
  k = split(list, parts, / +/)
  for (i = 1; i <= k; i++) {
    p = parts[i]
    if (p == "")
      continue
    if (p in bases)
      out = out bases[p]
    else if (p ~ /^[ACGTURYSWKMBDHVN]+$/)
      out = out p
    else
      fail("\"" p "\" is neither a named sequence nor bases")
  }
  if (out == "")
    fail("an empty sequence")
  return cstr(out)
}

function flush() {
  if (kind == "sequence") {
    if (f["bases"] == "")
      fail("the sequence " name " has no bases")
    bases[name] = f["bases"]
  } else if (kind == "kit") {
    if (f["start"] == "" && f["end"] == "")
      fail("the kit " name " has no adapter")
    names[n] = name
    rows[n] = "  {" cstr(name) ", " bases_of(f["start"]) ",\n   " \
              bases_of(f["end"]) "},"
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
    fail("no kit in the file")
  # the names are included once; the table is included by the one
  # translation unit which asks for it, which may be after them
  print "/* Generated from metadata/adapters.yaml by"
  print "   scripts/adapters-table.awk. Do not edit: edit the metadata. */"
  print ""
  print "#ifndef ADAPTERS_TABLE_H"
  print "#define ADAPTERS_TABLE_H"
  print ""
  printf "#define ADAPTERS_N_KITS %d\n", n
  print ""
  list = ""
  comma = ""
  for (i = 0; i < n; i++) {
    list = list (i == 0 ? "" : ", ") names[i]
    comma = comma (i == 0 ? "" : ",") names[i]
  }
  print "/* the kits, for the option --adapters */"
  printf "#define ADAPTERS_KITS %s\n", cstr(list)
  print ""
  print "/* the default: all the kits */"
  printf "#define ADAPTERS_DEFAULT %s\n", cstr(comma)
  print ""
  print "#endif"
  print ""
  print "#if defined(ADAPTERS_TABLE_DEFINE) && !defined(ADAPTERS_TABLE_DEFINED)"
  print "#define ADAPTERS_TABLE_DEFINED"
  print "static const adapters_kit_t ADAPTERS_TABLE[ADAPTERS_N_KITS] = {"
  for (i = 0; i < n; i++)
    print rows[i]
  print "};"
  print "#endif"
}
