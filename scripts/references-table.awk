# Turn metadata/references.yaml into the lists of the program.
#
# Writes a header with, for each set of references, its version and the
# accessions of the set as one comma-separated string (see
# builtin-references.c). Run by the Makefile.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/references-table.awk \
#          metadata/references.yaml > out.h

function flush(   parts, i, k, p) {
  if (kind == "set") {
    if (f["version"] == "")
      fail("the set " name " has no version")
    if (name !~ /^[A-Z]$/)
      fail("a set is named by the type of influenza, found \"" name "\"")
    set = name
    sets[n_sets] = name
    version[name] = f["version"]
    n_sets++
  } else if (kind == "group") {
    if (set == "")
      fail("the group " name " is before any set")
    k = split(f["accessions"], parts, / +/)
    for (i = 1; i <= k; i++) {
      p = parts[i]
      if (p == "")
        continue
      if (p !~ /^[A-Z]+[0-9]+(\.[0-9]+)?$/)
        fail("\"" p "\" is not an accession")
      if (p in seen)
        fail("the accession " p " is listed twice")
      seen[p] = 1
      accessions[set] = accessions[set] (accessions[set] == "" ? "" : ",") p
      n_acc[set]++
    }
  } else if (kind != "") {
    fail("unknown record \"" kind "\"")
  }
  kind = ""
  name = ""
  split("", f)
}

BEGIN { n_sets = 0; set = "" }

END {
  if (n_sets == 0)
    fail("no set in the file")
  print "/* Generated from metadata/references.yaml by"
  print "   scripts/references-table.awk. Do not edit: edit the metadata. */"
  print ""
  print "#ifndef REFERENCES_TABLE_H"
  print "#define REFERENCES_TABLE_H"
  for (i = 0; i < n_sets; i++) {
    s = sets[i]
    if (n_acc[s] == 0)
      fail("the set " s " has no accession")
    print ""
    printf "/* %d accessions */\n", n_acc[s]
    printf "#define BUILTIN_REFERENCES_%s_VERSION %s\n", s, cstr(version[s])
    printf "static const char *const BUILTIN_ACCESSIONS_%s =\n  %s;\n", s,
           cstr(accessions[s])
  }
  print ""
  print "#endif"
}
