# Turn metadata/primers.yaml into the tables of the program.
#
# Writes a header with the named sequences, the built-in primers and the
# built-in schemes (see primers.h and scheme.c). Run by the Makefile; it
# needs nothing but awk, so that building flucut still needs only a C
# compiler and make.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/primers-table.awk \
#          metadata/primers.yaml > out.h

# the macro of the named sequence <s>
function seq_macro(s,   out) {
  out = toupper(s)
  gsub(/[^A-Z0-9]/, "_", out)
  return "PRIMERS_" out
}

# the bases of <list>, a space-separated list of named sequences and
# literal bases
function bases_of(list,   parts, i, k, out, p) {
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
  return out
}

# the segments of <list> as the bits of primer_segments_t
function segments_of(list,   parts, i, k, out, p) {
  if (list == "all")
    return "PRIMER_SEGMENTS_ALL"
  out = ""
  k = split(list, parts, / +/)
  for (i = 1; i <= k; i++) {
    p = parts[i]
    if (p == "")
      continue
    if (p !~ /^(PB2|PB1|PA|HA|NP|NA|M|NS)$/)
      fail("unknown segment " p)
    out = out (out == "" ? "" : " | ") "PRIMER_SEGMENT(SEGMENT_" p ")"
  }
  if (out == "")
    fail("no segment")
  return out
}

function end_of(s) {
  if (s == "left")
    return "PRIMER_END_LEFT"
  if (s == "right")
    return "PRIMER_END_RIGHT"
  fail("the end must be left or right, found \"" s "\"")
}

function virus_of(s) {
  if (s == "A")
    return "VIRUS_A"
  if (s == "B")
    return "VIRUS_B"
  fail("the virus must be A or B, found \"" s "\"")
}

# the names of <list> as C strings, comma-separated
function names_of(list,   parts, i, k, out) {
  k = split(list, parts, / +/)
  out = ""
  for (i = 1; i <= k; i++)
    if (parts[i] != "")
      out = out (out == "" ? "" : ", ") cstr(parts[i])
  if (out == "")
    fail("no primer")
  return out
}

# the record read so far, written out
function flush(   seq, term, at) {
  if (kind == "sequence") {
    if (f["bases"] == "")
      fail("the sequence " name " has no bases")
    bases[name] = f["bases"]
    seq_names[n_seq] = name
    n_seq++
  } else if (kind == "primer") {
    seq = bases_of(f["seq"])
    term = bases_of(f["terminus"])
    at = length(seq) - length(term)
    # the terminus is the end of the oligonucleotide, which is what
    # primer_t.terminus_from says; checked here, where the two are
    # written apart
    if (at < 0 || substr(seq, at + 1) != term)
      fail("the terminus of " name " is not the end of its sequence: " \
           seq " does not end in " term)
    primers[n_prim] = \
      "  {" cstr(name) ", " cstr(seq) ", " length(seq) ", " at ", " \
      length(term) ",\n   " end_of(f["end"]) ", " virus_of(f["virus"]) \
      ", " segments_of(f["segments"]) "},"
    n_prim++
  } else if (kind == "scheme") {
    schemes[n_sch] = \
      "  {" cstr(name) ",\n   {" names_of(f["left"]) ", NULL},\n" \
      "   {" names_of(f["right"]) ", NULL}},"
    n_sch++
  } else if (kind != "") {
    fail("unknown record \"" kind "\"")
  }
  kind = ""
  name = ""
  split("", f)
}

BEGIN { n_seq = 0; n_prim = 0; n_sch = 0 }

END {
  if (n_prim == 0)
    fail("no primer in the file")
  print "/* Generated from metadata/primers.yaml by"
  print "   scripts/primers-table.awk. Do not edit: edit the metadata. */"
  print ""
  # the macros are included once; each table is included by the one
  # translation unit which asks for it, which may be after the macros
  print "#ifndef PRIMERS_TABLE_H"
  print "#define PRIMERS_TABLE_H"
  print ""
  print "/* the conserved termini and the tails the primers are made of */"
  for (i = 0; i < n_seq; i++)
    printf "#define %s %s\n", seq_macro(seq_names[i]),
           cstr(bases[seq_names[i]])
  print ""
  print "#endif"
  print ""
  print "#if defined(PRIMERS_TABLE_DEFINE) && !defined(PRIMERS_TABLE_DEFINED)"
  print "#define PRIMERS_TABLE_DEFINED"
  printf "#define PRIMERS_N_BUILTIN %d\n", n_prim
  print "static const primer_t PRIMERS_BUILTIN[PRIMERS_N_BUILTIN] = {"
  for (i = 0; i < n_prim; i++)
    print primers[i]
  print "};"
  print "#endif"
  print ""
  print "#if defined(SCHEMES_TABLE_DEFINE) && !defined(SCHEMES_TABLE_DEFINED)"
  print "#define SCHEMES_TABLE_DEFINED"
  printf "#define SCHEMES_N_BUILTIN %d\n", n_sch
  print "static const builtin_scheme_t SCHEMES_BUILTIN[SCHEMES_N_BUILTIN] = {"
  for (i = 0; i < n_sch; i++)
    print schemes[i]
  print "};"
  print "#endif"
}
