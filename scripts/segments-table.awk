# Turn metadata/segments.yaml into the tables of the program.
#
# Writes a header with the number of segments, their gene names and the
# words which name each of them (see segments.h). Run by the Makefile.
#
# Usage: awk -f scripts/metadata-records.awk -f scripts/segments-table.awk \
#          metadata/segments.yaml > out.h

function flush(   parts, i, k, p) {
  if (kind == "segment") {
    if (f["number"] + 0 != n + 1)
      fail("the segments must be numbered in order, " name " is " \
           f["number"] " and not " (n + 1))
    genes[n] = name
    n++
    k = split(f["words"], parts, / +/)
    if (k == 0)
      fail("the segment " name " has no word")
    for (i = 1; i <= k; i++) {
      p = parts[i]
      if (p == "")
        continue
      if (p in seen)
        fail("the word " p " names more than one segment")
      seen[p] = 1
      words[n_words] = "  {" cstr(p) ", " (f["number"] + 0) "},"
      n_words++
    }
    # a name written with a hyphen or a space: the word after it is
    # checked as well, which a negative segment asks for
    if (f["split"] != "") {
      words[n_words] = "  {" cstr(f["split"]) ", -" (f["number"] + 0) "},"
      n_words++
    }
  } else if (kind != "") {
    fail("unknown record \"" kind "\"")
  }
  kind = ""
  name = ""
  split("", f)
}

BEGIN { n = 0; n_words = 0 }

END {
  if (n == 0)
    fail("no segment in the file")
  print "/* Generated from metadata/segments.yaml by"
  print "   scripts/segments-table.awk. Do not edit: edit the metadata. */"
  print ""
  print "#ifndef SEGMENTS_TABLE_H"
  print "#define SEGMENTS_TABLE_H"
  print ""
  print "/* the segments, numbered as in influenza A */"
  printf "#define SEGMENTS_N %d\n", n
  for (i = 0; i < n; i++) {
    name = genes[i]
    gsub(/[^A-Za-z0-9]/, "_", name)
    printf "#define SEGMENT_%s %d\n", toupper(name), i + 1
  }
  print ""
  print "#endif"
  print ""
  print "#if defined(SEGMENTS_TABLE_DEFINE) && !defined(SEGMENTS_TABLE_DEFINED)"
  print "#define SEGMENTS_TABLE_DEFINED"
  print "static const char *const SEGMENTS_GENES[SEGMENTS_N] = {"
  line = " "
  for (i = 0; i < n; i++)
    line = line " " cstr(genes[i]) ","
  print line
  print "};"
  print ""
  print "/* the words which name the gene of a segment, a negative"
  print "   segment asking for the word after it to be checked too */"
  printf "#define SEGMENTS_N_GENE_WORDS %d\n", n_words
  print "static const segment_word_t"
  print "    SEGMENTS_GENE_WORDS[SEGMENTS_N_GENE_WORDS] = {"
  for (i = 0; i < n_words; i++)
    print words[i]
  print "};"
  print "#endif"
}
