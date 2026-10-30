#!/usr/bin/env python3
"""
Write the API reference chapter of the manual (reference.rst) from the
headers of the library: for each header, its description (the first
comment after the copyright line) and its declarations (functions, types,
macros), without the documentation comments.

Usage (from the repository root): docs/manual/make_reference.py
"""
import os
import re

COMPONENTS = [("core", "Core"), ("ds", "Data structures"),
              ("io", "Large files and compression"),
              ("seq", "Sequences and reads")]
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")

def description(text):
  """the comment before the include guard, other than the copyright
  line, as one paragraph"""
  guard = re.search(r"#ifndef \w+_H", text)
  head = text[:guard.start()] if guard else ""
  for c in re.findall(r"/\*(.*?)\*/", head, re.S):
    if c.strip().startswith("(c)"):
      continue
    para = []
    for l in c.strip().splitlines():
      if not l.strip():
        break
      para.append(l.strip())
    return " ".join(para)
  return ""

def declarations(text):
  """the declarations, without comments, preprocessor guards, includes,
  test functions and the internal functions (_ or __ prefix, _ suffix)"""
  text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
  text = re.sub(r"//[^\n]*", "", text)
  out = []
  macros = set()
  lines = text.splitlines()
  i = 0
  while i < len(lines):
    l = lines[i]
    s = l.strip()
    if not s or s.startswith(("#include", "#ifndef", "#endif", "#else",
                              "#if ", "#ifdef")) or \
       re.match(r"#define [A-Z_]+_H\s*$", s):
      i += 1
      continue
    block = [l.rstrip()]
    if s.startswith("#define"):
      # the macro name and parameters only
      full = s
      while lines[i].rstrip().endswith("\\"):
        i += 1
        full = full.rstrip("\\").rstrip() + " " + lines[i].strip()
      i += 1
      m = re.match(r"(#define\s+(\w+)(\([^)]*\))?)", full)
      if m and not m.group(2).startswith("_") and \
         m.group(2) not in macros:
        macros.add(m.group(2))
        out.append(m.group(1))
      continue
    # a declaration or definition: until ";" or the end of a body
    depth = l.count("{") - l.count("}")
    while (depth > 0 or not (block[-1].rstrip().endswith(";") or
           (block[-1].rstrip().endswith("}") and depth == 0))) and \
            i + 1 < len(lines):
      i += 1
      block.append(lines[i].rstrip())
      depth += lines[i].count("{") - lines[i].count("}")
    i += 1
    decl = "\n".join(block)
    name = re.search(r"(\w+)\s*\(", decl)
    if name and (name.group(1).endswith(("_test", "_")) or
                 name.group(1).startswith(("_", "test_"))):
      continue
    if "static inline" in decl:
      decl = decl.split("{")[0].rstrip() + ";"
    out.append(decl)
  return out

def main():
  rst = [".. _reference:", "", "API reference", "-------------", "",
         "The declarations of the headers of the library (generated "
         "from the headers by ``docs/manual/make_reference.py``; the "
         "documentation of each function is in the header, before its "
         "declaration). The unit test functions (``<module>_test()``) "
         "and the internal functions (``_`` or ``__`` prefix, called by "
         "macros, or ``_`` suffix, called by inline functions) are not "
         "listed.", ""]
  for comp, title in COMPONENTS:
    rst += [title + " (``" + comp + "``)",
            "~" * (len(title) + len(comp) + 8), ""]
    d = os.path.join(ROOT, "src", comp)
    for h in sorted(f for f in os.listdir(d) if f.endswith(".h")):
      text = open(os.path.join(d, h)).read()
      rst += ["**" + h + "**", ""]
      desc = description(text)
      if desc:
        rst += [desc, ""]
      decls = declarations(text)
      if decls:
        rst += [".. code-block:: c", ""]
        for decl in decls:
          rst += ["   " + l for l in decl.splitlines()
                  if l.strip()]
        rst += [""]
  path = os.path.join(ROOT, "docs", "manual", "reference.rst")
  open(path, "w").write("\n".join(rst).rstrip() + "\n")

if __name__ == "__main__":
  main()
