# (c) 2026 Giorgio Gonnella
"""Highlighting of the pseudocode blocks: procedure names, step numbers,
keywords and comments."""

from pygments.lexer import RegexLexer, bygroups
from pygments.style import Style
from pygments.token import (Comment, Keyword, Name, Number, Text,
                            Whitespace)

KEYWORDS = (r"for each|elseif|while|return|else|for|and|not|if|or|"
            r"NIL|TRUE|FALSE")


class PseudocodeLexer(RegexLexer):
  name = "pseudocode"
  aliases = ["pseudocode"]
  tokens = {
      "root": [
          # header of a procedure: NAME(arguments)
          (r"^([A-Z][A-Z0-9_]*)(\(.*\))$",
           bygroups(Name.Function, Text)),
          # step number, at the start of a line
          (r"^( *)(\d+)(?=  )", bygroups(Whitespace, Number)),
          (r"//.*$", Comment),
          (r"\b(%s)\b" % KEYWORDS, Keyword),
          # (the "to" of a counted loop)
          (r"(?<= )to(?= )(?=[^\n]*$)(?<=\d to)|(?<=[0-9a-zA-Z] )to"
           r"(?= [^p])", Keyword),
          (r"\s+", Whitespace),
          (r".", Text),
      ],
  }


class PseudocodeStyle(Style):
  """Few colours: those of the title page."""
  background_color = "#F7F8FA"
  styles = {
      Comment: "#6B7280",
      Keyword: "bold #03215A",
      Name.Function: "bold #F84A2A",
      Number: "#6B7280",
  }
