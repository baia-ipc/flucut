# (c) 2026 Giorgio Gonnella
# Configuration of the Sphinx builder for the PDF of the methods
# (methods.md in this directory, the only copy of the text). Build it
# with "make methods" in the root of the repository (requires Sphinx with
# myst-parser, matplotlib for the figures, and XeLaTeX).

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pseudocode_lexer import PseudocodeLexer  # noqa: E402

project = "flucut methods"
author = "Giorgio Gonnella"
copyright = "2026, Giorgio Gonnella"

extensions = ["myst_parser"]
master_doc = "methods"
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "figure_scripts"]
myst_heading_anchors = 2
highlight_language = "pseudocode"
pygments_style = "pseudocode_lexer.PseudocodeStyle"
numfig = False


def adapt_markdown(app, docname, source):
  """What the PDF needs differently from the Markdown file: the table
    of contents is generated, and HTML fragments have no meaning."""
  text = source[0]
  # the hand-made contents list (up to the first rule) and the rules
  text = re.sub(r"\nContents\n.*?\n---\n", "\n", text, count=1,
                flags=re.S)
  text = re.sub(r"(?m)^---\n", "", text)
  text = re.sub(r"<sub>(.*?)</sub>", r"{sub}`\1`", text)
  text = text.replace("&nbsp;", " ")
  # the source files of a section and the procedures that call it,
  # block quotes in the Markdown file: as boxes
  text = re.sub(r"(?m)^> \*\*(Source|Called from):\*\* (.*(?:\n> .*)*)",
                lambda m: "```{admonition} %s\n" % m.group(1) +
                m.group(2).replace("\n> ", "\n") + "\n```", text)
  # the figures: the PDF copy of each, with its caption
  text = re.sub(r"!\[(.*?)\]\(figures/(\w+)\.svg\)\n\n\*Figure: (.*?)\*\n",
                r"```{figure} figures/\2.pdf\n:alt: \1\n\n\3\n```\n", text)
  # each numbered section (A1, B2, ...) begins on a new page
  text = re.sub(r"(?m)^(## [A-E]\d+\. )",
                lambda m: "```{raw} latex\n\\clearpage\n```\n\n" +
                m.group(1), text)
  # a block of pseudocode is kept on one page: if the room left on the
  # page is not enough, it begins on the next one (a block longer than
  # a page begins at the top of one)
  def keep_together(m):
    height = min(10.2 * (m.group(0).count("\n") + 1) + 16, 640)
    return ("```{raw} latex\n\\Needspace*{%dpt}\n```\n\n" % height +
            m.group(0))
  text = re.sub(r"(?m)^```\n[A-Z][A-Z0-9_]*\(.*?^```$", keep_together,
                text, flags=re.S)
  # a document which is not part of the PDF: its name, not a link
  text = re.sub(r"\[(methods_\w+\.md)\]\(\1\)", r"`\1`", text)
  source[0] = text


def setup(app):
  app.add_lexer("pseudocode", PseudocodeLexer)
  app.connect("source-read", adapt_markdown)


latex_engine = "xelatex"
latex_documents = [
  (master_doc, "methods.tex", "flucut methods",
   author, "manual"),
]
# title page as that of the flucut manual (docs/manual), with its
# illustration and colours
latex_additional_files = ["../manual/figures/title.jpg"]
latex_maketitle = r"""
\begin{titlepage}
\begin{tikzpicture}[remember picture, overlay]
  \node[anchor=north west, inner sep=0]
    at ($(current page.north west)+(-5.6,0)$)
    {\includegraphics[height=18cm]{title.jpg}};
  \fill[white, path fading=east] (current page.north west) rectangle
    ($(current page.north west)+(9.5,-18)$);
  \fill[white, path fading=north]
    ($(current page.north west)+(0,-18.05)$) rectangle
    ($(current page.north east)+(0,-14)$);
  \begin{scope}[every node/.style={anchor=west, inner sep=0, align=left},
                 font=\titlefont]
    \node[text=fcblue, font=\titlefont\bfseries\small]
      at ($(current page.north west)+(2,-12.6)$)
      {ALGORITHMS \textbullet\ PSEUDOCODE};
    \fill[fcorange] ($(current page.north west)+(2,-13.15)$) rectangle
      ++(1.3,0.07);
    \node[font=\titlefont\bfseries\fontsize{52}{56}\selectfont]
      at ($(current page.north west)+(1.85,-14.5)$)
      {\textcolor{fcnavy}{flu}\textcolor{fcorange}{cut}
       \textcolor{fcnavy}{methods}};
    \node[text=fcink, font=\titlefont\LARGE]
      at ($(current page.north west)+(2,-16.3)$)
      {Technical description of the algorithms};
    \node[anchor=south west, text=fcnavy, font=\titlefont\large]
      at ($(current page.south west)+(2,2.5)$)
      {\textbf{@AUTHOR@}\\[0.2em]
       {\color{fcmuted}\normalsize BAIA, Institut Pasteur du Cambodge}\\[1.2em]
       {\color{fcmuted}\small \today}};
  \end{scope}
\end{tikzpicture}
\end{titlepage}
""".replace("@AUTHOR@", author)

latex_elements = {
  "papersize": "a4paper",
  "pointsize": "11pt",
  # a figure stays where the text has it
  "figure_align": "H",
  "extraclassoptions": "openany,oneside",
  # the longest lines of pseudocode have about 120 characters
  "sphinxsetup": "hmargin={1.8cm,1.8cm}, vmargin={2.3cm,2.3cm}, "
                 "verbatimwithframe=false, VerbatimColor={HTML}{F7F8FA}, "
                 # the boxes of the source files and of the callers of a section
                 "div.note_border-width=0.4pt, div.note_border-radius=3pt, "
                 "div.note_border-TeXcolor={HTML}{C9CED6}, "
                 "div.note_background-TeXcolor={HTML}{F3F6FB}, "
                 "div.note_padding=4pt, div.note_TeXextras=\\footnotesize, "
                 "div.note_title-background-TeXcolor={HTML}{EEF1F5}, "
                 "div.note_title-foreground-TeXcolor={HTML}{03215A}",
  "preamble": r"""
\usepackage{needspace}
\usepackage{enumitem}
\setlistdepth{9}
\usepackage{tikz}
\usetikzlibrary{calc,fadings}
\definecolor{fcnavy}{HTML}{03215A}
\definecolor{fcorange}{HTML}{F84A2A}
\definecolor{fcblue}{HTML}{0A72F5}
\definecolor{fcink}{HTML}{2B2B2B}
\definecolor{fcmuted}{HTML}{6B7280}
\newfontfamily\titlefont{FiraSans}[Extension=.otf,
  UprightFont=*-Regular, BoldFont=*-Bold]
% a monospaced font with the mathematical symbols of the pseudocode
\setmonofont{DejaVuSansMono}[Extension=.ttf, BoldFont=*-Bold,
  ItalicFont=*-Oblique, BoldItalicFont=*-BoldOblique, Scale=0.85]
\fvset{fontsize=\fontsize{7.8}{10.2}\selectfont}
% (a symbol which that font lacks)
\usepackage{newunicodechar}
\newfontfamily\symbolfallback{FreeMono}[Extension=.otf, Scale=0.85]
\newunicodechar{ℓ}{{\symbolfallback ℓ}}
% the sections carry their own labels (A1, B2, ...): no numbers
\setcounter{secnumdepth}{-1}
% room for a procedure name in a narrow column of a table
\setlength{\tymin}{105pt}
""",
  "maketitle": latex_maketitle,
  "fncychap": r"\usepackage[Bjornstrup]{fncychap}",
  "printindex": "",
}
