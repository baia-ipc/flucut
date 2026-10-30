# (c) 2026 Giorgio Gonnella
# Configuration of the Sphinx documentation builder for the flucut
# manual. Build the PDF with "make manual" in the root of the repository
# (requires Sphinx and a LaTeX installation), or the HTML version with
# "sphinx-build -M html docs/manual docs/manual/_build".

project = "flucut"
author = "Giorgio Gonnella"
copyright = "2026, Giorgio Gonnella"
# the version of the program (src/options.h)
version = "0.1"
release = "0.1"

extensions = []
templates_path = []
exclude_patterns = ["_build"]
master_doc = "index"
source_suffix = ".rst"
highlight_language = "none"
pygments_style = "sphinx"

html_theme = "alabaster"
html_title = "flucut manual"

latex_engine = "pdflatex"
latex_documents = [
  ("index", "flucut.tex", "flucut manual", author, "manual"),
]
# title page: the illustration (white on the left) over the upper part of
# the page, faded into the page, with the title over its left side; the
# colors are those of the logo (docs/logo)
latex_additional_files = ["figures/title.jpg"]
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
                 font=\fontfamily{FiraSans-TLF}\selectfont]
    \node[text=fcblue, font=\fontfamily{FiraSans-TLF}\fontseries{b}%
           \selectfont\small]
      at ($(current page.north west)+(2,-12.6)$)
      {USER MANUAL \textbullet\ RELEASE @RELEASE@};
    \fill[fcorange] ($(current page.north west)+(2,-13.15)$) rectangle
      ++(1.3,0.07);
    \node[font=\fontfamily{FiraSans-TLF}\fontseries{b}%
           \selectfont\fontsize{52}{56}\selectfont]
      at ($(current page.north west)+(1.85,-14.5)$)
      {\textcolor{fcnavy}{flu}\textcolor{fcorange}{cut}};
    \node[text=fcink, font=\fontfamily{FiraSans-TLF}\selectfont\LARGE]
      at ($(current page.north west)+(2,-16.3)$)
      {Primer trimming for influenza A\\[0.1em]before the assembly with IRMA};
    \node[text=fcmuted, font=\fontfamily{FiraSans-TLF}\selectfont\large]
      at ($(current page.north west)+(2,-18)$)
      {Uni12/Uni13 amplicons, internal priming};
    \node[anchor=south west, text=fcnavy,
           font=\fontfamily{FiraSans-TLF}\selectfont\large]
      at ($(current page.south west)+(2,2.5)$)
      {\textbf{@AUTHOR@}\\[0.2em]
       {\color{fcmuted}\normalsize BAIA, Institut Pasteur du Cambodge}\\[1.2em]
       {\color{fcmuted}\small \today}};
  \end{scope}
\end{tikzpicture}
\end{titlepage}
""".replace("@RELEASE@", release).replace("@AUTHOR@", author)

latex_elements = {
  "papersize": "a4paper",
  "pointsize": "11pt",
  "preamble": r"""
\usepackage{enumitem}
\setlistdepth{9}
\usepackage{pmboxdraw}
\usepackage{tikz}
\usetikzlibrary{calc,fadings}
\definecolor{fcnavy}{HTML}{03215A}
\definecolor{fcorange}{HTML}{F84A2A}
\definecolor{fcblue}{HTML}{0A72F5}
\definecolor{fcink}{HTML}{2B2B2B}
\definecolor{fcmuted}{HTML}{6B7280}
""",
  "maketitle": latex_maketitle,
  "fncychap": r"\usepackage[Bjornstrup]{fncychap}",
  "printindex": "",
}
latex_show_urls = "footnote"
