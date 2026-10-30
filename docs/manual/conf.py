# Configuration of the Sphinx documentation builder for the libcbaia
# manual. Build the PDF with "make manual" in the root of the repository
# (requires Sphinx and a LaTeX installation), or the HTML version with
# "sphinx-build -M html docs/manual docs/manual/_build".

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

project = "libcbaia"
author = "Giorgio Gonnella"
copyright = "2026, Giorgio Gonnella"
version = ""
release = ""

extensions = []
templates_path = []
exclude_patterns = ["_build", "examples"]
master_doc = "index"
source_suffix = ".rst"
highlight_language = "none"
numfig = True
pygments_style = "cbaia_style.CbaiaStyle"

html_theme = "alabaster"
html_title = "libcbaia manual"

latex_engine = "pdflatex"
latex_documents = [
    ("index", "libcbaia.tex", "libcbaia manual", author, "manual"),
]
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "releasename": "",
    "preamble": r"""
\usepackage{enumitem}
\setlistdepth{9}
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{2192}{\ensuremath{\rightarrow}}
""",
    "fncychap": r"\usepackage[Bjornstrup]{fncychap}",
    "printindex": "",
    "extraclassoptions": "openany",
}
latex_show_urls = "footnote"
