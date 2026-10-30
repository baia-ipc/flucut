"""
Pygments style of the manual: the style of Sphinx, but with the
preprocessor directives (e.g. the #define of the macros in the API
reference) as declarations, not as comments. In the LaTeX output, the style
of a token includes that of its parent types (Comment.Preproc: also that of
Comment), thus the comments are not in italic either (only colored).
"""
from pygments.token import Comment
from sphinx.pygments_styles import SphinxStyle

class CbaiaStyle(SphinxStyle):
  styles = dict(SphinxStyle.styles)
  styles[Comment] = "#408090"
  styles[Comment.Preproc] = "#06287e"
