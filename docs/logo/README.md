# Logo

The logos were generated with AI image tools and provided by the author
(to be disclosed where they are used, e.g. in a paper).

| file | what it is |
|---|---|
| `flucut-logo-original.png` | the square logo as received, 1536 x 1024, on an off-white background |
| `flucut-logo.png` | the same with a transparent background, cropped (1260 x 765) |
| `flucut-logo-640.png` | 640 pixels wide (light theme) |
| `flucut-logo-dark-640.png` | the same for dark backgrounds: the navy "flu" recolored to white (derived, not a separate design: the near-white strands become dark) |
| `flucut-banner-original.png` | the wide logo as received, 2172 x 724, on an off-white background |
| `flucut-banner.png` | the same with a transparent background, cropped (1965 x 569) |
| `flucut-banner-640.png` | 640 pixels wide, at the top of the README (light theme) |
| `flucut-banner-dark-original.png` | the wide logo as received, designed on a black background, 2172 x 724 |
| `flucut-banner-dark.png` | the same with a transparent background, cropped |
| `flucut-banner-dark-640.png` | 640 pixels wide, at the top of the README (dark theme) |
| `flucut-banner.ans` | the wide logo as colored text (Unicode Braille characters, 24-bit ANSI colors on black, 77 columns, 9 lines), as received from ASCIILogo.com in a shell script, without the script lines; embedded in the program and shown at the start of a run in a terminal (`src/banner.c`) |

The transparent and reduced versions are made by `make_logos.py` (Python,
Pillow, NumPy; the method is described in the script):

    python3 docs/logo/make_logos.py docs/logo flucut-logo \
        --name-box 0 0.62 1 1
    python3 docs/logo/make_logos.py docs/logo flucut-banner --dark-original

In short: the background (white, or black for the dark original) is
made transparent pixel by pixel ("color to alpha"), so that the glows
stay partly transparent; the image is cropped to the logo and reduced to
640 pixels wide; the README shows the wide logo 480 pixels wide, so
that it stays sharp on high-resolution screens.

The logo of the HTML report (`src/report-logo-data.inc`: the wide logo
and its dark variant, as WebP data URIs, 520 pixels wide) is made by
`scripts/make_report_logo.py` (from Lumipore; Pillow):

    python3 -I scripts/make_report_logo.py docs/logo/flucut-banner.png \
        docs/logo/flucut-banner-dark.png > src/report-logo-data.inc
