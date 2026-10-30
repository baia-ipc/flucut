# Logo

| file | what it is |
|---|---|
| `flucut-logo-original.png` | the logo as received, 1536 x 1024, on an off-white background (AI-generated image, provided by the author) |
| `flucut-logo.png` | the same with a transparent background, cropped to the logo (1260 x 765) |
| `flucut-logo-640.png` | `flucut-logo.png` 640 pixels wide, for the README (light theme) |
| `flucut-logo-dark-640.png` | the same for the dark theme: the navy "flu" recolored to white |

How the transparent and reduced versions were made (Python, Pillow and
NumPy), from `flucut-logo-original.png`:

1. the white background is turned into transparency ("color to
   alpha"): the alpha of a pixel is the largest distance of its
   channels from white, and its color is unblended from white
   accordingly, so that the glows stay partially transparent; alphas
   below 0.02 (the noise of the background) are set to 0;
2. the image is cropped to the pixels with alpha above 8 (of 255),
   with a margin of 16 pixels;
3. it is resized to 640 pixels wide (Lanczos);
4. for the dark variant, the pixels of the lower 38% of the image
   (the name) which are dark and bluish (blue above red by more than
   15, mean of the channels below 140, alpha above 20) are made white,
   keeping their alpha; then it is resized as above.

On a dark background, the near-white strands and the white stripe of
the cutter, which became transparent, show as dark: the dark variant is
derived, not a separate design.
