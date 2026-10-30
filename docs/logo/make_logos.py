# (c) 2026 Giorgio Gonnella
"""Make the transparent, reduced and dark versions of the logos.

From <name>-original.png (an image on a white background) make, in the
same directory:
  <name>.png            the white background made transparent ("color
                        to alpha": the alpha of a pixel is the largest
                        distance of its channels from white, its color
                        unblended from white, so the glows stay partly
                        transparent; alphas below 0.02 set to 0),
                        cropped to the pixels with alpha above 8 (of
                        255), with a margin of 16 pixels
  <name>-640.png        the same, 640 pixels wide (Lanczos)
  <name>-dark-640.png   for dark backgrounds: with --dark-original,
                        from <name>-dark-original.png (the logo designed
                        on a black background), its black made
                        transparent in the same way (the alpha of a pixel
                        is its largest channel), cropped and reduced;
                        otherwise the navy pixels of the name (dark and
                        bluish: blue above red by more than 15, mean of
                        the channels below 140, alpha above 20) made
                        white, keeping their alpha, in the region given
                        by --name-box
  <name>-dark.png       with --dark-original, the full-size version

Usage:
  python3 docs/logo/make_logos.py <dir> <name> \
      (--name-box X0 Y0 X1 Y1 | --dark-original)
(the box as fractions of the width and height of the cropped image)
"""

import argparse
import os

import numpy as np
from PIL import Image

WIDTH = 640


def to_alpha(path, background=255):
  """The image of <path>, its white (or black: <background> 0)
  background made transparent and cropped (see the module
  documentation)."""
  im = np.asarray(Image.open(path).convert("RGB")).astype(float)
  alpha = np.abs(background - im).max(axis=2) / 255.0
  safe = np.where(alpha > 0, alpha, 1)
  rgb = np.clip((im - background * (1 - alpha[..., None]))
                / safe[..., None], 0, 255)
  alpha = np.where(alpha < 0.02, 0, alpha)
  img = Image.fromarray(np.dstack([rgb, alpha * 255]).astype(np.uint8),
                        "RGBA")
  box = img.getchannel("A").point(lambda v: 255 if v > 8 else 0).getbbox()
  m = 16
  box = (max(box[0] - m, 0), max(box[1] - m, 0),
         min(box[2] + m, img.width), min(box[3] + m, img.height))
  return img.crop(box)


def reduced(img):
  return img.resize((WIDTH, round(img.height * WIDTH / img.width)),
                    Image.LANCZOS)


def dark(img, name_box):
  """The navy pixels of the name in <name_box> made white."""
  d = np.asarray(img).astype(float)
  h, w = d.shape[:2]
  x0, y0, x1, y1 = name_box
  region = np.zeros((h, w), dtype=bool)
  region[int(y0 * h):int(y1 * h), int(x0 * w):int(x1 * w)] = True
  r, g, b, a = d[..., 0], d[..., 1], d[..., 2], d[..., 3]
  navy = region & (b > r + 15) & (r + g + b < 3 * 140) & (a > 20)
  d[navy, 0:3] = 255
  return Image.fromarray(d.astype(np.uint8), "RGBA")


def main():
  parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
  parser.add_argument("dir")
  parser.add_argument("name")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("--name-box", type=float, nargs=4)
  group.add_argument("--dark-original", action="store_true")
  args = parser.parse_args()
  base = os.path.join(args.dir, args.name)
  img = to_alpha(base + "-original.png")
  img.save(base + ".png", optimize=True)
  reduced(img).save(base + "-640.png", optimize=True)
  if args.dark_original:
    dimg = to_alpha(base + "-dark-original.png", background=0)
    dimg.save(base + "-dark.png", optimize=True)
  else:
    dimg = dark(img, args.name_box)
  reduced(dimg).save(base + "-dark-640.png", optimize=True)
  print(base, img.size)


if __name__ == "__main__":
  main()
