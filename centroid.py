#!/usr/bin/env python3
"""
Centroiding a point source in a DS9-displayed FITS image using the 1D-Gaussian algorithm from Photutils.

Steps:
 1. Read off a rough integer (x, y) position of your star from DS9.
 2. Load the same FITS file in Python, extract a small cutout around (x, y), subtract the median background.
 3. Run centroid_1dg on the background-subtracted cutout.
 4. Translate the resulting (x, y) back to full-image coordinates.
 5. Plot the cutout with both centroids and write a DS9 region file.
 
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from photutils.centroids import centroid_1dg
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch


# Path to your FITS file
filename = "./vestafits2june/063011.fits"

# Rough integer coordinates from DS9 (read off by clicking on the star)
# Replace these with the values you saw in DS9. (AKA guess it)
x_guess = 510
y_guess = 503

# Half‐width of the box (in pixels) around (x_guess, y_guess). 
# A 20×20 box corresponds to half_width = 10.
half_width = 10

# 2. Load the FITS image

hdul = fits.open(filename)
data_full = hdul[0].data.astype(float)
hdul.close()

# 3. Define the cutout boundaries

# Convert the guesses to integers
x0 = int(x_guess)
y0 = int(y_guess)

# Compute min/max indices, ensuring they stay within array bounds
y_min = max(0, y0 - half_width)
y_max = min(data_full.shape[0], y0 + half_width)
x_min = max(0, x0 - half_width)
x_max = min(data_full.shape[1], x0 + half_width)

# Extract the 2D cutout array
cutout = data_full[y_min:y_max, x_min:x_max]

# 4. Subtract median background from cutout


median_bkg = np.median(cutout)
data_cut = cutout - median_bkg

# 5. Compute the 1D‐Gaussian centroid


# centroid_1dg returns (x_centroid, y_centroid) relative to the cutout,
# where (0, 0) is the lower-left corner of data_cut.
xc_rel, yc_rel = centroid_1dg(data_cut)

# Translate back to full-image coordinates
x_centroid = xc_rel + x_min
y_centroid = yc_rel + y_min

print(f"Refined centroid (full frame): x = {x_centroid:.2f}, y = {y_centroid:.2f}")

norm = ImageNormalize(data_cut, interval=ZScaleInterval(), stretch=LinearStretch())

plt.figure(figsize=(6, 6))
plt.imshow(data_cut, origin="lower", cmap="viridis", norm=norm)
plt.colorbar(label="(counts − median)")

# Convert the integer‐pixel guess into cutout‐coordinates:
x_guess_cut = x_guess - x_min
y_guess_cut = y_guess - y_min

# Plot the integer‐pixel guess (red "×") and 1D‐Gaussian centroid (blue "+")
plt.scatter(x_guess_cut, y_guess_cut,
            s=100, c="red", marker="x", label="DS9 integer guess")
plt.scatter(xc_rel, yc_rel,
            s=100, c="blue", marker="+", label="1D-Gaussian centroid")

plt.legend(loc="upper right")
plt.xlabel("X (pixels in cutout)")
plt.ylabel("Y (pixels in cutout)")
plt.title("Centroiding around the star")
plt.tight_layout()
plt.show()


reg_text = (
    "# Region file format: DS9 version 4.1\n"
    "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\"\n"
    "image\n"
    f"point({x_centroid:.4f},{y_centroid:.4f})\n"
)

with open("centroid.reg", "w") as f:
    f.write(reg_text)

print("Wrote DS9 region file → centroid.reg")
