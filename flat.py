#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: flat.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import pathlib
import numpy as np
import matplotlib.pyplot as plt 
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch

def create_median_flat(
    flat_list,
    bias_filename,
    median_flat_filename,
    dark_filename=None,
):
    """This function must:

    - Accept a list of flat file paths to combine as flat_list. Make sure all
      the flats are for the same filter.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in flat_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each flat image.
    - Optionally you can pass a dark frame filename as dark_filename and subtract
      the dark frame from each flat image (remember to scale the dark frame by the
      exposure time of the flat frame).
    - Use a sigma clipping algorithm to combine all the bias-corrected flat frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Create a normalised flat divided by the median flat value.
    - Save the resulting median flat frame to a FITS file with the name
      median_flat_filename.
    - Return the normalised median flat frame as a 2D numpy array.

    """
    # get bias
    bias = fits.getdata(bias_filename).astype('f4')
    # open files, make sure flats are for the same filter. create a list of 2d arrs
    flats = [] # list of 2d arr to store flats data
    # make sure our flats have the g filter
    for path in flat_list:
        flatHeader = fits.getheader(path)
        if flatHeader['FILTER'] == 'g':
            data = fits.getdata(path)
            flats.append(data.astype('f4'))
    # subtract bias from each flat
    flatsMinusBias = flats - bias
    # SIGMAA SIGMAA BOIII
    flatsSC = sigma_clip(flatsMinusBias, cenfunc='median', sigma=3, axis=0)
    # combine and normalize
    median_flat = np.ma.median(flatsSC, axis=0) # combine
    median_flat_trim = median_flat[100:-100, 100:-100] # trim
    median_value = np.ma.median(median_flat_trim) # median value
    median_flat = median_flat / median_value # normalize
    # create and save the resulting median flat frame in fits
    flat_hdu = fits.PrimaryHDU(data=median_flat.data, header=flatHeader)
    flat_hdu.header['COMMENT'] = 'Normalized flat-field image'
    flat_hdu.header['BIASFILE'] = (str(bias_filename), 'Bias image used to subtract bias level')
    flat_hdu.writeto(median_flat_filename, overwrite=True)
    # return 2d arr
    return median_flat


def plot_flat(
    median_flat_filename,
    ouput_filename="median_flat.png",
    profile_ouput_filename="median_flat_profile.png",
):
    """This function must:

    - Accept a normalised flat file path as median_flat_filename.
    - Read the flat file.
    - Plot the flat frame using matplotlib.imshow with reasonable vmin and vmax
      limits. Save the plot to the file specified by output_filename.
    - Take the median of the flat frame along the y-axis. You'll end up with a
      1D array.
    - Plot the 1D array using matplotlib.
    - Save the plot to the file specified by profile_output_filename.

    """
    # get normalized flat
    flat = fits.getdata(median_flat_filename).astype('f4')
    # plot flat using plt imshow
    norm = ImageNormalize(flat, interval=ZScaleInterval(), stretch=LinearStretch())
    plt.imshow(flat, origin='lower', norm=norm, cmap='YlOrBr_r') # cannot set vmin,vmax since a norm instance is passed
    plt.savefig(ouput_filename)
    plt.close()
    # median along y-axis
    profile = np.ma.median(flat, axis=1)
    plt.plot(profile)
    plt.savefig(profile_ouput_filename)
    return None
