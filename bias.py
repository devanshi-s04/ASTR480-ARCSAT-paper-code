#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: bias.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import pathlib
import numpy as np


def create_median_bias(bias_list, median_bias_filename):
    """This function must:

    - Accept a list of bias file paths as bias_list.
    - Read each bias file and create a list of 2D numpy arrays.
    - Use a sigma clipping algorithm to combine all the bias frames using
      the median and removing outliers outside 3-sigma for each pixel.
    - Save the resulting median bias frame to a FITS file with the name
      median_bias_filename.
    - Return the median bias frame as a 2D numpy array.

    """
    # bias_list: Python List of Path objs. 

    # itr over list, read each file, append 2d arrs
    biasImg = []
    for bias in bias_list:
        biasData = fits.getdata(bias)
        biasImg.append(biasData.astype('f4'))
    
    # ligma
    biasImg_masked = sigma_clip(biasImg, cenfunc='median', sigma=3, axis=0)
    median_bias = np.ma.median(biasImg_masked, axis=0)

    # create a new FITS file from the resulting median bias frame.
    primary = fits.PrimaryHDU(data=median_bias.data, header=fits.Header())
    hdul = fits.HDUList([primary])
    hdul.writeto(median_bias_filename, overwrite=True)

    # return img
    return median_bias
