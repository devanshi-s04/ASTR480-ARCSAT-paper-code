#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: dark.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import pathlib

def create_median_dark(dark_list, bias_filename, median_dark_filename):
    """This function must:

    - Accept a list of dark file paths to combine as dark_list.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in dark_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each dark image.
    - Divide each dark image by its exposure time so that you get the dark current
      per second. The exposure time can be found in the header of the FITS file.
    - Use a sigma clipping algorithm to combine all the bias-corrected dark frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Save the resulting dark frame to a FITS file with the name median_dark_filename.
    - Return the median dark frame as a 2D numpy array.

    """
    # read all img in dark_list, create a list of 2d arrs + save header
    headerList = [] # list to store the headers
    darkImg = [] # list to store 2d arrs
    darkCurrent = [] # list to store dark current values (dark / exptime)
     # lorfoop
    for path in dark_list:
        darkData = fits.getdata(path)
        darkImg.append(darkData.astype('f4'))
        header = fits.getheader(path)
        headerList.append(header) 
    # get bias
    bias = fits.getdata(bias_filename).astype('f4') 
    # dark - bias for each dark
    darkMinusBias = darkImg - bias # array arithmetic works?
    # dark current per sec
    darkBiasData = [] # list to store dark/exptime
    for darkBias, header in zip(darkMinusBias, headerList):
        darkBiasData.append(darkBias / header['EXPTIME'])
    # ligma
    darkSC = sigma_clip(darkBiasData, cenfunc='median', sigma=3, axis=0)
    # combine 
    median_dark = np.ma.median(darkSC, axis=0)
    # Save resulting dark frame to {median_dark_filename}.fits
    # using the last header from saved list
    dark_hdu = fits.PrimaryHDU(data=median_dark.data, header=headerList[-1])
    dark_hdu.header['EXPTIME'] = 1
    dark_hdu.header['COMMENT'] = 'Combined dark image with bias subtracted'
    dark_hdu.header['BIASFILE'] = (bias_filename, 'Bias image used to subtract bias level')
    dark_hdu.writeto(median_dark_filename, overwrite=True)
    # return 2d arr
    return median_dark
