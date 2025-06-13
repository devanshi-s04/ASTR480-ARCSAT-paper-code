#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: science.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import pathlib
import numpy as np
from astroscrappy import detect_cosmics

def reduce_science_frame(
    science_filename,
    median_bias_filename,
    median_flat_filename,
    median_dark_filename,
    reduced_science_filename="reduced_science.fits",
):
    """This function must:

    - Accept a science frame filename as science_filename.
    - Accept a median bias frame filename as median_bias_filename (the one you created
      using create_median_bias).
    - Accept a median flat frame filename as median_flat_filename (the one you created
      using create_median_flat).
    - Accept a median dark frame filename as median_dark_filename (the one you created
      using create_median_dark).
    - Read all files.
    - Subtract the bias frame from the science frame.
    - Subtract the dark frame from the science frame. Remember to multiply the
      dark frame by the exposure time of the science frame. The exposure time can
      be found in the header of the FITS file.
    - Correct the science frame using the flat frame.
    - Optionally, remove cosmic rays.
    - Save the resulting reduced science frame to a FITS file with the filename
      reduced_science_filename.
    - Return the reduced science frame as a 2D numpy array.

    """
    # get sci, bias, flat, dark
    bias = fits.getdata(median_bias_filename).astype('f4')
    flat = fits.getdata(median_flat_filename).astype('f4')
    dark = fits.getdata(median_dark_filename).astype('f4')
    sci = fits.getdata(science_filename).astype('f4')
    sciHeader = fits.getheader(science_filename) # for later

    # subtract bias from sci
    sci -= bias
    # subtract dark from sci
    dark *= sciHeader['EXPTIME'] # scale dark current to sci exptime
    sci -= dark
    # correction with flat
    sci /= flat

    # (optional) remove cosmic rays
    mask, reduced_science = detect_cosmics(sci)
    # save fits
    sci_hdu = fits.PrimaryHDU(data=reduced_science, header=sciHeader)
    sci_hdu.header['COMMENT'] = 'Reduced science image'
    sci_hdu.header['BIASFILE'] = (median_bias_filename, 'Bias image used to subtract bias level')
    sci_hdu.header['DARKFILE'] = (median_dark_filename, 'Dark image used to subtract dark current')
    sci_hdu.header['FLATFILE'] = (median_flat_filename, 'Flat image used for FPN correction')
    mask_hdu = fits.ImageHDU(data=mask.astype('uint8'), name='COSMICRAY_MASK') #secondary HDU
    hdul = fits.HDUList([sci_hdu, mask_hdu])
    hdul.writeto(reduced_science_filename, overwrite=True)
    # return cleaned data
    return reduced_science
