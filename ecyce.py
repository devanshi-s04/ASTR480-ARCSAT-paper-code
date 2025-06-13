#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: ecyce.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import pathlib
from bias import create_median_bias
from dark import create_median_dark
from flat import create_median_flat, plot_flat
from science import reduce_science_frame
from astropy.stats import sigma_clipped_stats


def run_reduction(myDir, sciName):
    """This function must run the entire CCD reduction process. You can implement it
    in any way that you want but it must perform a valid reduction for the two
    science frames in the dataset using the functions that you have implemented in
    this module. Then perform aperture photometry on at least one of the science
    frames, using apertures and sky annuli that make sense for the data.

    No specific output is required but make sure the function prints/saves all the
    relevant information to the screen or to a file, and that any plots are saved to
    PNG or PDF files.

    """
    ############ PARAMS ############
    # myDir: directory (of each obs date)
    ### myDir/cal: calibration frames
    ### myDir/sci: science images
    ### mydir/sci_reduced: directory to store reduced images
    # sciName: str prefix of sci files one wishes to reduce
    
    cal_path = myDir / 'cal'
    sci_path = myDir / 'sci'
    red_path = myDir / 'sci_reduced'
    
    # take a set of bias, create median bias
    bias_files = sorted(cal_path.glob("Bias*.fit*")) # list of paths
    median_bias = create_median_bias(bias_files, 'median_bias.fits') # saves .fits to current dir. median_bias is a 2d array
    
    # take a set of darks and create median dark
    dark_files = sorted(cal_path.glob("Dark*.fit*"))
    median_dark = create_median_dark(dark_files, 'median_bias.fits','median_dark.fits')
    
    # take a set of flats and create median flat
    flat_files = sorted(cal_path.glob("domeflat_g*.fit*"))
    median_flat = create_median_flat(flat_files, 'median_bias.fits', 'median_flat.fits')
    # plot flat
    plot_flat('median_flat.fits') # flat img and profile saved with default filenames
    
    # get science images and reduce
    sci_files = sorted(sci_path.glob(f"{sciName}*.fit*"))
    for j, sci in enumerate(sci_files):
        _ = reduce_science_frame(sci_files[j], 'median_bias.fits', 'median_flat.fits', 'median_dark.fits', reduced_science_filename= red_path/f'reduced_science_{j:02}.fits') # saved files in same order as input sci imgs
    
    # LATER! calculate gain & readout noise
    # gain = ptc.calculate_gain(flat_files[0:2]) # picking the first two just for convenience
    # ron = ptc.calculate_readout_noise(bias_files[0:2], gain)
    # display
    #print(f'Estimated Gain:\t{gain:.2f}\nEstimated Readout Noise [e-]:\t{ron:.2f}')

    #print("Done!")
    
    pass
