#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: demo.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
### Danbi thinks this is important, but its 11:13pm and I am tired, this is due in 30 minutes...
### Python Script to run CCD reduction on team data
# I recommend running this by chunks on console or jupyter notebook
# Need some manual file cleanup in between steps
# but if you decide to use this script, make sure to run this code in /src or change the paths
# and to modify this code to reduce the specific data you want


from ecyce import run_reduction
import pathlib

# data directories in /src
dir1 = pathlib.Path('20250515') 
dir2 = pathlib.Path('20250602')

# run reductions
run_reduction(dir1, 'vesta_g') # switch the string to 'markarian421' and 'tcrb_g' if you need

# move median calibration files to a separate directory if you want to keep them
# I have already stored median frames in cal_median subdir corresponding to the date

run_reduction(dir2, 'Vesta1s') # switch to 'Vesta5s' after moving reduced sci images to a seperate directory. There's also one frame each of 2s and 3s exposure if you need

