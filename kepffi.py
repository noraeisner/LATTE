from __future__ import print_function, absolute_import, division

import os
import astropy
import numpy as np
import pandas as pd
import seaborn as sb

import requests
import lightkurve as lk
from os.path import exists
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from astropy.stats import BoxLeastSquares

from glob import glob
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astroquery.mast import Catalogs
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.patches import Rectangle
from lightkurve import TessTargetPixelFile

from tess_stars2px import tess_stars2px_function_entry
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons

tic = '55525572'
indir = "./LATTE_output"
sectors_all = [2,3,4,5]
sectors = [8]
ra = 72.69401
dec = -60.905461
noshow = False


from argparse import ArgumentParser


if __name__ == '__main__':
    ap = ArgumentParser(description='Lightcurve Analysis Tool for Transiting Exoplanets')
    ap.add_argument('--new-data', action='store_true')
    ap.add_argument('--tic', type=str, help='the target TIC id, if blank will be prompted to enter', default = 'no')
    ap.add_argument('--sector', type=str, help='the sector(s) to look at', default = 'no')
    ap.add_argument('--targetlist', type=str, help='the link to the target file list', default='no')
    ap.add_argument('--noshow', action='store_false', help='if you want to NOT show the plots write --noshow in the command line')
    ap.add_argument('--o', action='store_true', help='if you call this old files will be overwriten in the non-interactive version')
    ap.add_argument('--auto', action='store_false', help='automatic aperture selection')
    ap.add_argument('--nickname', type=str, help='give the target a memorable name', default='no')
    ap.add_argument('--FFI', action='store_true', help='is this an FIIs?')
    ap.add_argument('--north', action='store_false', help='write "north" if you want all the images to be aligned North')
    
    args = ap.parse_args()

    args.append()