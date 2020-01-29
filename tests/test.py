
import os
import csv
import astropy
import numpy as np
import pandas as pd
import seaborn as sb

import requests
import lightkurve as lk
from os.path import exists
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from astroplan import FixedTarget
from astropy.stats import BoxLeastSquares
from astroplan.plots import plot_finder_image
from astropy.stats import median_absolute_deviation
from lightkurve import TessTargetPixelFile

# test the downloading of the data
# get the outir (or indir) path
#syspath = str(os.path.abspath(utils.__file__))[0:-14]
#print (syspath)
#
#lcfile = ['https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018319095959-s0005-0000000055525572-0125-s_lc.fits']    
#
#response = requests.get(lcfile[0])
#
## open the file using the response url  
#lchdu  = pf.open(response.url) # this needs to be a URL - not a filefile
#
#lchdu.writeto('/Users/Nora/Documents/research/TESS/planethunters/code/LATTE/tests/tic55525572_lc.fits')
#exit()

#file = 'https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018319095959-s0005-0000000055525572-0125-s_tp.fits'
#
#tpf = pf.open(file)   # open the file
#
#print (tpf)
#print (type(tpf))
#
#tpf.writeto('/Users/Nora/Documents/research/TESS/planethunters/code/LATTE/tests/tic55525572_tp.fits')















