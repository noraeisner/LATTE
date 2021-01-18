# -----------------------------
'''

This test case is carried out with a pre-downloaded lighcurve file from TIC 55525572 Sector 5.

Testing whether the plots are produced as expected. 

python -m unittest tests/test_plotting.py
'''

import os
import sys
import time 
import datetime
import unittest
from dateutil import parser


import warnings
warnings.filterwarnings("ignore")


import LATTE.LATTEutils as utils

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(utils.__file__))[0:-14]

indir = "./LATTE/test_output"		


# -------------
# test with these input parameters (this data is already downloaded so it's only testing that the plotting part works)
tic = '55525572'
sector = '5'
sectors = [5]
transit_list = [1454.7]
transit_sec = '5'
# -------------

# create a mock 'argparser' becasuse the argparse function doesn't work within the unittest module 
# for now use all of the default values. 
class Namespace:
	def __init__(self, **kwargs):  # arbitrary number of input values
		self.__dict__.update(kwargs)

global args
args = Namespace(new_data=False, tic='no',sector='no', targetlist='no', 
	noshow=True, o=False, auto=False,nickname='noname', FFI=False, save=True,north=False,new_path=False,mpi=False)

# test it with the default arguments (except noshow - we don't want the plots to show up in the test run)
# --------------

class TestDataPlotting(unittest.TestCase):
	
	def test_plot(self):

		# download data needed to generate the plots
		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = utils.download_data(indir,sector, tic, binfac = 5, test = './LATTE/tests/tic55525572_lc.fits')
		X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = utils.download_tpf_mast(indir, transit_sec, transit_list, tic, test = './LATTE/tests/tic55525572_tp.fits')
		
		# this function downloads the data and saved plots that show the aperture sizes (for the large and small apertures)
		TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list = utils.download_tpf_lightkurve(indir, transit_list, sectors, tic, test = './LATTE/tests/tic55525572_tp.fits')
		# ---------

		# create the plots that the normal routine makes
		utils.plot_full_md(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, transit_list, args)
		utils.plot_centroid(tic, indir, alltime12, allx1, ally1, allx2, ally2, transit_list, args)
		utils.plot_background(tic, indir, alltime, allfbkg, transit_list, args)
		utils.plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, transit_list, args)
		utils.plot_aperturesize(tic, indir, TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, transit_list, args)

		# --------
		# now check that the plots were actually made!

		# in order to be able to run this code multiple times, assert that the files were created in the last minuet - if the code takes lomger than that to run something is goimg wrong...?
		# without checking when the file is made, the files would have to be manually deleted after every test

		# these are the paths of the files that should have bee created
		full_LC_path = '{}/55525572/55525572_fullLC_md.png'.format(indir)
		full_centroid = '{}/55525572/55525572_centroids.png'.format(indir)
		bkg_path = '{}/55525572/55525572_background.png'.format(indir)
		pixel_path = '{}/55525572/55525572_individual_pixel_LCs_0.png'.format(indir)

		ap_LC_path = '{}/55525572/55525572_aperture_size.png'.format(indir)
		apertures_path = '{}/55525572/55525572_apertures_0.png'.format(indir)

		# Get file's Last modification time stamp only in terms of seconds since epoch 
		time_created_full_LC= os.path.getmtime(full_LC_path)
		time_created_centroid = os.path.getmtime(full_centroid)
		time_created_bkg = os.path.getmtime(bkg_path)
		time_created_pixel = os.path.getmtime(pixel_path)
		 
		time_created_ap_LC = os.path.getmtime(ap_LC_path)
		time_created_apertures = os.path.getmtime(apertures_path)

		# Convert seconds since epoch to readable timestamp
		t_create_full_LC = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_full_LC)))
		t_create_centroid = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_centroid)))		
		t_create_bkg = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_bkg)))
		t_create_pixel = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_pixel)))
		
		t_create_ap_LC = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_ap_LC)))
		t_create_apertures = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_apertures)))
				
		# the time now (needed to get the time since creation)
		t_now = (datetime.datetime.now())

		# -------
		# time difference in minutes
		time_since_creation_full_LC =  ((t_now - t_create_full_LC).seconds / 60)
		time_since_creation_centroid =  ((t_now - t_create_centroid).seconds / 60) 
		time_since_creation_bkg =  ((t_now - t_create_bkg).seconds / 60) 
		time_since_creation_pixel =  ((t_now - t_create_pixel).seconds / 60) 

		time_since_creation_bkg =  ((t_now - t_create_ap_LC).seconds / 60) 
		time_since_creation_pixel =  ((t_now - t_create_apertures).seconds / 60) 

		# check that the expected files were created less than a minute ago
		self.assertLess(time_since_creation_full_LC, 1, "No (new) full LC plot was made in the last five minutes") # a less than b
		self.assertLess(time_since_creation_centroid, 1, "No (new) centroid plot was made in the last five minutes") # a less than b
		self.assertLess(time_since_creation_bkg, 1, "No (new) background plot was made in the last five minutes") # a less than b
		self.assertLess(time_since_creation_pixel, 1, "No (new) pixel level LC plot was made in the last five minutes") # a less than b

		self.assertLess(time_since_creation_bkg, 1, "No (new) background plot was made in the last five minutes") # a less than b
		self.assertLess(time_since_creation_pixel, 1, "No (new) pixel level LC plot was made in the last five minutes") # a less than b

if __name__ == '__main__':

	unittest.main()



