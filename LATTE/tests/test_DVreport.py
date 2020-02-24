# -----------------------------
'''
A test to see whether the DV reports are generated as we would expect them to be. 
This uses previously generated input files and therfore doesn't test the generatino of the files, but only the compilation of the report. 
'''

import os
import sys
import time 
import datetime
import unittest

from dateutil import parser

import warnings
warnings.filterwarnings("ignore")

#custom modules to import
from LATTE import LATTE_DV as ldv

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(ldv.__file__))[0:-14]

indir = "./LATTE/test_output"		

# -------------
# test with these input parameters (this data is already downloaded so it's only testing that the plotting part works)
tic = '55525572'
sector = '5'
sectors = [5]
sectors_all = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13]
transit_list = [1454.7]
target_ra = 72.6941
target_dec = -60.9055

tessmag = 9.8200
teff  = 5824
srad  = 1.9325
tpf_corrupt  = False
astroquery_corrupt = False
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

class TestDVreport(unittest.TestCase):
	
	def test_DVreport(self):

		ldv.LATTE_DV(tic, indir, syspath, transit_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, [0], [0], tpf_corrupt, astroquery_corrupt, FFI = False,  bls = False, model = False, mpi = None, test = './LATTE/tests/')

		# now that it has been run, test to make sure that the report was created. 
		# check that a NEW report was made, and not just an old one that exists from a previous test run 
		DV_path = '{}/55525572/DV_report_55525572.pdf'.format(indir)

		# Get file's Last modification time stamp only in terms of seconds since epoch 
		time_created_DV = os.path.getmtime(DV_path)

		# Convert seconds since epoch to readable timestamp
		t_create_DV = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_DV)))
		
		# the time now (needed to get the time since creation)
		t_now = (datetime.datetime.now())

		# -------
		# time difference in minutes
		time_since_creation_DV =  ((t_now - t_create_DV).seconds / 60)

		# make sure that the new DV report was made less than 0.1 minutes ago - we want to make sure it's a new file and not am old one.
		self.assertLess(time_since_creation_DV, 0.1, "No full LC plot was generated in the last 60 seconds") # a less than b


if __name__ == '__main__':

	unittest.main()



