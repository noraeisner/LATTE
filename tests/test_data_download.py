
'''

This script is designed to test that the data is downloaded correctly.


run with: python -m unittest tests/test_data_download.py

# NOTE: requires intenet connection

'''

import os
import unittest
import requests
import warnings

warnings.filterwarnings("ignore")
from LATTE import LATTEutils

# get the outir (or indir) path
syspath = str(os.path.abspath(LATTEutils.__file__))[0:-14]

with open("{}/_config.txt".format(syspath), 'r') as f:
	outdir = str(f.readlines()[-1])

class TestTESSpoint(unittest.TestCase):

	'''
	Test that TESS point returns the right RA and Dec

	Note: requires internet connection
	'''
	def test_testpoint(self):
		output = LATTEutils.tess_point(outdir, '55525572')

		self.assertEqual(output[1], 72.69404300000001, "TESS point RA is not correct")
		self.assertEqual(output[2], -60.905449, "TESS point Dec is not correct")


class TestNearestNeighbours(unittest.TestCase):

	'''
	test that the code can calculate the nearest neighbour information correctly - involves calculations so ensures that they are correct and work as expected.
	'''
	def test_nn(self):
		output = LATTEutils.nn_ticids(outdir, [5], '55525572')

		self.assertEqual(output[0], [55525572, 55525518, 358806415, 55557836, 55557135, 55522986], "nearest neighbour TIC IDs are incorrect")
		self.assertEqual(output[1], [0.0, 6.795727328666717, 10.622509290130298, 15.63959237521102, 16.881930554329756, 17.79841005439404], "nearest neighbour distances are incorrect")
		self.assertEqual(output[2], 72.6941, "RA is incorrect")
		self.assertEqual(output[3], -60.9055, "DEC is incorrect")


class TestDownloadLC(unittest.TestCase):
	'''
	test to ensure that the data download works as it should 
	- test that it return the right file withthe right values.

	NOTE: will only work if the internet connection works. 

	'''
	def test_downloadLC(self):
		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = LATTEutils.download_data(outdir, [5], '55525572', binfac = 5)


		self.assertEqual(float(alltime[0]), 1437.9785214992496, "alltime is incorrect")
		self.assertEqual(float(allflux[1000]), 1.0011740922927856, "allflux is incorrect")
		self.assertEqual(float(allflux_err[1000]), 0.0009966295910999179, "allflux_err is incorrect")
		self.assertEqual(float(all_md[0]), 1441.0243360822994, "all_md is incorrect")
		self.assertEqual(float(alltimebinned[0]), 1437.9812992567897, "alltimebinned is incorrect")
		self.assertEqual(float(allfluxbinned[1000]), 1.000087034702301, "alltimebinned is incorrect")
		
		self.assertEqual(float(allx1[0]), -0.006673532877357502, "allx1 is incorrect")
		self.assertEqual(float(allx2[0]), -0.008006655611097813, "allx2 is incorrect")
		self.assertEqual(float(ally1[0]), -0.018879381330179967, "ally1 is incorrect")
		self.assertEqual(float(ally2[0]), -0.02256738394498825, "ally2 is incorrect")
		self.assertEqual(float(alltimel2[0]), 1437.9924102871835, "alltimel2 is incorrect")
		self.assertEqual(float(allfbkg[1000]), 974.2296752929688, "allfbkg is incorrect")
		
		self.assertEqual(float(start_sec[0][0]), 1437.9785214992496, "start_sec is incorrect")
		self.assertEqual(float(end_sec[0][0]), 1464.2879709506096, "end_sec is incorrect")
		self.assertEqual(float(in_sec[0]), 5.0, "in_sec is incorrect")

		self.assertEqual(float(tessmag), 9.81999969, "tessmag is incorrect")
		self.assertEqual(float(teff), 5823.70019531, "teff is incorrect")
		self.assertEqual(float(srad), 1.93253005, "srad is incorrect")


if __name__ == '__main__':

	unittest.main()

