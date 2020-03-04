
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

		self.assertAlmostEqual(output[1], 72.69404300000001,places=5,  msg ="TESS point RA is not correct")
		self.assertAlmostEqual(output[2], -60.905449, places=4, msg ="TESS point Dec is not correct")


class TestNearestNeighbours(unittest.TestCase):

	'''
	test that the code can calculate the nearest neighbour information correctly - involves calculations so ensures that they are correct and work as expected.
	'''
	def test_nn(self):
		output = LATTEutils.nn_ticids(outdir, [5], '55525572')

		self.assertEqual(output[0], [55525572, 55525518, 358806415, 55557836, 55557135, 55522986], msg="nearest neighbour TIC IDs are incorrect")

		self.assertEqual(output[1][0], 0.0, msg="nearest neighbour distances are incorrect")
		self.assertAlmostEqual(output[1][1], 6.795727328666717, places=5, msg="nearest neighbour distances are incorrect")
						
		self.assertAlmostEqual(output[1][2], 10.622509290130298,places=5, msg="nearest neighbour distances are incorrect")
		self.assertAlmostEqual(output[1][3], 15.63959237521102, places=5, msg="nearest neighbour distances are incorrect")
						

		self.assertAlmostEqual(output[2], 72.6941, places=3, msg="RA is incorrect")
		self.assertAlmostEqual(output[3], -60.9055, places=3, msg="DEC is incorrect")


class TestDownloadLC(unittest.TestCase):
	'''
	test to ensure that the data download works as it should 
	- test that it return the right file withthe right values.

	NOTE: will only work if the internet connection works. 

	'''
	def test_downloadLC(self):
		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = LATTEutils.download_data(outdir, [5], '55525572', binfac = 5)


		self.assertAlmostEqual(float(alltime[0]), 1437.9785214992496, places=5, msg="alltime is incorrect")
		self.assertAlmostEqual(float(allflux[1000]), 1.0011740922927856, places=5, msg="allflux is incorrect")
		self.assertAlmostEqual(float(allflux_err[1000]), 0.0009966295910999179, places=5,msg= "allflux_err is incorrect")
		self.assertAlmostEqual(float(all_md[0]), 1441.0243360822994, places=5, msg="all_md is incorrect")
		self.assertAlmostEqual(float(alltimebinned[0]), 1437.9812992567897, places=5, msg="alltimebinned is incorrect")
		self.assertAlmostEqual(float(allfluxbinned[1000]), 1.000087034702301, places=5, msg="alltimebinned is incorrect")
		
		self.assertAlmostEqual(float(allx1[0]), -0.006673532877357502,places=5,  msg="allx1 is incorrect")
		self.assertAlmostEqual(float(allx2[0]), -0.008006655611097813, places=5, msg="allx2 is incorrect")
		self.assertAlmostEqual(float(ally1[0]), -0.018879381330179967, places=5,msg= "ally1 is incorrect")
		self.assertAlmostEqual(float(ally2[0]), -0.02256738394498825, places=5, msg="ally2 is incorrect")
		self.assertAlmostEqual(float(alltimel2[0]), 1437.9924102871835, places=5, msg="alltimel2 is incorrect")
		self.assertAlmostEqual(float(allfbkg[1000]), 974.2296752929688, places=5, msg="allfbkg is incorrect")
		
		self.assertAlmostEqual(float(start_sec[0][0]), 1437.9785214992496, places=5, msg="start_sec is incorrect")
		self.assertAlmostEqual(float(end_sec[0][0]), 1464.2879709506096, places=5, msg="end_sec is incorrect")
		self.assertAlmostEqual(float(in_sec[0]), 5.0, msg="in_sec is incorrect")

		self.assertAlmostEqual(float(tessmag), 9.81999969, places=5, msg="tessmag is incorrect")
		self.assertAlmostEqual(float(teff), 5823.70019531, places=5, msg="teff is incorrect")
		self.assertAlmostEqual(float(srad), 1.93253005, places=5, msg="srad is incorrect")


if __name__ == '__main__':

	unittest.main()

