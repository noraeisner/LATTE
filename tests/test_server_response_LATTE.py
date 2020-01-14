
'''

This script is designed to test the LATTE python code to ensure that all the major function are working as they should. 

'''

import unittest
import requests

import warnings
warnings.filterwarnings("ignore")

import sys

LATTEutils = __import__("/Users/Nora/Documents/research/TESS/planethunters/code/LATTE/LATTE/LATTEutils")

# test the downloading of the data
indir = "/Users/Nora/Documents/research/TESS/planethunters/code/LATTE/_config.txt"

with open("_config.txt", 'r') as f:
	indir = str(f.readlines()[-1])
				

def download_LC_data():
	# function to reach the external server to download the scipt - this is for the first sector only
	LC_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_1_lc.sh"
	r_LC = requests.get(LC_url) # create HTTP response object
	
	return r_LC.status_code

def download_TP_data():
	# function to reach the external server to download the scipt - this is for the first sector only
	TP_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_1_tp.sh"
	r_TP = requests.get(TP_url) # create HTTP response object 

	return r_TP.status_code	

def download_target_list():
	# function to reach the external server to download the scipt - this is for the first sector only
	target_list = "https://tess.mit.edu/wp-content/uploads/all_targets_S001_v1.txt"
	r_target_list = requests.get(target_list) # create HTTP response object

	return r_target_list.status_code	


class TestServerResponse(unittest.TestCase):
	
	def test_LC_request_response(self):
		# Call the service, which will send a request to the server.
		responseLC = download_LC_data()
		responseTP = download_TP_data()
		responseTL = download_target_list()

		# If the request is sent successfully, then I expect a response to be returned.
		self.assertEqual(responseLC, 200, "LC data Download link does not work - can't connect to the server")
		self.assertEqual(responseTP, 200, "TP data Download link does not work - can't connect to the server")
		self.assertEqual(responseTL, 200, "Target list data Download link does not work - can't connect to the server")


class TestTESSpoint(unittest.TestCase):

	def test_testpoint(self):
		output = tess_point(indir, '55525572')

		self.assertEqual(output[0], [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13], "TESS sectors are not correct")
		self.assertEqual(output[1], 72.69404300000001, "TESS point RA is not correct")
		self.assertEqual(output[2], -60.905449, "TESS point Dec is not correct")


class TestNearestNeighbours(unittest.TestCase):

	def test_nn(self):
		output = nn_ticids(indir, [5], '55525572')

		self.assertEqual(output[0], [55525572, 55525518, 358806415, 55557836, 55557135, 55522986], "nearest neighbour TIC IDs are incorrect")
		self.assertEqual(output[1], [0.0, 6.795727328666717, 10.622509290130298, 15.63959237521102, 16.881930554329756, 17.79841005439404], "nearest neighbour distances are incorrect")
		self.assertEqual(output[2], 72.6941, "RA is incorrect")
		self.assertEqual(output[3], -60.9055, "DEC is incorrect")


class TestDownloadLC(unittest.TestCase):

	def test_downloadLC(self):
		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir, [5], '55525572', binfac = 5)


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

