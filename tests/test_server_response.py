
'''

This script is designed to test the LATTE python code to ensure that the server response is working. 

NOTE: requires connection to the internet and to the archive database. If these tests fail, first check that these connections 
are working by accessing the website: https://archive.stsci.edu/missions/tess/download_scripts/ in a browser.

'''
import os
import unittest
import requests

import warnings
warnings.filterwarnings("ignore")

import sys
import LATTE.LATTEutils as utils

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(utils.__file__))[0:-14]

with open("{}/_config.txt".format(syspath), 'r') as f:
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


if __name__ == '__main__':

	unittest.main()

