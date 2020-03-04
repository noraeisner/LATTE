
# -----------------------------
'''
This code makes sure that everytime that the returned data is what it is expected to be. 

This test case is carried out with a pre-downloaded lighcurve file from TIC 55525572 Sector 5.


run using: python -m unittest tests/test_data_products.py

'''

# -------------


import os
import mock
import unittest
import requests
import argparse
from argparse import ArgumentParser

import warnings
warnings.filterwarnings("ignore")

import sys
import LATTE.LATTEutils as utils

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(utils.__file__))[0:-14]


indir = "./LATTE/output"        

# test with these input parameters (this data is already downloaded so it's only testing that the data extraction from the FITS files works)
tic = '55525572'
sector = '5'
sectors = [5]
transit_list = [1454.7]
transit_sec = '5'
# -------------

class TestDataImport_LC(unittest.TestCase):
    '''
    Test the extraction of the information from the LC file (which is already on the system)
    '''

    def test_LC_request_response(self):
        # Call the service, which will send a request to the server.
        alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = utils.download_data(indir,sector, tic, binfac = 5, test = './LATTE/tests/tic55525572_lc.fits')

        # If the request is sent successfully, then I expect a response to be returned.
        self.assertAlmostEqual(float(alltime[100]),       float(1438.1174094083246), places=5, msg='LC data is not what it should be.')
        self.assertAlmostEqual(float(allflux[100]),       float(1.0003796815872192), places=5)
        self.assertAlmostEqual(float(allflux_err[100]),   float(0.0011117355898022652), places=5)
        self.assertAlmostEqual(float(all_md[0]),          float(1441.0243360822994), places=5)
        self.assertAlmostEqual(float(alltimebinned[100]), float(1438.6757392292352), places=5)
        self.assertAlmostEqual(float(allfluxbinned[100]), float(0.9995232224464417), places=5)
        self.assertAlmostEqual(float(allx1[100]),         float(-0.010814569022272735), places=5)
        self.assertAlmostEqual(float(allx2[100]),         float(-0.011804798617959023), places=5)
        self.assertAlmostEqual(float(ally1[100]),         float(-0.024266568269581512), places=5)
        self.assertAlmostEqual(float(ally2[100]),         float(-0.02981671877205372), places=5)
        self.assertAlmostEqual(float(alltime12[100]),     float(1438.1312982026025), places=5)


class TestDataImport_TP(unittest.TestCase):
    '''
    Test the extraction of the information from the TP file (which is already on the system)
    '''

    def test_LC_request_response(self):

        # Call the service, which will send a request to the server.
        X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = utils.download_tpf_mast(indir, transit_sec, transit_list, tic, test = './LATTE/tests/tic55525572_tp.fits')

        # If the request is sent successfully, then I expect a response to be returned.
        self.assertAlmostEqual(float(X1_list[0][0][0]),float(23.402481079101562), places=5, msg ='TP data is not what it should be.')
        self.assertAlmostEqual(float(oot_list[0][0]),float(0.0))
        self.assertAlmostEqual(float(intr_list[0][0]),float(0.0))
        self.assertAlmostEqual(float(bkg_list[0][0][0]),float(29.239688873291016), places=5)
        self.assertAlmostEqual(float(apmask_list[0][0][0]),float(0.0))
        self.assertAlmostEqual(float(arrshape_list[0][0]),float(18944.0))
        self.assertAlmostEqual(float(t_list[0][0]),float(1437.9924102871835), places=5,)
        self.assertAlmostEqual(float(T0_list[0]),float(1454.7))


class TestDataImport_TP_lighkurve(unittest.TestCase):
    '''
    Test the extraction of the information from the TP file (which is already on the system)
    '''

    def test_LC_request_response(self):

        # Call the service, which will send a request to the server.
        TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list = utils.download_tpf_lightkurve(indir, transit_list, sectors, tic, test = './LATTE/tests/tic55525572_tp.fits')
        
        print(float(TESS_unbinned_t_l[0]))
        print(float(TESS_binned_t_l[0]))
        print(float(small_binned_t_l[0]))
        print(float(TESS_unbinned_l[0]))
        print(float(TESS_binned_l[0]))
        print(float(small_binned_l[0]))

        
        # If the request is sent successfully, then I expect a response to be returned.
        self.assertAlmostEqual(float(TESS_unbinned_t_l[0]), float(1437.9924102871835), msg='TP lighkurve data is not what it should be.')
        self.assertAlmostEqual(float(TESS_binned_t_l[0]), float(1437.9972713631325), places=5)
        self.assertAlmostEqual(float(small_binned_t_l[0]), float(1437.9972713631325), places=5)
        self.assertAlmostEqual(float(TESS_unbinned_l[0]), float(0.9975904130977179), places=5)
        self.assertAlmostEqual(float(TESS_binned_l[0]), float(0.9986087992155438), places=5)
        self.assertAlmostEqual(float(small_binned_l[0]), float(0.9971978343909532), places=5)



# ---------------------
if __name__ == "__main__":

    unittest.main()

