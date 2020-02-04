import os
import sys
import time 
import datetime
import unittest
import numpy as np
from dateutil import parser


import warnings
warnings.filterwarnings("ignore")


import LATTE.LATTEutils as utils

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(utils.__file__))[0:-14]

indir = "./test_output"    

# -----------------------------
'''
Test whether the pyaneti modeling is working correctly. Pyaneti is not available in the current verison of the code. 

-- NOTE: this is testing whether it is working for single transit events and not multi transit
         -- for multi transit, there is a prior on the data - check these are working. 

'''

# -------------
# test with these input parameters (this data is already downloaded so it's only testing that the plotting part works)
tic = '55525572'
sector = '5'
sectors = [5]
transit_list = [1454.7]
transit_sec = '5'

mstar = 1 # this can be assumed to 1 for the modelling
teff  = 5824
srad  = 1.9325
# -------------


class TestPyaneti(unittest.TestCase):
    '''
    Test the extraction of the information from the TP file (already on the system)
    '''    
    def test_pyaneti(self):


        # get data needed to run the BLS - use other unittest to make sure that te data handling works as expected to get to this point. 
        

        if os.path.exists("{}/pyaneti_LATTE.py".format(syspath)):

            print ("Running Pyaneti modelling - this could take a while so be patient...")
        
            transit_list_model =  ("{}".format(str(np.asarray(transit_list)))[1:-1]) # change the list into a string and get rid of the brackets
            # the code is usually run through the command line so call it using the os.system function.
            
            os.system("python3 {}/pyaneti_LATTE.py {} {} {} {} {} {} {}".format(syspath, tic, indir, syspath, mstar, teff, srad, transit_list_model))
        
        else:
            print ("Pyaneti has not been installed so you can't model anything yet. Contact Nora or Oscar for the LATTE version of the Pyaneti code.")


if __name__ == '__main__':

    unittest.main()



