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

# -----------------------------
'''
Test whether the BLS algorithm works correctly and is able to plot. 

'''

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

# test it with the default arguments (except noshow - we don't want the plots to show up in the test run)
global args
args = Namespace(new_data=False, tic='no',sector='no', targetlist='no', 
    noshow=True, o=False, auto=False,nickname='noname', FFI=False, save=True,north=False,new_path=False,mpi=False)
# --------------

class TestBoxLeastSquareTest(unittest.TestCase):
    '''
    Test the extraction of the information from the TP file (already on the system)
    '''    
    def test_BLS(self):

        # get data needed to run the BLS - use other unittest to make sure that te data handling works as expected to get to this point. 

        alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = utils.download_data(indir,sector, tic, binfac = 5, test = './LATTE/tests/tic55525572_lc.fits')

        # run the BLS
        # the plotting is called from within this 
        bls_stats1, bls_stats2 = utils.data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, args)

        # this function should output a plot tahts hows the detrending, and two BLS plots. 
        # Ensure that these are generated and that the number are what we expect them to be. 

        #these are the paths of the files that should have bee created
        BLS1_path = '{}/55525572/55525572_bls_first.png'.format(indir)
        BLS2_path = '{}/55525572/55525572_bls_second.png'.format(indir)

        # Get file's Last modification time stamp only in terms of seconds since epoch 
        time_created_BLS1= os.path.getmtime(BLS1_path)
        time_created_BLS2= os.path.getmtime(BLS2_path)

        # Convert seconds since epoch to readable timestamp
        t_create_BLS1 = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_BLS1)))
        t_create_BLS2 = parser.parse(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_created_BLS2)))
      
        # the time now (needed to get the time since creation)
        t_now = (datetime.datetime.now())

        # -------
        # time difference in minutes
        time_since_creation_BLS1 =  ((t_now - t_create_BLS1).seconds / 60)
        time_since_creation_BLS2 =  ((t_now - t_create_BLS2).seconds / 60)

        self.assertLess(time_since_creation_BLS1, 1, "No BLS plot generated in the last 60 seconds") # a less than b

        # check that the output numbers make sense
        self.assertAlmostEqual(float(bls_stats1[0]), float(16.910000000000014), places=5)
        self.assertAlmostEqual(float(bls_stats1[1]), float(0.3901880448498858), places=5)

        self.assertAlmostEqual(float(bls_stats2[0]), float(0.51))
        self.assertAlmostEqual(float(bls_stats2[1]), float( 0.305186334480843), places=5)


if __name__ == '__main__':

    unittest.main()



