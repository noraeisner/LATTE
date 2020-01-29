
'''

This script is designed to test the LATTE python code to ensure that all the major function are working as they should. 

'''
import os
import mock
import unittest
import requests
import numpy as np
import seaborn as sb
import astropy.io.fits as pf
from argparse import ArgumentParser

from sklearn.decomposition import PCA
from scipy.interpolate import interp1d

import warnings
warnings.filterwarnings("ignore")

import sys
import LATTE.LATTEutils as utils

# test the downloading of the data
# get the outir (or indir) path
syspath = str(os.path.abspath(utils.__file__))[0:-14]

#with open("{}/_config.txt".format(syspath), 'r') as f:
#    indir = str(f.readlines()[-1])


indir = "./output"        
# -----------------------------
'''
This code makes sure that everytime that the data is called to plot something, the data is still what I expect it to be. 

This test case is carried out with a pre-downloaded lighcurve file from TIC 55525572 Sector 5.


Run each plotting function and make sure that it produces the output file that I would expect it to produce. 
 - don't check the actual file as this is computer dependednt and therefore not a robust test. 
'''

# ------------------------------

tic = '55525572'
sector = '5'
transit_list = [1454.7]
transit_sec = '5'

def download_data(indir,sector, tic, binfac = 5):


    # plot using the seaborn library
    sb.set(style='ticks')
    sb.set_color_codes('deep')
    
    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
    

    # file on local computer - this test only works on my computer...
    dwload_link = ['./tests/tic55525572_lc.fits']    

    # define all the empty lists to append to in order to return the data that will be requrides later on in the script
    alltimebinned = []
    allfluxbinned = []
    
    allx1 = []
    allx2 = []
    ally1 = []
    ally2 = []
    alltimel2 = []
    allfbkg = []
    
    start_sec = []
    end_sec = []
    in_sec = []
    
    alltime = []
    allflux = []
    allflux_err = []
    all_md = []
    
    # loop through all the download links - all the data that we want to access
    for lcfile in dwload_link:
        
        # open the file using the response url  
        lchdu  = pf.open(lcfile) # this needs to be a URL - not a filefile
        
        # open and view columns in lightcurve extension
        lcdata = lchdu[1].data
        lchdu[1].columns

        f02 = lcdata['PDCSAP_FLUX'] # Presearch Data Conditioning 
        f02_err = lcdata['PDCSAP_FLUX_ERR']
        quality = lcdata['QUALITY']  # quality flags as determined by the SPOC pipeline 
        time    = lcdata['TIME']
        f0     = lcdata['SAP_FLUX']  # 
        fbkg     = lcdata['SAP_BKG']  # background flux 
        
        med = np.nanmedian(f02)  # determine the median flux (ignore nan values)
        f1 = f02/med  # normalize by dividing byt the median flux
        f1_err = f02_err/med  # normalise the errors on the flux
        
        x1      = lcdata['MOM_CENTR1']  # CCD column position of target’s flux-weighted centroid 
        x1      -= np.nanmedian(x1)
        y1      = lcdata['MOM_CENTR2']  
        y1      -= np.nanmedian(y1)
        x2      = lcdata['POS_CORR1'] # The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects.
        x2      -= np.nanmedian(x2)
        y2      = lcdata['POS_CORR2']
        y2      -= np.nanmedian(y2)
        l       = (quality>0)   # good quality data
        l2      = (quality<=0)  # bad quality data
        
        sec     = int(lchdu[0].header['SECTOR'])  # the TESS observational sector

        tessmag = lchdu[0].header['TESSMAG']  # magnitude in the FITS header
        teff    = lchdu[0].header['TEFF']     # effective temperature in the FITS header (kelvin)
        srad    = lchdu[0].header['RADIUS']   # stellar radius in the FITS header (solar radii) 

        flux     = lcdata['SAP_FLUX']

        # store the sector we are looking at
        in_sec.append(sec)

        # binned data
        N       = len(time)
        n       = int(np.floor(N/binfac)*binfac)
        X       = np.zeros((2,n))
        X[0,:]  = time[:n]
        X[1,:]  = f1[:n]
        Xb      = rebin(X, (2,int(n/binfac)))

        time_binned    = Xb[0]
        flux_binned    = Xb[1]

        # the time of the momentum dumps are indicated by the quality flag
        mom_dump = np.bitwise_and(quality, 2**5) >= 1

        # store the relevant information in the given list
        alltime.append(list(time)) 
        allflux.append(list(f1)) 
        allflux_err.append(list(f1_err))
        all_md.append(list(time[mom_dump]))
        
        alltimebinned.append(list(time_binned))
        allfluxbinned.append(list(flux_binned))
        
        allx1.append(list(x1[l2]))
        allx2.append(list(x2[l2]))
        ally1.append(list(y1[l2]))
        ally2.append(list(y2[l2]))
        alltimel2.append(list(time[l2]))
        
        allfbkg.append(fbkg)
    
        start_sec.append([time[0]])
        end_sec.append([time[-1]])
    
    # flatten lists of lists
    alltime = [val for sublist in alltime for val in sublist]
    allflux = [val for sublist in allflux for val in sublist]
    allflux_err = [val for sublist in allflux_err for val in sublist]
    all_md = [val for sublist in all_md for val in sublist]

    alltimebinned = [val for sublist in alltimebinned for val in sublist]
    allfluxbinned = [val for sublist in allfluxbinned for val in sublist]
    
    allx1 = [val for sublist in allx1 for val in sublist]
    allx2 = [val for sublist in allx2 for val in sublist]
    ally1 = [val for sublist in ally1 for val in sublist]
    ally2 = [val for sublist in ally2 for val in sublist]
    alltimel2 = [val for sublist in alltimel2 for val in sublist]
    allfbkg = [val for sublist in allfbkg for val in sublist]
   

    return alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad

def download_tpf_mast(indir, transit_sec, transit_list, tic):
    '''
    Download the TPF LCs for the target star for all the indicated sectors. Not using Lightkurve
    
    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_sec  :  list or str
        list of the sectors that have a transit in them. If 'all', all the sectors in whic the target appears will be downloaded
    tic : str
        TIC (Tess Input Catalog) ID of the target
    transit_list  :  int
        list of all the marked transits

    Returns
    -------
    X1_list  :  list
        flux vs time for each pixel
    X4_list  :  list
        PCA corrected flux vs time for each pixel
    oot_list  :  list
        out of transit mask
    intr_list  :  list
        in transit mask
    bkg_list  :  list
        the flux that was used to normalise each pixel - i.e. what is used to make the background plot colour for each pixel.
    apmask_list  :  list
        aperture masks from the pipeline
    arrshape_list  :  list
        shape of the array
    t_list  :  list
        time arrays
    T0_list  :  list
        list of the peaks
    tpf_filt_list   : list 
        list of the filteres tpfs
    '''

    X1_list = []
    X4_list = []
    oot_list = []
    intr_list = []
    bkg_list = []
    apmask_list = []
    arrshape_list = []
    t_list = []
    T0_list = []
    tpf_filt_list = []


    dwload_link_tp = ['./tests/tic55525572_tp.fits']


    for file in dwload_link_tp:

        tpf = pf.open(file)   # open the file

        for T0 in transit_list:

            X1 = tpf[1].data['FLUX']
            arrshape = X1.shape
            
            Q = tpf[1].data['QUALITY']
            t = tpf[1].data['TIME']
            
            aperture = tpf[2].data
            apmask = aperture >= np.nanmax(aperture)


            if (T0 > np.nanmin(t)) and (T0 < np.nanmax(t)):

                lkeep = (Q==0)
                bkg = X1[lkeep,:]
                bkg = bkg.mean(axis = 0)
                    
                bkg_list.append(bkg)
                apmask_list.append(apmask)
                arrshape_list.append(arrshape)

                s = X1.shape
                X1 = X1.reshape(s[0],s[1]*s[2])
                
                lkeep = np.isfinite(X1.sum(axis=1)) * (Q==0) * (X1.sum(axis=1)>0)
                X1 = X1[lkeep,:]
                
                X2 = np.zeros_like(X1)
                M = np.zeros(X1.shape[1])
                S = np.zeros(X1.shape[1])
                for n in range(X1.shape[1]):
                    a = X1[:,n]
                    x, m, s = norm(a)
                    X2[:,n]=x
                    M[n] = m
                    S[n] = s
                
                ncomp = 4
                pca = PCA(n_components=ncomp)
                trends = pca.fit_transform(X2)
                weights = pca.components_
                
                X3 = np.copy(X2)
                for n in range(X2.shape[1]):
                    for m in range(ncomp):
                        X3[:,n] -= trends[:,m] * weights[m,n]
                        

                X4 = np.zeros_like(X3)
                for n in range(X2.shape[1]):
                    x = X3[:,n]
                    a = unnorm(x, M[n], S[n])
                    X4[:,n] = a
                
                t=tpf[1].data['TIME'][lkeep]
                
                tpf_filt = X4.reshape(tpf[1].data['FLUX'][lkeep,:,:].shape)


                oot = (abs(T0-t) < 0.56) * (abs(T0-t) < 0.3) 
                intr = abs(T0-t) < 0.1
        
                X1_list.append(X1) # not corrected
                X4_list.append(X4) #  PDC corrected
                oot_list.append(oot)  # out of transit filter
                intr_list.append(intr)  # in transit filter
                t_list.append(t)
                T0_list.append(T0)
                tpf_filt_list.append(tpf_filt)


    return X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list

class TestDataImport_LC(unittest.TestCase):
    
    def test_LC_request_response(self):
        # Call the service, which will send a request to the server.
        alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir,sector, tic, binfac = 5)

        # If the request is sent successfully, then I expect a response to be returned.
        self.assertEqual(float(alltime[100]),       float(1438.1174094083246))
        self.assertEqual(float(allflux[100]),       float(1.0003796815872192))
        self.assertEqual(float(allflux_err[100]),   float(0.0011117355898022652))
        self.assertEqual(float(all_md[0]),          float(1441.0243360822994))
        self.assertEqual(float(alltimebinned[100]), float(1438.6757392292352))
        self.assertEqual(float(allfluxbinned[100]), float(0.9995232224464417))
        self.assertEqual(float(allx1[100]),         float(-0.010814569022272735))
        self.assertEqual(float(allx2[100]),         float(-0.011804798617959023))
        self.assertEqual(float(ally1[100]),         float(-0.024266568269581512))
        self.assertEqual(float(ally2[100]),         float(-0.02981671877205372))
        self.assertEqual(float(alltimel2[100]),     float(1438.1312982026025))

def norm(a):
    '''
    function to normalise the data - used in download_tpf_mast
    '''
    m = np.median(a)
    x = a-m
    s = np.std(x)
    x = x/s
    return x, m, s

def unnorm(x,m,s):
    '''
    function to un-normalise the data - used in download_tpf_mast
    '''
    y = x * s
    a = y + m
    return a


class TestDataImport_TP(unittest.TestCase):
    
    def test_LC_request_response(self):
        # Call the service, which will send a request to the server.
        X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = download_tpf_mast(indir, transit_sec, transit_list, tic)

        # If the request is sent successfully, then I expect a response to be returned.
        self.assertEqual(float(X1_list[0][0][0]),float(23.402481079101562))
        self.assertEqual(float(oot_list[0][0]),float(0.0))
        self.assertEqual(float(intr_list[0][0]),float(0.0))
        self.assertEqual(float(bkg_list[0][0][0]),float(29.239688873291016))
        self.assertEqual(float(apmask_list[0][0][0]),float(0.0))
        self.assertEqual(float(arrshape_list[0][0]),float(18944.0))
        self.assertEqual(float(t_list[0][0]),float(1437.9924102871835))
        self.assertEqual(float(T0_list[0]),float(1454.7))


# ----------------------
# not that the data has been iported, test that it can make the plots that you want it to make
# ----------------------

def create_parser():
    ap = ArgumentParser(description='Lightcurve Analysis Tool for Transiting Exoplanets')
    ap.add_argument('--new-data', action='store_true')
    ap.add_argument('--tic', type=str, help='the target TIC id, if blank will be prompted to enter', default = 'no')
    ap.add_argument('--sector', type=str, help='the sector(s) to look at', default = 'no')
    ap.add_argument('--targetlist', type=str, help='the link to the target file list', default='no')
    ap.add_argument('--noshow', action='store_true', help='if you want to NOT show the plots write --noshow in the command line')
    ap.add_argument('--o', action='store_true', help='if you call this old files will be overwriten in the non-interactive version')
    ap.add_argument('--auto', action='store_true', help='automatic aperture selection')
    ap.add_argument('--nickname', type=str, help='give the target a memorable name', default='noname')
    ap.add_argument('--FFI', action='store_true', help='is this an FIIs?')
    ap.add_argument('--save', help='is this an FIIs?', default = True)
    ap.add_argument('--north', action='store_true', help='write "north" if you want all the images to be aligned North')
    ap.add_argument('--new-path', action = 'store_true', help ='enter to change the output path')

    args = ap.parse_args()

    args.noshow = True

    return args



def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    args = parser.parse_args()
    print(args)  # NOTE: this is how you would check what the kwargs are if you're unsure
    return args.accumulate(args.integers)


@mock.patch('argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(accumulate=sum, integers=[1,2,3]))
def test_command(mock_args):
    res = main()
    assert res == 6, "1 + 2 + 3 = 6"


if __name__ == "__main__":
    print(main())



    
#class TestDataPlottig(unittest.TestCase):
#    
#    def setUp(self):
#        self.args = create_parser()
#
#    def test_plot(self):
#        parser = self.parser.parse_args(['--noshow'])
#        
#        alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir,sector, tic, binfac = 5)
#        X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = download_tpf_mast(indir, transit_sec, transit_list, tic)
#        
#        #utils.plot_full_md(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, transit_list, args)
#        

#class TestDataPlottig(unittest.TestCase):
#
#
#    def test_plot(self):



        ## ---------
        #utils.plot_full_md(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, transit_list, args)
        #
        #utils.plot_centroid(tic, indir, alltime12, allx1, ally1, allx2, ally2, transit_list, args)
        #
        ##utils.plot_aperturesize(tic, indir, TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, transit_list, args)
        #
        #utils.plot_background(tic, indir, alltime, allfbkg, transit_list, args)
    
        #utils.plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, transit_list, args)


if __name__ == '__main__':


    unittest.main(args)











