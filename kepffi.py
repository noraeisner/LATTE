from __future__ import print_function, absolute_import, division

import os
import astropy
import numpy as np
import pandas as pd
import seaborn as sb

import requests
import lightkurve as lk
from os.path import exists
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from astropy.stats import BoxLeastSquares

from glob import glob
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astroquery.mast import Catalogs
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.patches import Rectangle
from lightkurve import TessTargetPixelFile

from tess_stars2px import tess_stars2px_function_entry
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons

tic = '55525572'
indir = "./LATTE_output"
sectors_all = [2,3,4,5]
sectors = [8]
ra = 72.69401
dec = -60.905461
noshow = False


def norm(a):
    '''
    function to normalise the data - used in download_data_tpf
    '''
    m = np.median(a)
    x = a-m
    s = np.std(x)
    x = x/s
    return x, m, s

def unnorm(x,m,s):
    '''
    function to un-normalise the data - used in download_data_tpf
    '''
    y = x * s
    a = y + m
    return a


def download_data_FFI_interact(indir,sector, sectors_all, tic, save = False):
    '''
    Download the LCs for the target star for all the indicated sectors from the FFIs. This uses the lighkurve package.
    
    Parameters
    ----------
    indir : str
        path to where the files will be saved.
    sector  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    tic : str
        TIC (Tess Input Catalog) ID of the target
    binfac  :  int
        The factor by which the data should be binned. Default = 5 (which is what is shown on PHT)

    Returns
    -------
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    allflux_err  :  list
        normalized flux errors (not binned)
    allline  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    allx1  :  list
        CCD column position of target’s flux-weighted centroid. In x direction
    allx2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In x direction
    ally1  :  list
        CCD column position of target’s flux-weighted centroid. In y direction
    ally2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In y direction
    alltimel2  :  list
        time used for the x and y centroid position plottin
    allfbkg  :  list
        background flux
    start_sec  :  list
        times of the start of the sector
    end_sec  :  list
        times of the end of the sector
    in_sec  :  list
        the sectors for which data was downloaded
    tessmag  :  list
        TESS magnitude of the target star
    teff  :  list
        effective temperature of the tagret star (K)
    srad  :  list
        radius of the target star (solar radii)

    '''

    future_sectors = list(set(sectors_all) - set(sector))

    # dowload the data 

    searchtic = 'TIC' + tic

    
    start_sec = []
    end_sec = []
    in_sec = []
    
    alltime_list = []
    allline = []

    X1_list = []
    X1flux_list = []
    X4_list = []
    apmask_list = []
    arrshape_list = []
    tpf_filt_list = []
    t_list = []
    bkg_list = []
    tpf_list = []

    if len(sector) > 1:
        print ("Warning: Downloading data from multiple FFI sectors may take a while to run.")
    
    # open file which shows the momentum dumps
    momentum_dumps_list = "{}/data/tess_mom_dumps.txt".format(indir)
    mom_df = pd.read_csv(momentum_dumps_list, comment = '#', delimiter = '\t')

    for sec in sector:

        # import the data
        print ("Importing FFI data sector {}...".format(sec), end =" ")
        search_result = lk.search_tesscut(searchtic, sector=sec)
        tpf = search_result.download(cutout_size=15)
        
        # get rid of the 'bad' quality data - use the data flags and only take data where the quality = 0. 
        quality = tpf.quality
        tpf = tpf[quality == 0]


        tpf_list.append(tpf)

        print ("Done.\n")

        # extract the information and perform PCA
        print ("Start PCA analysis...", end =" ")
        
        X1 = tpf.flux
        X1flux_list.append(X1)

        arrshape_list.append(X1.shape)
        
        bkg = X1
        bkg = bkg.mean(axis = 0)

        bkg_list.append(bkg)

        # reshape the array in order to perform PCA on it.
        s = X1.shape
        X1 = X1.reshape(s[0],s[1]*s[2])
        
        lkeep = np.isfinite(X1.sum(axis=1)) * (X1.sum(axis=1)>0)
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
        
        ncomp = 5
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
        
        t=tpf.time[lkeep]
        
        alltime = t

         # ------------------------------------------------

        # add all the information to lists that are returned at the end. 
        tpf_filt_list.append(X4.reshape(tpf.flux[lkeep,:,:].shape))

        in_sec.append(sec)

        alltime_list.append(list(alltime))

        allline.append(list((mom_df.loc[mom_df['sec'] == int(sec)])['time']))
        
        start_sec.append([alltime[0]])
        end_sec.append([alltime[-1]])
    
        X1_list.append(X1) #  not corrected
        X4_list.append(X4) #  PCA corrected
        t_list.append(np.array(alltime))  # add it here because this list isn't flattened but the other one is


    alltime_list = [val for sublist in alltime_list for val in sublist]
    allline = [val for sublist in allline for val in sublist]

    del mom_df

    return alltime_list, allline, start_sec, end_sec, in_sec, X1_list, X1flux_list,  X4_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list


def interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, noshow):

    alltime_list, allline, start_sec, end_sec, in_sec, X1_list, X1flux_list,  X4_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI_interact(indir, sectors, sectors_all, tic, save = False)
    
    def close(event):
        plt.close('all')
    

    def extract_LC(aperture):
        ax[0].cla()
        flux = X4[:,aperture.flatten()].sum(axis=1)
        m = np.nanmedian(flux)
        return t, flux/m
          
    def extract_LC2(aperture):
        flux = X4[:,aperture.flatten()].sum(axis=1)
        m = np.nanmedian(flux)
        return t, flux/m


    allfbkg = []
    allfbkg_t = []

    allflux = []
    allflux_flat = []
    allflux_small = []
    apmask_list = []


    for i, X4 in enumerate(X4_list): 
        tpf = tpf_list[i]
        t = t_list[i]
        X1 = X1flux_list[i]
        sec = sectors[i]

        global mask
        global mask2
        global aperture
        global aperture2

        mask = []
        mask2 = []

        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
        aperture2 = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
                

        def onclick(event):
            global mask
            global mask2
            global aperture
            global aperture2

            events = ((int(event.xdata+0.5), int(event.ydata+0.5)))

            if event.inaxes in [ax[1]]:
                print ("event 1")
                [p.remove() for p in reversed(ax[1].patches)]

                if (len(mask) > 0) and (events in list(mask)): 
                    mask = [x for x in mask if x != events] 
                else:
                    mask.append(events)
            
                sqcol = '#ffffee'
                alpha = 0.5

                for pixel in mask:
                    m = int(pixel[0])
                    n = int(pixel[1])
                    r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
                    ax[1].add_patch(r)

                fig.canvas.draw()

            if event.inaxes in [ax[2]]:
                print ("event 2")
                [p.remove() for p in reversed(ax[2].patches)]
            
                if (len(mask2) > 0) and (events in list(mask2)): 
                    mask2 = [x for x in mask2 if x != events] 
                else:
                    mask2.append(events)
            
                sqcol = '#ffffee'
                alpha = 0.5
            
                for pixel2 in mask2:
                    m = int(pixel2[0])
                    n = int(pixel2[1])
                    r2 = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='#c33c7d', facecolor='#e7b1cb', alpha = 0.6)
                    ax[2].add_patch(r2)
                
                fig.canvas.draw()
                
                # update the extraction aperture

            aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
            aperture2 = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)

            for coord in mask:
                aperture[coord] = True

            for coord2 in mask2:
                aperture2[coord2] = True

            ax[0].plot(extract_LC(aperture)[0], extract_LC(aperture)[1],marker='o',color = '#054950', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#054950')
            ax[0].plot(extract_LC2(aperture2)[0], extract_LC2(aperture2)[1],marker='x',color = '#c94f8a', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#c94f8a')
    
            ax[0].set_xlabel("Time")
            ax[0].set_ylabel("Normalized Flux")
            
            fig.canvas.draw_idle()
            
            plt.draw()

        fig, ax = plt.subplots(1,3, figsize=(10,3), gridspec_kw={'width_ratios': [3, 1, 1]})
        
        plt.tight_layout()
        
        ax[0].set_title("Sector {}".format(sec))
        ax[1].set_axis_off()
        ax[2].set_axis_off()
        
        im = ax[1].imshow(X1.mean(axis = 0))
        im2 = ax[2].imshow(X1.mean(axis = 0))

        
        ax[0].set_xlabel("Time")
        ax[0].set_ylabel("Normalized Flux")
        
        ax[1].set_title("Large Aperture")
        ax[2].set_title("Small Aperture")
        fig.canvas.mpl_connect('button_press_event', onclick)
        
        fig.subplots_adjust(left=0.08, bottom=0.2, top=0.90, wspace = 0.1)

        ax[1].text(1.1,-0.2, "Select the pixels for the large and small apertures by \n clicking on the two images", size=9, ha="center", transform=ax[1].transAxes)
        ax[0].set_title("Sector {}".format(sec))
        
        ebx = plt.axes([0.73, 0.05, 0.13, 0.06])
        exit = Button(ebx, 'Close', color='orange')
        exit.on_clicked(close)
        
        plt.show()
        
        # --------------
        # have now chosen an aperture size - now need to extract the LC for the smaller aperture size
        
        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
        for coord in mask[:-1]:
            aperture[coord] = True

        apmask_list.append(aperture)

        target_mask_small = aperture2
        target_mask = aperture 
    

        mask_plot = tpf.plot(aperture_mask=target_mask.T, mask_color='k')
        plt.savefig('{}/{}/{}_mask.png'.format(indir, tic, tic), format='png')
        plt.close('all')

        mask_plot = tpf.plot(aperture_mask=target_mask_small.T, mask_color='k')
        plt.savefig('{}/{}/{}_mask_small.png'.format(indir, tic, tic), format='png')
        plt.close('all')

        flux = X4[:,target_mask.flatten()].sum(axis=1) 
        flux_small = X4[:,target_mask_small.flatten()].sum(axis=1)
        
        m = np.nanmedian(flux)
        m_small = np.nanmedian(flux_small)
        
        # normalize by dividing by the median value
        flux = flux/m
        flux_small = flux_small/m_small
        
        print ("Done.\n")
    
        # -------- flatten the normal lighcurve --------
        print ("Flatten LC...", end =" ")
    
        l = np.isfinite(flux/m)
    
        fr_inj = flux/m
        alltime  = t
    
        T_dur = 0.01  #may need to change!!
        
        nmed = int(720*3*T_dur)
        nmed = 2*int(nmed/2)+1 # make it an odd number 
        ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
        # first number is three time transit durations, the second quite small (10,20 )
    
        l = np.isfinite(ff)
        
        g = interp1d(alltime[l],ff[l],bounds_error=False,fill_value=np.nan)
        ff = g(alltime)
        fr = fr_inj / ff
    
        plt.figure(figsize=(16,5))
        plt.plot(alltime,fr_inj + 0.95,'.')
        plt.plot(alltime,ff + 0.95,'.')
        plt.plot(alltime,fr,'o', color = 'r')
    
        plt.savefig('{}/{}/{}_fit_test.png'.format(indir, tic, tic), format='png')
        plt.clf()
        plt.close()
    
        # ---------------------------------------------
        print ("Done.\n")
        
        print ("Extract background...", end =" ")
        # -------- extract the backrgound flux -----------
        background_mask = ~tpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
        n_background_pixels = background_mask.sum()
    
        background_lc_per_pixel = tpf.to_lightcurve(aperture_mask=background_mask) / n_background_pixels
        n_target_pixels = target_mask.sum()
        background_estimate_lc = background_lc_per_pixel * n_target_pixels
    
        allfbkg.append(background_estimate_lc.flux) # uncorrected background
        allfbkg_t.append(background_estimate_lc.time) 

        print ("Done.\n")
         # ------------------------------------------------
    
        # add all the information to lists that are returned at the end. y
        allflux_flat.append(list(fr))
        allflux.append(list(flux))
        allflux_small.append(list(flux_small))

    allflux_flat = [val for sublist in allflux_flat for val in sublist]
    allflux_small = [val for sublist in allflux_small for val in sublist]
    allflux = [val for sublist in allflux for val in sublist]
    allfbkg = [val for sublist in allfbkg for val in sublist]
    allfbkg_t = [val for sublist in allfbkg_t for val in sublist]
 
    # we're not back to where the same place as with interact_LATTE_FFI

    return alltime, allflux, allflux_small, allflux_flat, allline, allfbkg,allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list


interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, noshow)


