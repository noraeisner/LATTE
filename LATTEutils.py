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
from astroplan import FixedTarget
from astropy.stats import BoxLeastSquares
from astroplan.plots import plot_finder_image
from astropy.stats import median_absolute_deviation

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

# custom modules
import filters
import LATTEbrew as brew

# ------ interact -------
def interact_LATTE(tic, indir, sectors_all, sectors, ra, dec, noshow):
    '''
    Function to run the Interactive LATTE code using the matplotlib interactive tool.
    Calls the plot where the transit-event times can be identifies and the plotting/modeling options specified.
    
    Parameters
    ----------
    tic  :   str
        target TIC ID
    indir  :  str
        path to directory where all the plots and data will be saved. 
    sectors_all  :   list
        all the sectors in which the target has been/ will be observed
    sectors  :  list
        the sectors which will be analysed

    Returns
    -------
        runs the brew_LATTE code...
    '''

    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
    

    print ("Start data download.....", end =" ")
    alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir, sectors, tic)
    print ("Done.\n")
    
    plt.close('all')
    # -------------------------
    # Plot the interactive plot
    # -------------------------
    
    fig, ax = plt.subplots(2, 1, figsize=(10,7))
    plt.tight_layout()
    
    # Adjust tbplots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.25)
    
    fluxmin = np.nanmin(allflux)
    fluxmax = np.nanmax(allflux)
    
    # Draw the initial plot
    # The 'line' variable is used for modifying the line later

    def cutout(transit):
        mask_binned = (np.array(alltimebinned) > transit-1) & (np.array(alltimebinned) < transit+1)
        mask = (np.array(alltime) > transit-1) & (np.array(alltime) < transit+1)
        
        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux), np.array(alltimebinned)[mask_binned], np.array(allfluxbinned)[mask_binned], np.array(alltimebinned), np.array(allfluxbinned)]
    
    def binning(binfac):
        # binned data
        N      = len(alltime)
        n      = int(np.floor(N/binfac)*binfac)
        X      = np.zeros((2,n))
        X[0,:]  = alltime[:n]
        X[1,:]  = allflux[:n]
        Xb    = rebin(X, (2,int(n/binfac)))
        
        time_binned = Xb[0]
        flux_binned = Xb[1]
    
        return [time_binned, flux_binned]
    
    transit = np.nanmean(alltimebinned)
    binfac = 7
    
    [line_full] = ax[0].plot(alltime, allflux , marker='o',lw = 0, markersize = 4, color = 'orange', alpha = 0.8, label = 'unbinned', markerfacecolor='white')
    [line_full_binned] = ax[0].plot(binning(binfac)[0], binning(binfac)[1],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, label = 'binning = 7', markerfacecolor='k')
    
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',lw = 0, markersize = 4, color = 'orange', alpha = 0.8, label = 'unbinned', markerfacecolor='white')
    [line_binned] =  ax[1].plot(cutout(transit)[4], cutout(transit)[5],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, label = 'binning = 7', markerfacecolor='k')
    
    # Define an axes area and draw a slider in it
    transit_slider_ax  = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=transit, color='teal')
    
    scale_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
    scale_slider = Slider(scale_slider_ax, 'Y-Axis Scale', 0.99, 1.01, valinit=1, color='silver')

    
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])
    
    ax[1].set_xlim([np.nanmean(alltime)-1, np.nanmean(alltime)+1])
    ax[1].set_ylim([fluxmin, fluxmax])
    
    # Define an action for modifying the line when any slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])
    
        line_binned.set_xdata(cutout(transit_slider.val)[4])
        line_binned.set_ydata(cutout(transit_slider.val)[5])
    
        fig.canvas.draw_idle()
    
    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)
    
    def update_axis(val):   
        ax[1].set_xlim([transit_slider.val - 1,transit_slider.val + 1])
        
        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)
    
    def update_yaxis(val):  
    
        med = 1
        diff = abs(med - (fluxmin * scale_slider.val))
    
        ax[0].set_ylim([med - diff ,med + diff])
        ax[1].set_ylim([med - diff ,med + diff])
    
    
    transit_slider.on_changed(update_axis)
    scale_slider.on_changed(update_yaxis)
    transit_slider.on_changed(sliders_on_changed)

    # Determine whether to save the values the plots or not left, bottom, width, height

    var_ax = fig.add_axes([0.025, 0.3, 0.1, 0.15])
    save_var = CheckButtons(var_ax, ('Simple', 'BLS', 'model', 'Save', 'DVR'), (False, False, False, True, False))
    
    simple = False
    BLS = False
    model = False
    save = True
    DV = False

    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        BLS = status[1]
        model = status[2]
        save = status[3]
        DV = status[4]

    # Add a set of radio buttons for changing color. slider = [left, bottom, width, height]
    binning_ax = fig.add_axes([0.025, 0.5, 0.10, 0.15])
    binning_radios = RadioButtons(binning_ax, ('2', '5', '7', '10'), active=0)
    
    def binning_button(label):
        line_full_binned.set_xdata(binning(int(label))[0])
        line_full_binned.set_ydata(binning(int(label))[1])
        fig.canvas.draw_idle()
    
    binning_radios.on_clicked(binning_button)
    
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf
    
    ax[0].tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = 12)
    ax[0].vlines(allline, minf-1,minf + height*0.3 , colors = 'r', label = "Momentum Dump")
    
    ax[1].tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("BJD-2457000", fontsize = 12)
    ax[1].set_ylabel("Normalised Flux", fontsize = 12)
    ax[1].vlines(allline, minf-0.5,minf-0.5 + height*0.3 , colors = 'r', label = "Momentum Dump")
    
    #plt.text(4.5, -3.15, "(Press Enter, then close window)", fontsize=10, verticalalignment='center')
    plt.text(0.05, 1.1, "Binning Factor", fontsize=10, verticalalignment='center')
    
    initial_text = ""
    
    transit_times = []

    def submit(text):
        ydata = eval(text)
        transit_times.append(ydata)
    
    axbox = plt.axes([0.25, 0.04, 0.50, 0.04])
    text_box = TextBox(axbox, 'Enter transit-event times', initial=initial_text)
    text_box.on_submit(submit)
    

    ebx = plt.axes([0.77, 0.04, 0.13, 0.04])
    exit = Button(ebx, 'Close', color='orange')

    def close(event):
        plt.close('all')

    exit.on_clicked(close)

    plt.show()
    
    end_status = save_var.get_status()

    simple = end_status[0]
    BLS = end_status[1]
    model = end_status[2]
    save = end_status[3]
    DV = end_status[4]


    if len(transit_times) == 0:
        print ("\n Warning: You can't continue without entering a transit-time.\n")
        print ("Exit.\n")
        raise SystemExit
        #exit[0]
    
    if type(transit_times[-1]) == tuple:
        peak_list = list(transit_times[-1])
        peak_list = [float(i) for i in peak_list]
    
    else:
        peak_list = [float(i) for i in transit_times]
        peak_list = [peak_list[-1]] #if you entered it twice

    print ("Transits you have entered:    {}   \n".format(str(peak_list))[1:-1])
    print ("Check that these are the transits that you want")
    
    
    # END OF INTERACTIVE PART OF CODE

    #  -----  BREW  ------
    brew.brew_LATTE(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, show = noshow)



def interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, noshow):

    alltime_list, allline, start_sec, end_sec, in_sec, X1_list, X1flux_list,  X4_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI_interact(indir, sectors, sectors_all, tic, save = False)
    
    def close(event):
        if (np.sum(aperture) > 0 ) * (np.sum(aperture2) > 0 ):
            plt.close('all')
        else:
            print ("Must select at least one pixel per aperture!")
            ax[2].text(1.1,-0.29, "Select at least one pixel per aperture!", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            fig.canvas.draw()

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
        
        # place a random aperture in the centre to start off with - just as an example.
        for i in range(int(len(aperture)/2 -1), int(len(aperture)/2 + 1)):
            for j in range(int(len(aperture)/2 -1), int(len(aperture)/2 + 1)):
                mask.append((i,j))

        def onclick(event):
            global mask
            global mask2
            global aperture
            global aperture2

            events = ((int(event.xdata+0.5), int(event.ydata+0.5)))

            if event.inaxes in [ax[1]]:
                [p.remove() for p in reversed(ax[1].patches)]

                if (len(mask) > 0) and (events in list(mask)): # if the square has already been selected, get rid of it. 
                    mask = [x for x in mask if x != events] 
                else:
                    mask.append(events) # otherwise it's a new event and it should be added.
            
                sqcol = '#ffffee'
                alpha = 0.5

                for pixel in mask:
                    m = int(pixel[0])
                    n = int(pixel[1])
                    r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
                    ax[1].add_patch(r)

                fig.canvas.draw()

            if event.inaxes in [ax[2]]:
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

        sqcol = '#ffffee'
        alpha = 0.5


        # ----------
        # start off with the plotting of the larger aperture
        for pixel in mask:
            m = int(pixel[0])
            n = int(pixel[1])
            r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
            ax[1].add_patch(r)

        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)

        for coord in mask:
            aperture[coord] = True

        ax[0].plot(extract_LC(aperture)[0], extract_LC(aperture)[1],marker='o',color = '#054950', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#054950')
    
        # ----------

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
        
        ebx = plt.axes([0.73, 0.025, 0.13, 0.06])
        exit = Button(ebx, 'Close', color='orange')
        exit.on_clicked(close)
        
        plt.show()
        
        # --------------
        # have now chosen an aperture size - now need to extract the LC for the smaller aperture size
        
        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
        for coord in mask[:-1]:
            aperture[coord] = True

        apmask_list.append(aperture)

        if (np.sum(aperture2) == 0):
            target_mask_small = aperture
        else:
            target_mask_small = aperture2

        if (np.sum(aperture2 == 0)) * (np.sum(aperture == 0)):
            print ("WARNIGN: you have not selected any apertures!")

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
        
        T_dur = 0.7  #The transit duration - may need to change!! 
        
        nmed = int(48*3*T_dur)
        nmed = 2*int(nmed/2)+1 # make it an odd number 
        ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
        # first number (nmed) is three time transit durations, the second quite small (10,20 )
        
        l = np.isfinite(ff)
        
        g = interp1d(alltime[l],ff[l],bounds_error=False,fill_value=np.nan)
        ff = g(alltime)
        
        fr = fr_inj / ff
        
        
        # --- do some sigma clipping to make the LC look better ---- 
        MAD = median_absolute_deviation(fr)
        madrange = (5 * MAD * 1.456)
        ymask = (fr < 1 + madrange) * (fr > 1 - madrange) 
        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  


        fix, ax = plt.subplots(3,1,figsize=(16,12))
        
        ax[0].plot(alltime,fr_inj, '.',label = 'Uncorrected')
        ax[0].plot(alltime,ff,'.',label = 'Model fit')
        ax[0].legend(fontsize = 16, loc = 1)
        
        ax[1].plot(alltime[~ymask],fr[~ymask], 'k.', markersize = 10, label = 'Clipped')
        ax[1].plot(alltime[ymask],fr[ymask], 'r.',label = 'Corrected')
        ax[1].legend(fontsize = 16, loc = 1)
        
        ax[2].plot(alltime[ymask],fr[ymask], '.', color = 'navy', label = 'Clipped + corrected')
        ax[2].legend(fontsize = 16, loc = 1)
        
        ax[2].set_xlabel("Time", fontsize = 16)
        ax[0].set_ylabel("Flux", fontsize = 16)
        ax[1].set_ylabel("Flux", fontsize = 16)
        ax[2].set_ylabel("Flux", fontsize = 16)
        
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


# ------ interact FFI -------
def interact_LATTE_FFI(tic, indir, sectors_all, sectors, ra, dec, noshow, FFIap = True):
    '''
    Function to run the Interactive LATTE code using the matplotlib interactive tool.
    Calls the plot where the transit-event times can be identifies and the plotting/modeling options specified.
    
    Parameters
    ----------
    tic  :   str
        target TIC ID
    indir  :  str
        path to directory where all the plots and data will be saved. 
    sectors_all  :   list
        all the sectors in which the target has been/ will be observed
    sectors  :  list
        the sectors which will be analysed

    Returns
    -------
        runs the brew_LATTE code...
    '''
    if FFIap == True:
        alltime0, allflux_list, allflux_small, allflux0, allline, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, noshow)
    
    else:
        print ("Start data download.....", end =" ")
        alltime0, allflux_list, allflux_small, allflux0, allline, allfbkg, allfbkg_t,start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI(indir, sectors, sectors_all, tic, save = True)
        print ("Done.\n")
    
    plt.close('all')
    # -------------------------
    # Plot the interactive plot
    # -------------------------
    # allflux is the corrected one - these arent actually binned but they are 30 mins cadence.
    

    # ------- median absolute deviation in order to determine the clipping
    MAD = median_absolute_deviation(allflux0)
    madrange = (5 * MAD * 1.456)
    
    # make a mask for the points we want to get rid of
    ymask = (allflux0 < 1 + madrange) * (allflux0 > 1 - madrange)
    #  ------

    # need to tidy this code up a bit later !!
    alltime = np.array(alltime0)[ymask]
    allflux = np.array(allflux0)[ymask]

    alltimebinned = alltime
    allfluxbinned = allflux

    fig, ax = plt.subplots(2, 1, figsize=(10,7))
    plt.tight_layout()
    
    # Adjust tbplots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.25)
    
    fluxmin = np.nanmin(allfluxbinned)
    fluxmax = np.nanmax(allfluxbinned)
    
    # Draw the initial plot
    # The 'line' variable is used for modifying the line later
    def cutout(transit):
        mask_binned = (np.array(alltimebinned) > transit-2) & (np.array(alltimebinned) < transit+2)
        mask = (np.array(alltime) > transit-2) & (np.array(alltime) < transit+2)
        
        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux), np.array(alltimebinned)[mask_binned], np.array(allfluxbinned)[mask_binned], np.array(alltimebinned), np.array(allfluxbinned)]
    
    transit = np.nanmean(alltimebinned)
    binfac = 7
    
    [line_full] = ax[0].plot(alltime, allflux , marker='o',lw = 0, markersize = 4, color = '#003941', alpha = 0.8, label = 'unbinned', markerfacecolor='#003941')
    #[line_full_binned] = ax[0].plot(binning(binfac)[0], binning(binfac)[1],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, label = 'binning = 7', markerfacecolor='k')
    
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',lw = 0, markersize = 4, color = '#003941', alpha = 0.8, label = 'unbinned', markerfacecolor='#003941')
    [line_binned] =  ax[1].plot(cutout(transit)[4], cutout(transit)[5],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, label = 'binning = 7', markerfacecolor='k')
    
    # Define an axes area and draw a slider in it
    transit_slider_ax  = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=transit, color='teal')
    
    scale_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
    scale_slider = Slider(scale_slider_ax, 'Y-Axis Scale', 0.99, 1.01, valinit=1, color='silver')

    
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])
    
    ax[1].set_xlim([np.nanmean(alltime)-2, np.nanmean(alltime)+2])
    ax[1].set_ylim([fluxmin, fluxmax])
    
    # Define an action for modifying the line when any slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])
    
        line_binned.set_xdata(cutout(transit_slider.val)[4])
        line_binned.set_ydata(cutout(transit_slider.val)[5])
    
        fig.canvas.draw_idle()
    
    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)
    
    def update_axis(val):   
        ax[1].set_xlim([transit_slider.val - 2,transit_slider.val + 2])
        
        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)
    
    def update_yaxis(val):  
    
        med = 1
        diff = abs(med - (fluxmin * scale_slider.val))
    
        ax[0].set_ylim([med - diff ,med + diff])
        ax[1].set_ylim([med - diff ,med + diff])
    
    
    transit_slider.on_changed(update_axis)
    scale_slider.on_changed(update_yaxis)
    transit_slider.on_changed(sliders_on_changed)

    # Determine whether to save the values the plots or not left, bottom, width, height

    var_ax = fig.add_axes([0.025, 0.3, 0.1, 0.15])
    save_var = CheckButtons(var_ax, ('Simple', 'BLS', 'model', 'Save', 'DVR'), (False, False, False, True, False))
    
    simple = False
    BLS = False
    model = False
    save = True
    DV = False

    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        BLS = status[1]
        model = status[2]
        save = status[3]
        DV = status[4]


    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf
    
    ax[0].tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = 12)
    ax[0].vlines(allline, minf-1,minf + height*0.3 , colors = 'r', label = "Momentum Dump")
    
    ax[1].tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("BJD-2457000", fontsize = 12)
    ax[1].set_ylabel("Normalised Flux", fontsize = 12)
    ax[1].vlines(allline, minf-0.5,minf-0.5 + height*0.3 , colors = 'r', label = "Momentum Dump")
    
    #plt.text(4.5, -3.15, "(Press Enter, then close window)", fontsize=10, verticalalignment='center')
    #plt.text(0.05, 1.1, "Binning Factor", fontsize=10, verticalalignment='center')
    
    initial_text = ""
    
    transit_times = []

    def submit(text):
        ydata = eval(text)
        transit_times.append(ydata)
    
    axbox = plt.axes([0.25, 0.04, 0.50, 0.04])
    text_box = TextBox(axbox, 'Enter transit-event times', initial=initial_text)
    text_box.on_submit(submit)
    

    ebx = plt.axes([0.77, 0.04, 0.13, 0.04])
    exit = Button(ebx, 'Close', color='orange')

    def close(event):
        plt.close('all')

    exit.on_clicked(close)

    plt.show()
    
    end_status = save_var.get_status()

    simple = end_status[0]
    BLS = end_status[1]
    model = end_status[2]
    save = end_status[3]
    DV = end_status[4]


    if len(transit_times) == 0:
        print ("\n Warning: You can't continue without entering a transit-time.\n")
        print ("Exit.\n")
        raise SystemExit
    
    if type(transit_times[-1]) == tuple:
        peak_list = list(transit_times[-1])
        peak_list = [float(i) for i in peak_list]
    
    else:
        peak_list = [float(i) for i in transit_times]
        peak_list = [peak_list[-1]] #if you entered it twice

    print ("Transits you have entered:    {}   \n".format(str(peak_list))[1:-1])
    print ("Check that these are the transits that you want")
    

    # END OF INTERACTIVE PART OF CODE

    #  -----  BREW  ------

    brew.brew_LATTE_FFI(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime0, allflux_list, allflux_small, allflux0, allline, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, show = noshow)

# -----------------------------
# Download the data acess files 
# -----------------------------
# The Functions needed to get the files that know how to acess the data


def data_files(indir):
    '''
    Function to download all of the data that we want to the local computer.
    '''

    if not os.path.exists("{}/data/tesscurl_sector_all_lc.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_lc.sh".format(indir),'w') as f:
            f.write("#all LC file links")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will download all of the available sectors starting with sector 1")
        
    else:
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_lc.sh > {0}/data/temp.txt'.format(indir))
        
        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]
        
        first_sec = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here
            
    
    if not os.path.exists("{}/data/tesscurl_sector_all_tp.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_tp.sh".format(indir),'w') as f:
            f.write("#all LC file links")
        first_sec_tp = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will download all of the available sectors starting with sector 1")
        
    else:
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_tp.sh > {0}/data/temp_tp.txt'.format(indir))
        
        with open("{}/data/temp_tp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]
        
        first_sec_tp = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here
            
    
    for sec in range(first_sec+1,27): # 26 because that's how many TESS sectors there will be in total
    
        LC_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh".format(sec)
        r_LC = requests.get(LC_url) # create HTTP response object
            
        if r_LC.status_code == 404:
            print ("Data only available up to Sector {} -- try downloading more data later".format(sec))
            break
    
            
        with open("{}/data/tesscurl_sector_all_lc.sh".format(indir), 'ab') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_LC.content)
                print("finished adding files for sector {}".format(sec))
                #write the contents of the response (r.content)
                # to a new file in binary mode.    
    
                
        with open("{}/data/tesscurl_sector_{}_lc.sh".format(indir,sec), 'wb') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_LC.content)
    
                
    for sec in range(first_sec_tp+1,27): # 26 because that's how many TESS sectors there will be in total
    
        TP_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_tp.sh".format(sec)
        r_TP = requests.get(TP_url) # create HTTP response object 
            
        if r_TP.status_code == 404:
            print ("TP data only available up to Sector {} -- try downloading more data later".format(sec))
            break
    
        with open("{}/data/tesscurl_sector_all_tp.sh".format(indir), 'ab') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_TP.content)
                #rint("finished adding sector {} for TP".format(sec))
                #write the contents of the response (r.content)
                # to a new file in binary mode.    
    
                
        with open("{}/data/tesscurl_sector_{}_tp.sh".format(indir,sec), 'wb') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_TP.content)

def nn_files(indir):
    '''
    Function to download all of the TPF data that we want to the local computer.
    '''
    if not os.path.exists("{}/data/all_targets_list.txt".format(indir)):
        with open("{}/data/all_targets_list.txt".format(indir),'w') as f:
            f.write("#all targets file links")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will download all of the available sectors starting with sector 1")
    
    else:
        files = np.sort(glob('{}/data/all_targets_S*'.format(indir)))
        
        exist = []
        for f in files:
            exist.append(int(f[-9:-7]))  # get the first sector number (last that has already been downloaded)
            
        first_sec = (np.max(exist))
        
            
    for sec in range(first_sec+1,27): # 26 because that's how many TESS sectors there will be in total
    
        if sec < 10:
            download_sector = "00{}".format(sec)
        else:
            download_sector = "0{}".format(sec)
        
        target_list = "https://tess.mit.edu/wp-content/uploads/all_targets_S{}_v1.txt".format(download_sector)
    
        r_target_list = requests.get(target_list) # create HTTP response object
        
        if r_target_list.status_code == 404:
            print ("Target lists only available up to Sector {} -- try downloading more data later".format(sec))
            break
        
        
        with open("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), 'wb') as f:
            '''
            Saving recieved content as a png file in binary format
            '''
            f.write(r_target_list.content)

            
        with open("{}/data/all_targets_list.txt".format(indir), 'ab') as f:
            '''
            Saving recieved content as a png file in binary format
            '''
            if sec == 1:
                f.write(r_target_list.content)
    
            else:
                start = str(r_target_list.content).find('t  Dec')
                f.write(r_target_list.content[start-3:])
                print("finished adding TP sector {}".format(sec))

def TOI_TCE_files(indir):
    '''
    Function to download the files that list all the known TOI's and TCEs.
    This is useful to display on the DV report.
    '''

    # ------ TOIs ------
    TOI_url = "https://tev.mit.edu/data/collection/193/csv/5/"
    r_TOI = requests.get(TOI_url) # create HTTP response object
        
    if r_TOI.status_code == 404:
        print ("Can't download the TOI list at the moment. Has the URL changed?")
    
    with open("{}/data/TOI_list.txt".format(indir),'wb') as f:
           f.write(r_TOI.content)


    # ------ TCE ------

    if not os.path.exists("{}/data/tesscurl_sector_all_dv.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_dv.sh".format(indir),'w') as f:
            f.write("#all LC file links")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will download all of the DV report links for all available sectors starting with sector 1... ")
        
    else:
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_dv.sh > {0}/data/temp.txt'.format(indir))
        
        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]
        
        first_sec = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here
          

    for sec in range(first_sec+1,27): # 27 because that's how many TESS sectors there will be in total
    
        TCE_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_dv.sh".format(sec)
        r_TCE = requests.get(TCE_url) # create HTTP response object
            
        if r_TCE.status_code == 404:
            print ("DV report data only available up to Sector {} -- try downloading more data later".format(sec))
            break
    
        with open("{}/data/tesscurl_sector_all_dv.sh".format(indir), 'ab') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_TCE.content)
                print("finished adding DV links for sector {}".format(sec))

def momentum_dumps_info(indir):
    '''
    function to the create a list of all of the momentum dump times for each sector - only needed in the FFIs
    '''
    print ("store the times of the momentum dumps for each sector - only needed when looking at the FFIs")

    if not os.path.exists("{}/data/tess_mom_dumps.txt".format(indir)):
        with open("{}/data/tess_mom_dumps.txt".format(indir),'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['sec', 'time'])
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will determine the times of the momentum dumps for all sectors and store them... ")
        
    else:
        os.system('tail -n 1 {0}/data/tess_mom_dumps.txt > {0}/data/temp.txt'.format(indir))
        
        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]
            first_sec = int(string[0:2])
            
    for sec in range(first_sec+1,27): # 27 because that's how many TESS sectors there will be in total
        
        try:
            lc_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_lc.sh'.format(indir, str(sec)), dtype = str)
            lcfile = (lc_sec[0][6])
            
            response = requests.get(lcfile)
            lchdu  = pf.open(response.url) # this needs to be a URL - not a file
            
                                   #Open and view columns in lightcurve extension
            lcdata = lchdu[1].data                 
            quality = lcdata['QUALITY']
            time    = lcdata['TIME']
        
            mom_dump_mask = np.bitwise_and(quality, 2**5) >= 1                               
            
            momdump = (list(time[mom_dump_mask]))
            sector = list([sec] * len(momdump))
                                   
            with open("{}/data/tess_mom_dumps.txt".format(indir), 'a') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(zip(sector,momdump))
            
        except:
            continue


# -----------------------

def tess_point(indir,tic):
    '''
    Use tess-point to find out in what sectors the dubject appears in.
    
    Parameters
    ----------
    indir : str
        path to where the files will be saved.
    tic : str
        TIC (Tess Input Catalog) ID of the target

    Returns
    -------
    outSec  : list
        list of all the sectors in which this subejct has been /  will be observed  

    '''

    #target_list = np.genfromtxt("{}/data/all_targets_list.txt".format(indir), dtype = str, usecols = (0,4,5)) 
    #header = ['tic', 'ra', 'dec']
    #pd_list = pd.DataFrame(target_list, columns=header)
    #
    ##print (tic)
    ##print (type(tic))
    ##print("print this valie: {}".format(pd_list.loc[pd_list['tic'] == str(55525572), 'ra']))

    #ra = pd_list.loc[pd_list['tic'] == str(tic), 'ra'].values[0]
    #dec = pd_list.loc[pd_list['tic'] == str(tic), 'dec'].values[0]
    #ra = float(ra)
    #dec = float(dec)
    
    #_, _, _, outSec, _, _, _, _, _ = tess_stars2px_function_entry(tic, ra, dec)
    
    if not exists('{}/tesspoint'.format(indir)):
        os.makedirs('{}/tesspoint'.format(indir))    

    os.system('python3 -m tess_stars2px -t {} > {}/tesspoint/{}_tesspoint.txt'.format(tic,indir,tic))

    df = pd.read_csv('{}/tesspoint/{}_tesspoint.txt'.format(indir,tic), comment = '#', delimiter = '|', names = ['TIC','RA','Dec','EclipticLong','EclipticLat','Sector','Camera','Ccd','ColPix', 'RowPix'])

    return list(df['Sector']), float(df['RA'][0]), float(df['Dec'][0])


def peak_sec(in_sec,start_sec, end_sec, peak_list):
    '''
    Use tess-point to find out in what sectors the dubject appears in.
    
    Parameters
    ----------
    in_sec : list
        the sectors that we are looking at.
    start_sec : list
        the start time of these sectors
    end_sec : list
        the end time of these sectors
    peak_list : list
        list of the marked peaks

    Returns
    -------
    peak_sec  : list
        list of the sectors in which the peaks appear.

    '''

    peak_sec = []
    for peak in peak_list:
        for n,start in enumerate(start_sec):
    
            if peak >= start[0] and peak <= end_sec[n][0]:
                
                peak_sec.append(in_sec[n])
    
    return list(set(peak_sec))

def nn_ticids(indir, peak_sec, tic):
    '''
    To find the TIC IDs of the 6 nearest neighbour stars.
    
    Parameters
    ----------
    indir : str
        path to where the files will be saved.
    peak_sec  : list
        list of the sectors in which the peaks appear.
    tic : str
        TIC (Tess Input Catalog) ID of the target

    Returns
    -------
    ticids  : list
        TIC IDs of the 6 TESS target pixels file stars that are closest to the target.
    target_ra  : float
        Right Ascension of the target
    target_dec  : float
        Declination of the target

    '''    
    neighbours_sector = peak_sec[0]

    if neighbours_sector < 10:
        download_sector = "00{}".format(neighbours_sector)
    else:
        download_sector = "0{}".format(neighbours_sector)


    #tic_list = pd.read_csv("{}/all_targets_S{}_v1.txt".format(indir, download_sector), sep = ',', comment='#').sort_values(['Camera', 'RA', 'Dec']).reset_index()
    
    tic_list = pd.read_table("{}/data/all_targets_S{}_v1.txt".format(indir,download_sector), sep='\t', lineterminator='\n', comment = '#', names = ['TICID', 'Camera', 'CCD', 'Tmag', 'RA', 'Dec']).sort_values(['Camera', 'RA', 'Dec']).reset_index()
    
    # the target
    target = tic_list.loc[tic_list['TICID'] == float(tic)]

    
    tic_idx = target.index[0]  # get the index of the target in the list

    target_ra = float(target['RA']) 
    target_dec = float(target['Dec'])
    
    # the closest targets
    tic_list_close = tic_list[tic_idx - 100:tic_idx + 101]
    
    # Calculated the angular separation to the stars to find the nearest neighbours
    def star_sep(row, ra, dec):
        
        ra2 = float(row['RA'])
        dec2 = float(row['Dec'])

        c1 = SkyCoord(ra*u.degree, dec*u.degree)
        c2 = SkyCoord(ra2*u.degree, dec2*u.degree)
        sep = c1.separation(c2)
        
        return sep.arcminute
    
    tic_list_close['dist'] = tic_list_close.apply(star_sep, args = (target_ra,target_dec), axis=1)
    
    # get the tics that are closest 
    closest_tic = tic_list_close.sort_values('dist')[0:6]
    ticids = closest_tic['TICID'].tolist()
    distance = closest_tic['dist'].tolist()

    return ticids, distance, target_ra, target_dec


# -----------------------
# download the data 
# -----------------------

# The functions to download the actual data
def download_data(indir,sector, tic, binfac = 5):
    '''
    Download the LCs for the target star for all the indicated sectors
    
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


    sb.set(style='ticks')
    sb.set_color_codes('deep')
    
    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
    
    dwload_link = []
    
    if sector == 'all':
        # locate the file string
        lc_all = np.genfromtxt('{}/data/tesscurl_sector_all_lc.sh'.format(indir), dtype = str)
    
        for i in lc_all:
            if str(tic) in str(i[6]):
                dwload_link.append(i[6])
    
    else:
        future_sectors = []
        # locate the file string
        for s in sector:
            try:
                lc_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_lc.sh'.format(indir, str(s)), dtype = str)
                
                for i in lc_sec:
                    if str(tic) in str(i[6]):
                        dwload_link.append(i[6])
            except:
                future_sectors.append(s)
    
        if len(future_sectors):
            print ("In the future, this TIC {} will be observed in sector(s) {}".format(tic, future_sectors))
    
    if len(dwload_link) == 0:
        print ("This TIC was not observed in Sector(s):   {}   .Try again with different sectors.".format(sector))
        raise SystemExit
        #exit[0]

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
    allline = []
    
    for lcfile in dwload_link:
        
        response = requests.get(lcfile)
    
        # open the file using the response url  
        lchdu  = pf.open(response.url) # this needs to be a URL - not a file
        
        #Open and view columns in lightcurve extension
        lcdata = lchdu[1].data
        lchdu[1].columns

        sapflux = lcdata['SAP_FLUX']
        f02 = lcdata['PDCSAP_FLUX']
        f02_err = lcdata['PDCSAP_FLUX_ERR']
        quality = lcdata['QUALITY']
        time    = lcdata['TIME']
        f0     = lcdata['SAP_FLUX']
        fbkg     = lcdata['SAP_BKG']
        
        med = np.nanmedian(f02)
        f1 = f02/med
        f1_err = f02_err/med
        
        x1      = lcdata['MOM_CENTR1']  # CCD column position of target’s flux-weighted centroid 
        x1      -= np.nanmedian(x1)
        y1      = lcdata['MOM_CENTR2']  
        y1      -= np.nanmedian(y1)
        x2      = lcdata['POS_CORR1'] # The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects.
        x2      -= np.nanmedian(x2)
        y2      = lcdata['POS_CORR2']
        y2      -= np.nanmedian(y2)
        l       = (quality>0)
        l2      = (quality<=0)
        
        sec     = int(lchdu[0].header['SECTOR'])

        tessmag = lchdu[0].header['TESSMAG']
        teff    = lchdu[0].header['TEFF']
        srad    = lchdu[0].header['RADIUS']

        flux     = lcdata['SAP_FLUX']

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

        #bad_bits = np.array([1,2,3,4,5,6,8,10,12])
        #value = 0
        #for v in bad_bits:
        #    value = value + 2**(v-1)

        #bad_data = np.bitwise_and(quality, value) >= 1

        #fluxcent_col = lcdata['MOM_CENTR1']
        #fluxcent_row = lcdata['MOM_CENTR2']

        mom_dump = np.bitwise_and(quality, 2**5) >= 1

        alltime.append(list(time)) #[~bad_data]
        allflux.append(list(f1)) #[~bad_data])fly
        allflux_err.append(list(f1_err))
        allline.append(list(time[mom_dump]))
        
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
    

    alltime = [val for sublist in alltime for val in sublist]
    allflux = [val for sublist in allflux for val in sublist]
    allflux_err = [val for sublist in allflux_err for val in sublist]
    allline = [val for sublist in allline for val in sublist]

    alltimebinned = [val for sublist in alltimebinned for val in sublist]
    allfluxbinned = [val for sublist in allfluxbinned for val in sublist]
    
    allx1 = [val for sublist in allx1 for val in sublist]
    allx2 = [val for sublist in allx2 for val in sublist]
    ally1 = [val for sublist in ally1 for val in sublist]
    ally2 = [val for sublist in ally2 for val in sublist]
    alltimel2 = [val for sublist in alltimel2 for val in sublist]
    
    allfbkg = [val for sublist in allfbkg for val in sublist]
   
    return alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltimel2, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad


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
        
        t = tpf.time[lkeep]
        
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


def download_data_FFI(indir,sector, sectors_all, tic, save = False):
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

    allfbkg = []
    allfbkg_t = []
    
    start_sec = []
    end_sec = []
    in_sec = []
    
    alltime_list = []
    allflux = []
    allflux_flat = []
    allflux_small = []
    allline = []

    X1_list = []
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
        arrshape_list.append(X1.shape)

        #identify the target mask - this will need to be imporved - a bit 'hacky' at the moment. 

        val = 1
        for i in range(0,50):
            
            target_mask = tpf.create_threshold_mask(threshold=val, reference_pixel='center')
            #print (np.sum(target_mask))
            #print ("val {}".format(val))
            if (np.sum(target_mask) < 9) and (np.sum(target_mask) > 7):
                threshhold = val
                break 
            else:
                if np.sum(target_mask) < 6:
                    if val < 1:
                        val -= 0.05
                    else:
                        val -= 0.5
                elif np.sum(target_mask) > 9:
                    if val > 20:
                        val += 5
                    elif val > 50:
                        val += 20
                    else:
                        val += 0.5
               

        val_small = val
        
        for i in range(0,40):
            
            target_mask_small = tpf.create_threshold_mask(threshold=val_small, reference_pixel='center')
            
            if (np.sum(target_mask_small) < 4) and (np.sum(target_mask_small) > 2):
                threshhold = val_small
                break 
            else:
                #print (np.sum(target_mask))
                
                if np.sum(target_mask_small) < 2:
        
                    if val_small < 1:
                        val_small -= 0.05
                    else:
                        val_small -= 0.5
                elif np.sum(target_mask_small) > 4:
        
                    if val_small > 20:
                        val_small += 5
                    elif val_small > 50:
                        val_small += 20
                    else:
                        val_small += 0.5
        
        plt.figure(figsize=(5,5))
        mask_plot = tpf.plot(aperture_mask=target_mask, mask_color='k')
        plt.savefig('{}/{}/{}_mask.png'.format(indir, tic, tic), format='png')
        plt.clf()
        plt.close()
        
        plt.figure(figsize=(5,5))
        mask_plot = tpf.plot(aperture_mask=target_mask_small, mask_color='k')
        plt.savefig('{}/{}/{}_mask_small.png'.format(indir, tic, tic), format='png')
        plt.clf()
        plt.close()

        bkg = X1
        bkg = bkg.mean(axis = 0)
        bkg_list.append(bkg)

        # -----
        apmask_list.append(target_mask)  # store the target mask used. 
        # -----

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


        T_dur = 0.7  #The transit duration - may need to change!! 
        
        nmed = int(48*3*T_dur)
        nmed = 2*int(nmed/2)+1 # make it an odd number 
        ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
        # first number (nmed) is three time transit durations, the second quite small (10,20 )
        
        l = np.isfinite(ff)
        
        g = interp1d(alltime[l],ff[l],bounds_error=False,fill_value=np.nan)
        ff = g(alltime)
        
        fr = fr_inj / ff
        
        # ------- median absolute deviation in order to determine the clipping
        MAD = median_absolute_deviation(fr)
        madrange = (5 * MAD * 1.456)
        
        # make a mask for the points we want to get rid of
        ymask = (fr < 1 + madrange) * (fr > 1 - madrange)
        #  ------

        if save == True:

            fix, ax = plt.subplots(3,1,figsize=(16,12))
            
            ax[0].plot(alltime,fr_inj, '.',label = 'Uncorrected')
            ax[0].plot(alltime,ff,'.',label = 'Model fit')
            ax[0].legend(fontsize = 16, loc = 1)
            
            ax[1].plot(alltime[~ymask],fr[~ymask], 'k.', markersize = 10, label = 'Clipped')
            ax[1].plot(alltime[ymask],fr[ymask], 'r.',label = 'Corrected')
            ax[1].legend(fontsize = 16, loc = 1)
            
            ax[2].plot(alltime[ymask],fr[ymask], '.', color = 'navy', label = 'Clipped + corrected')
            ax[2].legend(fontsize = 16, loc = 1)
            
            ax[2].set_xlabel("Time", fontsize = 16)
            ax[0].set_ylabel("Flux", fontsize = 16)
            ax[1].set_ylabel("Flux", fontsize = 16)
            ax[2].set_ylabel("Flux", fontsize = 16)
            
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

        # add all the information to lists that are returned at the end. 
        tpf_filt_list.append(X4.reshape(tpf.flux[lkeep,:,:].shape))

        in_sec.append(sec)

        alltime_list.append(list(alltime))
        allflux_flat.append(list(fr))
        allflux.append(list(flux))
        allflux_small.append(list(flux_small))

        allline.append(list((mom_df.loc[mom_df['sec'] == int(sec)])['time']))
        
        start_sec.append([alltime[0]])
        end_sec.append([alltime[-1]])
    
        X1_list.append(X1) #  not corrected
        X4_list.append(X4) #  PCA corrected
        t_list.append(np.array(alltime))  # add it here because this list isn't flattened but the other one is


    alltime_list = [val for sublist in alltime_list for val in sublist]
    allflux_flat = [val for sublist in allflux_flat for val in sublist]
    allflux_small = [val for sublist in allflux_small for val in sublist]
    allflux = [val for sublist in allflux for val in sublist]
    allline = [val for sublist in allline for val in sublist]
    allfbkg = [val for sublist in allfbkg for val in sublist]
    allfbkg_t = [val for sublist in allfbkg_t for val in sublist]
   
    del mom_df

    return alltime_list, allflux, allflux_small, allflux_flat, allline, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list


def download_data_neighbours(indir, sector, tics, distance, binfac = 5):
    
    sb.set(style='ticks')
    sb.set_color_codes('deep')
    
    try:
        lc_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_lc.sh'.format(indir, str(sector)), dtype = str)
    except:
        print ("Sector {}  has not yet been osberved - come back to this target later.".format(sector))

    alltime_neighbours = []
    allflux_neighbours = []
    allline_neighbours = []
    alltimebinned_neighbours = []
    allfluxbinned_neighbours = []
    outtics = []

    dwload_link = []
    
    for tic in tics:
        for i in lc_sec:  
            try:
                if str(tic) in str(i[6]):
                    dwload_link.append(i[6])
            except:
                print("{} was not observed in Sector {}".format(tic, sec))

    alltimebinned = []
    allfluxbinned = []

    start_sec = []
    end_sec = []
    
    alltime = []
    allflux = []
    allline = []

    tessmag_list = []
    
    for num,lcfile in enumerate(dwload_link):
        
        print ("Downloading nearest neighbour   {}   of   {}....".format(num + 1, len(dwload_link), end ='' ))

        response = requests.get(lcfile)
    
        # open the file using the response url  
        lchdu  = pf.open(response.url) # this needs to be a URL - not a file
        outtics.append(int(lchdu[0].header['TICID']))

        #Open and view columns in lightcurve extension
        lcdata = lchdu[1].data
        lchdu[1].columns
    
        sapflux = lcdata['SAP_FLUX']
        f02 = lcdata['PDCSAP_FLUX']
        quality = lcdata['QUALITY']
        time    = lcdata['TIME']
        f0     = lcdata['SAP_FLUX']
        fbkg     = lcdata['SAP_BKG']
    
        med = np.nanmedian(f02)
        f1 = f02/med
    
        x1      = lcdata['MOM_CENTR1']
        x1      -= np.nanmedian(x1)
        y1      = lcdata['MOM_CENTR2']
        y1      -= np.nanmedian(y1)
        x2      = lcdata['POS_CORR1']
        x2      -= np.nanmedian(x2)
        y2      = lcdata['POS_CORR2']
        y2      -= np.nanmedian(y2)
        l       = (quality>0)
        l2      = (quality<=0)
        
        sec     = int(lchdu[0].header['SECTOR'])
    
        flux     = lcdata['SAP_FLUX']
        
        tessmag = lchdu[0].header['TESSMAG']
        

        # binned data
        N       = len(time)
        n       = int(np.floor(N/binfac)*binfac)
        X       = np.zeros((2,n))
        X[0,:]  = time[:n]
        X[1,:]  = f1[:n]
        Xb      = rebin(X, (2,int(n/binfac)))
    
        time_binned    = Xb[0]
        flux_binned    = Xb[1]
    
        #bad_bits = np.array([1,2,3,4,5,6,8,10,12])
        #value = 0
        #for v in bad_bits:
        #    value = value + 2**(v-1)
    #
        #bad_data = np.bitwise_and(quality, value) >= 1
    #
        #fluxcent_col = lcdata['MOM_CENTR1']
        #fluxcent_row = lcdata['MOM_CENTR2']
    
        mom_dump = np.bitwise_and(quality, 2**5) >= 1
    
        alltime.append(list(time)) #[~bad_data]
        allflux.append(list(f1)) #[~bad_data])
        allline.append(list(time[mom_dump]))
        
        alltimebinned.append(list(time_binned))
        allfluxbinned.append(list(flux_binned))
        

        start_sec.append([time[0]])
        end_sec.append([time[-1]])
        tessmag_list.append(tessmag)
        
        print("Done.")
    

    return alltime, allflux, allline, alltimebinned, allfluxbinned, outtics, tessmag_list, distance

# WITH Lightkurve
def tpf_data(indir, sector, tic):

    
    dwload_link = []
    
    if sector == 'all':
        # locate the file string
        tpf_all = np.genfromtxt('{}/data/tesscurl_sector_all_tp.sh'.format(indir), dtype = str)
    
        for i in tpf_all:
            if str(tic) in str(i[6]):
                dwload_link.append(i[6])
        
    else:

        future_sectors = []
        # locate the file string
        for s in sector:
            try:
                tpf_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_tp.sh'.format(indir, str(s)), dtype = str)
            except:
                future_sectors.append(s)
            for i in tpf_sec:
                try:
                    if str(tic) in str(i[6]):
                        dwload_link.append(i[6])
                except:
                    print("{} was not observed in Sector {}".format(tic, sec))

        if len(future_sectors) > 0:
            print ("In the future, this TIC {} will be observed in sector(s) {}".format(tic, future_sectors))

    if len(dwload_link) == 0:
        print ("This TIC was not observed in Sector(s):   {}   .Try again with different sectors.".format(sector))


    TESS_unbinned_t_l = []
    TESS_binned_t_l = []
    small_binned_t_l= []
    TESS_unbinned_l= []
    TESS_binned_l= []
    small_binned_l= []
    tpf_list = []
 
    dwload_link_tp = dwload_link
    
    for idx,file in enumerate(dwload_link_tp):

        tpf = TessTargetPixelFile(file) # dowload the Target Pixel File
        tpf_list.append(tpf)

        try:

            lc = tpf.to_lightcurve()
            median_image = np.nanmedian(tpf.flux, axis=0)
            smaller_mask = median_image > np.nanpercentile(median_image, 50)
            #smaller_mask2 = median_image > np.nanpercentile(median_image, 75)
            
            # TESS binned
            TESS_unbinned = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask).flatten(window_length=100001)
            TESS_unbinned = TESS_unbinned.remove_outliers(6)
            
            # TESS binned
            TESS_binned = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask).flatten(window_length=100001)
            TESS_binned = TESS_binned.remove_outliers(6).bin(7)
            
            # Use a custom aperture binned
            small_binned = tpf.to_lightcurve(aperture_mask=smaller_mask).flatten(window_length=100001)
            small_binned = small_binned.remove_outliers(6).bin(7)
        
            #small_binned2 = tpf.to_lightcurve(aperture_mask=smaller_mask2).flatten(window_length=100001)
            #small_binned2 = small_binned2.remove_outliers(6).bin(7)
        
            TESS_unbinned_t = TESS_unbinned.time
            TESS_binned_t = TESS_binned.time
            small_binned_t = small_binned.time
             
            # ----------
            TESS_unbinned_t_l.append(TESS_unbinned_t)
            TESS_binned_t_l.append(TESS_binned_t)
            small_binned_t_l.append(small_binned_t)
            
            TESS_unbinned_l.append(TESS_unbinned.flux)
            TESS_binned_l.append(TESS_binned.flux)
            small_binned_l.append(small_binned.flux)
        
        except:
            continue
        
        #small_binned_t2 = small_binned2.time
        #plt.scatter(small_binned_t2, small_binned2.flux, s = 10, color = 'green', label = 'Small Aperture binned')
    
    TESS_unbinned_t_l = [val for sublist in TESS_unbinned_t_l for val in sublist]
    TESS_binned_t_l   = [val for sublist in TESS_binned_t_l for val in sublist]
    small_binned_t_l = [val for sublist in small_binned_t_l for val in sublist]
    
    TESS_unbinned_l = [val for sublist in TESS_unbinned_l for val in sublist]
    TESS_binned_l = [val for sublist in TESS_binned_l for val in sublist]
    small_binned_l = [val for sublist in small_binned_l for val in sublist]
    

    return TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list


# without Lightkurve
def download_data_tpf(indir, peak_sec, peak_list, tic):
    '''
    Download the TPF LCs for the target star for all the indicated sectors. Not using Lightkurve
    
    Parameters
    ----------
    indir : str
        path to where the files will be saved.
    peak_sec  :  list or str
        list of the sectors that have a transit in them. If 'all', all the sectors in whic the target appears will be downloaded
    tic : str
        TIC (Tess Input Catalog) ID of the target
    peak_list  :  int
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

    #print ("peak list {}".format(peak_sec))

    dwload_link_tp = []

    for sec in peak_sec: #the sector that this image is in

        tpf_all = np.genfromtxt('{}/data/tesscurl_sector_{}_tp.sh'.format(indir,sec), dtype = str)
        
        for i in tpf_all:
        
            if str(tic) in str(i[6]):
                dwload_link_tp.append(i[6])
    
    # download each file (i.e. each sector)

    for file in dwload_link_tp:

        tpf = pf.open(file)   # open the file

        for T0 in peak_list:

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


def data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, save = False, show = False):
    
    # make sure that there are no nan value sin the data - they cause everything to crash
    mask_binned = np.isfinite(alltimebinned) * np.isfinite(allfluxbinned)
    mask = np.isfinite(alltime) * np.isfinite(allflux)
    
    alltimebinned = np.array(alltimebinned)[mask_binned]
    allfluxbinned = np.array(allfluxbinned)[mask_binned]
    alltime = np.array(alltime)[mask]
    allflux = np.array(allflux)[mask]
    
    durations = np.linspace(0.05, 0.2, 10)
    model = BoxLeastSquares(alltimebinned, allfluxbinned)
    results = model.autopower(durations, frequency_factor=5.0)

    index = np.argmax(results.power)
    period = results.period[index]
    t0 = results.transit_time[index]
    duration = results.duration[index]
    

    # call the first round of plotting
    plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model, results,period,duration,t0, save = save, show = show)
    
    stats_period = period
    stats_t0 = t0
    stats_depth = model.compute_stats(period, duration, t0)['depth']
    stats_depth_phased = model.compute_stats(period, duration, t0)['depth_phased']
    stats_depth_half = model.compute_stats(period, duration, t0)['depth_half']
    stats_depth_odd = model.compute_stats(period, duration, t0)['depth_odd']
    stats_depth_even = model.compute_stats(period, duration, t0)['depth_even']


    # Find the in-transit points using a longer duration as a buffer to avoid ingress and egress
    in_transit = model.transit_mask(alltimebinned, period, 2*duration, t0)

    
    # Re-run the algorithm, and plot the results
    model2 = BoxLeastSquares(alltimebinned[~in_transit], allfluxbinned[~in_transit])
    results2 = model2.autopower(durations, frequency_factor=5.0)
    
    # Extract the parameters of the best-fit model
    index = np.argmax(results2.power)
    period2 = results2.period[index]
    t02 = results2.transit_time[index]
    duration2 = results2.duration[index]
    
    # call the second round of plotting - once the intitial transit has been removed
    plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model2, results2,period2,duration2,t02, in_transit = in_transit, save = save, show = show)

    stats2_period = period2
    stats2_t0 = t02
    stats2_depth = model2.compute_stats(period2, duration2, t0)['depth']
    stats2_depth_phased = model2.compute_stats(period2, duration2, t0)['depth_phased']
    stats2_depth_half = model2.compute_stats(period2, duration2, t0)['depth_half']
    stats2_depth_odd = model2.compute_stats(period2, duration2, t0)['depth_odd']
    stats2_depth_even = model2.compute_stats(period2, duration2, t0)['depth_even']
    
    return [stats2_period, stats2_t0, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [stats_period, stats_t0, stats2_depth, stats2_depth_phased, stats2_depth_half, stats2_depth_odd, stats2_depth_even]


def data_bls_FFI(tic, indir, alltime, allflux, save = False, show = False):
    
    # make sure that there are no nan value sin the data - they cause everything to crash

    mask = np.isfinite(alltime) * np.isfinite(allflux)
    
    alltime = np.array(alltime)[mask]
    allflux = np.array(allflux)[mask]
    
    durations = np.linspace(0.05, 0.2, 10)
    model = BoxLeastSquares(alltime, allflux)
    results = model.autopower(durations, frequency_factor=5.0)

    index = np.argmax(results.power)
    period = results.period[index]
    t0 = results.transit_time[index]
    duration = results.duration[index]
    
    # call the first round of plotting
    plot_bls_FFI(tic, indir, alltime, allflux, model, results, period, duration, t0, save = save, show = show)
    
    stats_period = period
    stats_t0 = t0
    stats_depth = model.compute_stats(period, duration, t0)['depth']
    stats_depth_phased = model.compute_stats(period, duration, t0)['depth_phased']
    stats_depth_half = model.compute_stats(period, duration, t0)['depth_half']
    stats_depth_odd = model.compute_stats(period, duration, t0)['depth_odd']
    stats_depth_even = model.compute_stats(period, duration, t0)['depth_even']

    # Find the in-transit points using a longer duration as a buffer to avoid ingress and egress
    in_transit = model.transit_mask(alltime, period, 2*duration, t0)

    
    # Re-run the algorithm, and plot the results
    model2 = BoxLeastSquares(alltime[~in_transit], allflux[~in_transit])
    results2 = model2.autopower(durations, frequency_factor=5.0)
    
    # Extract the parameters of the best-fit model
    index = np.argmax(results2.power)
    period2 = results2.period[index]
    t02 = results2.transit_time[index]
    duration2 = results2.duration[index]
    
    # call the second round of plotting - once the intitial transit has been removed
    plot_bls_FFI(tic, indir, alltime, allflux, model2, results2,period2,duration2,t02, in_transit = in_transit, save = save, show = show)

    stats2_period = period2
    stats2_t0 = t02
    stats2_depth = model2.compute_stats(period2, duration2, t0)['depth']
    stats2_depth_phased = model2.compute_stats(period2, duration2, t0)['depth_phased']
    stats2_depth_half = model2.compute_stats(period2, duration2, t0)['depth_half']
    stats2_depth_odd = model2.compute_stats(period2, duration2, t0)['depth_odd']
    stats2_depth_even = model2.compute_stats(period2, duration2, t0)['depth_even']
    
    return [stats2_period, stats2_t0, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [stats_period, stats_t0, stats2_depth, stats2_depth_phased, stats2_depth_half, stats2_depth_odd, stats2_depth_even]

# -----------------------
# plots
# -----------------------

def plot_nn(tic, indir,alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, peak_list, outtics, tessmag_list, distance, save = False, show = False):
    
    '''
    Plot the lighcurves of the 6 nearest neighbours to the target. 

    Distance to the targets is calculated using the RA and DEC of the targte and adding them in quadrature 
    
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
    alltime_nn  :  list
        times for all the nearest neighbours (not binned)
    allflux_nn  :  list
        normalized flux times for all the nearest neighbours (not binned)
    alltimebinned_nn  :  list
        binned time for all the nearest neighbours
    allfluxbinned_nn  :  list
        normalized binned flux for all the nearest neighbours
    peak_list  :  list
        list of all the marked transits
    outtics  :  list
        the tic IDs of the 6 nearest neighbours
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID

    Returns
    -------
        Plot of the 6 nearest neighbours. The vertical line indicated the location(s) of the marked transit. 

    '''

    fig, ax = plt.subplots(len(alltime_nn), 1, figsize=(13,8), sharex=True, gridspec_kw={'hspace': 0})
    plt.tight_layout()

    colors = ['r', 'darkorange', 'gold', 'seagreen', 'royalblue', 'navy','magenta' ]
    
    for i in range(0,len(alltime_nn)):
    
        for line in (peak_list):
            ax[i].axvline(line, color = 'k', linewidth = 2.2, alpha = 1, linestyle = '-')
        if str(outtics[i]) == str(tic): 
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "*** {}  Tmag = {:3f} ***".format(tic, tessmag_list[i]), marker = '.', ms = 2, linewidth = 0)

        else:
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "{}  Tmag = {:3f}   d = {:3f} arcsecs".format(outtics[i], tessmag_list[i], distance[i]), marker = '.', ms = 2, linewidth = 0)
        
        ax[i].plot(alltimebinned_nn[i], np.array(allfluxbinned_nn[i]), color = 'k', marker = '.', ms = 1, linewidth = 0)
              
        ax[i].legend(loc = 1)
    
    ax[0].set_title("LCs of Nearby Stars")
    ax[0].set_title("LCs of Nearby Stars")
    ax[len(alltime_nn) - 1].set_xlabel("Time (BJD-2457000)")
    ax[int(len(alltime_nn)/2)].set_ylabel("Normalised Flux")

    if save == True:
        plt.savefig('{}/{}/{}_nearest_neighbours.png'.format(indir, tic, tic), format='png')

    plt.xlim(np.nanmin(alltime_nn), np.nanmax(alltime_nn))

    if show == True:
        plt.show()
    else:
        plt.close()


def plot_cutout(image):
    """
    Plot image cut out of the target. 
    """
    plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r,
           vmax = np.percentile(image, 92),
           vmin = np.percentile(image, 5))

    plt.grid(axis = 'both',color = 'white', ls = 'solid')


def plot_centroid(tic, indir,alltime12, allx1, ally1, allx2, ally2, peak_list, save = False, show = False):
    '''
    Plot the x and y centroids around the time(s) of the marked transit.

    Download the LCs for the target star for all the indicated sectors
    
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
    alltimel2  :  list
        time used for the x and y centroid position plotting
    allx1  :  list
        CCD column position of target’s flux-weighted centroid. In x direction
    allx2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In x direction
    ally1  :  list
        CCD column position of target’s flux-weighted centroid. In y direction
    ally2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In y direction
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID
  
    Returns
    -------
        Plot of the centroid positions in the x and y direction. 
        The black points are the centroid position - moves if this is a blend.

    '''


    gs = len(peak_list) # the grid size
    
    centroid1 = [allx1, ally1]  # CCD column position of target’s flux-weighted centroid.
    centroid2 = [allx2, ally2]  # The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects.
    
    if gs == 1:
        plt.figure(figsize=(7, 7))
    else: 
        plt.figure(figsize=(5 * gs, 7))

    for g,peak in enumerate(peak_list):
        
        for i,cen1 in enumerate(centroid1):
            
            mask = (np.array(alltime12) < peak+2) & (np.array(alltime12) > peak-2)
            
            # ------ ALL X --------
            plt.subplot(2,gs,g+1) # define the plotting area

            minf = np.nanmin(np.concatenate((np.array(centroid1[0])[mask],np.array(centroid2[0])[mask]), axis=0) )
            maxf = np.nanmax(np.concatenate((np.array(centroid1[0])[mask],np.array(centroid2[0])[mask]), axis=0) )
    
            #maxf = np.nanmax(np.array(centroid1[0])[mask])
            height = maxf - minf
            if i == 0:
                plt.plot(np.array(alltime12)[mask],np.array(centroid2[0])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='Local Motion')
                plt.plot(np.array(alltime12)[mask],np.array(centroid1[0])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='Flux-weighted Centroid')  
            else:
                plt.plot(np.array(alltime12)[mask],np.array(centroid2[0])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='_nolegend_')
                plt.plot(np.array(alltime12)[mask],np.array(centroid1[0])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='_nolegend_') 
            
            plt.legend()
            plt.tight_layout()
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = 'orange')
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.title('x centroid, Transit {}'.format(g+1))    
            
            # ------ ALL Y --------
            plt.subplot(2,gs,gs+1+g) # define the plotting area
            minf = np.nanmin(np.concatenate((np.array(centroid1[1])[mask],np.array(centroid2[1])[mask]), axis=0) )
            maxf = np.nanmax(np.concatenate((np.array(centroid1[1])[mask],np.array(centroid2[1])[mask]), axis=0) )
                    
            height = maxf - minf

            if i == 0:
                plt.plot(np.array(alltime12)[mask],np.array(centroid2[1])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='Local Motion')
                plt.plot(np.array(alltime12)[mask],np.array(centroid1[1])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='Flux-weighted Centroid')  
            else:
                plt.plot(np.array(alltime12)[mask],np.array(centroid2[1])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='_nolegend_')
                plt.plot(np.array(alltime12)[mask],np.array(centroid1[1])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='_nolegend_') 

            #plt.legend()
            plt.tight_layout()
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = 'orange')
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.title('y centroid, Transit {}'.format(g+1))        

    if save == True:
        plt.savefig('{}/{}/{}_centroids.png'.format(indir, tic, tic), format='png')


    #206361691
    
def plot_aperturesize(tic, indir,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, peak_list, save = False, show = False, FFI = False):
    '''               tic,  indir,      alltime,         alltime,         alltime,      allflux_normal, allflux_normal, allflux_small,   peak_list

    LC plot around the time of transit-event extracted in two different aperture sizes. The LC is not corrected for any systematics. 

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
    TESS_unbinned_t_l   :  list
        time for LC extracted with the TESS SPOC pipeline aperture 
    TESS_binned_t_l   :  list
        binned time for LC extracted with the TESS SPOC pipeline aperture
    small_binned_t_l   :  list
        binned time for LC extracted with the aperture that is 50% smaller than the TESS SPOC pipeline aperture 
    TESS_unbinned_l   :  list
        flux extracted with the TESS SPOC pipeline aperture 
    TESS_binned_l   :  list
        binned flux extracted with the TESS SPOC pipeline aperture
    small_binned_l   :  list
        binned flux extracted with the aperture that is 50% smaller than the TESS SPOC pipeline aperture
    peak_list   :  list
        list of the marked transit events
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID
  
    Returns
    -------
        plot of the transit-event extracted in two different aperture sizes. EBs will exhibit different transit shapes with the different aperture sizes. 
        The time of the transit event is indicated by the vertical line. 
    
    '''
    if FFI == True:
        frame_width = 1
    else:
        frame_width = 0.5


    gs = len(peak_list)

    if gs == 1:

        plt.figure(figsize=(7,4))
        plt.tight_layout()
        for g,peak in enumerate(peak_list):
        
            mask_unb = (np.array(TESS_unbinned_t_l) < peak+frame_width) & (np.array(TESS_unbinned_t_l) > peak-frame_width)
            mask_bin = (np.array(TESS_binned_t_l) < peak+frame_width) & (np.array(TESS_binned_t_l) > peak-frame_width)
            mask_small = (np.array(small_binned_t_l) < peak+frame_width) & (np.array(small_binned_t_l) > peak-frame_width)
        
            if FFI == False:
                minf = np.nanmin(np.array(TESS_unbinned_l)[mask_unb])
                maxf = np.nanmax(np.array(TESS_unbinned_l)[mask_unb])
            else:
                minf = np.nanmin(np.array(TESS_binned_l)[mask_bin]) 
                maxf = np.nanmax(np.array(TESS_binned_l)[mask_bin])
                diff = maxf - minf
                minf = minf - (diff * 0.01)
                maxf = maxf + (diff * 0.01)

            if FFI == False:
                plt.scatter(TESS_unbinned_t_l, TESS_unbinned_l, s = 3, marker = 's',alpha = 0.4, color = 'black', label = 'TESS unbinned')
            
            plt.scatter(TESS_binned_t_l, TESS_binned_l, s = 11,  marker = 'o', alpha = 1, color = 'blue', label = 'TESS ap')
            plt.scatter(small_binned_t_l, small_binned_l, s = 12, marker = '>', alpha =1, color = 'red', label = 'Small ap')
            
            plt.title("Detrended LC with various aperture sizes for TIC {}".format(tic), fontsize = 12)
            plt.tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
            plt.tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
            plt.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
            
            #plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
            plt.xlim(peak-frame_width, peak+frame_width)
            plt.axvline(peak, color = 'orange', linestyle = '--')
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.title('Aperture Size Test, Transit {}'.format(g+1), fontsize = 12)
            

        for transit in peak_list:
            plt.axvline(transit, color = 'orange', linestyle = '--', linewidth = 2)
            
        plt.legend(fontsize = 13)

        if save == True:
            plt.savefig('{}/{}/{}_aperture_size.png'.format(indir, tic, tic), format='png')

        if show == True:
            plt.show()
        else:
            plt.close()
    
    else:   


        plt.figure(figsize=(gs*6,4))

        for g,peak in enumerate(peak_list):
            

            mask_unb = (np.array(TESS_unbinned_t_l) < peak+frame_width) & (np.array(TESS_unbinned_t_l) > peak-frame_width)
            mask_bin = (np.array(TESS_binned_t_l) < peak+frame_width) & (np.array(TESS_binned_t_l) > peak-frame_width)
            mask_small = (np.array(small_binned_t_l) < peak+frame_width) & (np.array(small_binned_t_l) > peak-frame_width)
            
            if np.sum(mask_unb) != 0:
    
                if FFI == False:
                    minf = np.nanmin(np.array(TESS_unbinned_l)[mask_unb])
                    maxf = np.nanmax(np.array(TESS_unbinned_l)[mask_unb])
                else:
                    minf = np.nanmin(np.array(TESS_binned_l)[mask_bin]) 
                    maxf = np.nanmax(np.array(TESS_binned_l)[mask_bin])
                    diff = maxf - minf
                    minf = minf - (diff * 0.01)
                    maxf = maxf + (diff * 0.01)


                plt.subplot(1,gs,g+1)
                
                if FFI == False:
                    plt.scatter(TESS_unbinned_t_l, TESS_unbinned_l, s = 3, marker = 's',alpha = 0.4, color = 'black', label = 'TESS Aperture unbinned')
                plt.scatter(TESS_binned_t_l, TESS_binned_l, s = 11,  marker = 'o', alpha = 1, color = 'blue', label = 'TESS Aperture binned')
                plt.scatter(small_binned_t_l, small_binned_l, s = 12, marker = '>', alpha =1, color = 'red', label = 'Small Aperture binned')

                if g > 0:
                    plt.tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)
                    plt.yticks([])
    
                plt.tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
                plt.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
                
                plt.tight_layout()
                plt.xlim(peak-frame_width, peak+frame_width)
                plt.axvline(peak, color = 'darkorange', linestyle = '--')
                plt.ylim(minf,maxf)
                plt.xlabel('Time (BJD-2457000)')
                plt.title('Aperture Size Test, Transit {}'.format(g+1), fontsize = 12)


        if save == True:
            plt.savefig('{}/{}/{}_aperture_size.png'.format(indir, tic, tic), format='png')

        if show == True:
            plt.show()
        else:
            plt.close()

def plot_background(tic, indir,alltime, allfbkg, peak_list, save = False, show = False):
    
    '''
    LC of the bakcground flux at the time of the transit event. 

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
    alltime  :  list
        times (not binned)
    allfbkg  :  list
        background flux
    peak_list   :  list
        list of the marked transit events
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID
  
    Returns
    -------
        plot of the background flux at the time of the marked transit-events.
        The backrgound flux should not exhibit any 'spikes' or 'unusual' behaviour at the time of the transit event. 
        The time of the transit event is indicated by the vertical line. 

    '''
      

    gs = len(peak_list) # the grid size

    if len(peak_list) == 1:
        plt.figure(figsize=(7,4))
        plt.tight_layout()

        for g,peak in enumerate(peak_list):
        
            mask_dd = (np.array(alltime) < peak+2) & (np.array(alltime) > peak-2)
            
            minf = np.nanmin(np.array(allfbkg)[mask_dd])
            maxf = np.nanmax(np.array(allfbkg)[mask_dd])
            height = maxf - minf
            
            plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = 'orange')
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.ylabel('Flux')
            plt.title('Background Flux, Transit {}'.format(g+1))
        
        if save == True:
            plt.savefig('{}/{}/{}_background.png'.format(indir, tic, tic), format='png')
        
        if show == True:
            plt.show()
        else:
            plt.close()

    else:

        gs = len(peak_list) # the grid size
        plt.figure(figsize=(6 * gs, 4))

        for g,peak in enumerate(peak_list):
    
            mask_dd = (np.array(alltime) < peak+2) & (np.array(alltime) > peak-2)
    
            if np.sum(mask_dd) != 0:
    
                minf = np.nanmin(np.array(allfbkg)[mask_dd])
                maxf = np.nanmax(np.array(allfbkg)[mask_dd])
                height = maxf - minf
    
                plt.subplot(1,gs,g+1)
                plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
                plt.xlim(peak-2, peak+2)
                plt.axvline(peak, color = 'darkorange')
                
                try:
                    plt.ylim(minf,maxf)
                except:
                    print('axis limit error (centroid positions)')
    
                plt.tight_layout()
                plt.xlabel('Time (BJD-2457000)')
                plt.ylabel('Flux')
                plt.title('Background Flux, Transit {}'.format(g+1))

        if save == True:
            plt.savefig('{}/{}/{}_background.png'.format(indir, tic, tic), format='png')

        if show == True:
            plt.show()
        else:
            plt.close()


def plot_TESS_stars(tic,indir,peak_list, peak_sec, tpf_list, save = False, show = False):
    
    '''
    Plot of the field fluxa round the target star showing nearby stars that are brighter than magnitude 15.
    
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
    peak_list   :  list
        list of the marked transit events
    peak_sec  :  list or str
        list of the sectors that are being analyse.
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID
  
    Returns
    -------
        Plot of the averaged flux per pixel around the target. The red star indicated the location of the target.
        The orange circles show the location of nearby stars with magnitudes brighter than 15 mag.

    '''


    # always check whether the file already exists... as to not waste computer power and time
    sector = str(peak_sec[0])

    starName = "TIC " + str(tic)
    radSearch = 4/60 #radius in degrees

    catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
    
    # ra and dec of the target star
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']


    # ----------
    survey = 'DSS2 Red'
    fig, ax = plt.subplots()
    plt.axis("off")
    
    args2 = {}
    args2.setdefault('linewidth', 1)
    args2.setdefault('color', 'red')
    
    # get the SDSS image and save it - this will appear in the report
    target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    target = FixedTarget(coord=target_coord, name="Survey = {}".format(survey))
    
    ax, hdu = plot_finder_image(target, survey = survey, reticle='True', reticle_style_kwargs = args2)
    
    plt.savefig('{}/{}/{}_SDSSstar_field.png'.format(indir, tic, tic), format='png', bbox_inches = 'tight', dpi = 100)
    plt.close()

    # ----------

    # Create a list of nearby bright stars (tess magnitude less than 14) from the rest of the data for later.
    bright = catalogData['Tmag'] < 17

    start = [np.float64(peak_list[0]) - 0.2]
    end = [np.float64(peak_list[0]) + 0.2]

    #background is just the array of the flux of all the pixels (used for backrgoudn in pixel level LC plot so mildly confusing and should change)
    for i, tpf in enumerate(tpf_list):

        # plt.subplot(row column number)

        if (start > np.nanmin(tpf.time) and start < np.nanmax(tpf.time)):
            fig, ax = plt.subplots(figsize=(5,5))
            plt.tight_layout()

            sector =  tpf.header['SECTOR']

            ax = plt.subplot(projection=tpf.wcs)

            plot_cutout(tpf.flux[0])

            ra_stars, dec_stars = catalogData[bright]['ra'], catalogData[bright]['dec']
            s = np.maximum((19 - catalogData[bright]['Tmag'])*5, 0)
            plt.scatter(ra_stars, dec_stars, s=s, transform=ax.get_transform('icrs'), color='orange', zorder=100)

            # plot the target
            plt.scatter(ra, dec, s= 200, transform=ax.get_transform('icrs'), marker = '*', color='red', zorder=100)

            plt.title("Sector {}".format(sector), fontsize = 12)

            plt.xlim(-0.5,10.5)
            plt.ylim(-0.5,10.5)

            if save == True:
                plt.savefig('{}/{}/{}_star_field.png'.format(indir, tic, tic), format='png')
    
            if show == True:
                plt.show()
            else:
                plt.close()
    

    return catalogData['Tmag'][0], catalogData['Teff'][0], catalogData['rad'][0], catalogData['mass'][0]


def plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list,t_list, peak_list, save = False, show = False, FFI = False):
    
    '''
    Plot the LC for each pixel around the time of the transit like event. 

    Parameters
    ----------
    indir : str
        path to where the files will be saved.
    peak_sec  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    tic : str
        TIC (Tess Input Catalog) ID of the target
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
    peak_list  :  int
        list of all the marked transits
    save (default = False)
        if save = True, the figure is saved in the directory which has the name of the TIC ID
  410153553
    Returns
    -------
        Plot of the normalised LC for each pixel around the time of the transit like event. 
        The pixel backrgound colour represents the average flux. 
        The time of the transit is highlighted in red/gold for each pixel LC.
    '''

    #only plot this for the first transit -- alter here to plot for all or others.

    for idx, X1 in enumerate(X1_list):

        mapimg = apmask_list[idx]
        X4 = X4_list[idx]
        oot = oot_list[idx]
        #intr = intr_list[n]
        bkg = bkg_list[idx]
        apmask = apmask_list[idx]
        arrshape = arrshape_list[idx]
        t = t_list[idx]
        peak = peak_list[idx]

        ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
        hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

        fig, ax = plt.subplots(arrshape[1], arrshape[2], sharex = True, sharey = False, gridspec_kw={'hspace': 0 ,'wspace': 0}, figsize=(8,8))
        
        plt.tight_layout()

        try:
            color = plt.cm.viridis(np.linspace(0, 1,int(np.nanmax(bkg))-int(np.nanmin(bkg))+1))
            simplebkg = False
        except:
            simplebkg = True

        for i in range(0,arrshape[1]):
            print ("{}   out of    {} ".format(i+1,arrshape[1] ))
            for j in range(0,arrshape[2]):
    
    
                apmask = np.zeros(arrshape[1:], dtype=np.int)
                apmask[i,j] = 1
                apmask = apmask.astype(bool)
                
                flux = X1[:,apmask.flatten()].sum(axis=1)
    
                m = np.nanmedian(flux[oot])
                
                normalizedflux = flux/m

                # bin the data 
                f1 = normalizedflux
                time = t

                if FFI == False:
                    binfac = 5
    
                    N       = len(time)
                    n       = int(np.floor(N/binfac)*binfac)
                    X       = np.zeros((2,n))
                    X[0,:]  = time[:n]
                    X[1,:]  = f1[:n]
                    Xb      = rebin(X, (2,int(n/binfac)))
        
                    # binned data
                    time_binned    =    np.array(Xb[0])
                    flux_binned  =   np.array(Xb[1])

                else:
                    # binned data -
                    time_binned    =    np.array(time)
                    flux_binned  =   np.array(flux)
    

                # create a mask that only looks at the times cut around the transit-event
                timemask = (time_binned < peak+0.7) & (time_binned > peak-0.7)
                
                #timemask = 
                intr = abs(peak-time_binned) < 0.1

                if simplebkg == True:
                    ax[i, j].set_facecolor(color = 'k')
                    linecolor = 'w'
                    transitcolor = 'gold'                   
                else:
                    ax[i, j].set_facecolor(color = color[int(bkg[i,j])-int(np.nanmin(bkg))])

                    if int(bkg[i,j])-abs(int(np.nanmin(bkg))) > ((np.nanmax(bkg))-abs(int(np.nanmin(bkg))))/2:
                        linecolor = 'k'
                        transitcolor = 'orangered'
                    else:
                        linecolor = 'w'
                        transitcolor = 'gold'
                

                ax[i, j].plot(time_binned[timemask],flux_binned[timemask], color = linecolor, marker = '.', markersize=1, lw = 0) 
                ax[i, j].plot(time_binned[intr],flux_binned[intr], color = transitcolor, marker = '.', markersize=1, lw = 0) 
                
                ax[i,j].set_yticklabels([])
                ax[i,j].set_xticklabels([])
    
        # ------------------    
        
        print ("\n Calculating the Aperture Mask...", end =" ")
        
        for i in range(0,len(ver_seg[1])):
            ax[ver_seg[0][i], ver_seg[1][i]].spines['right'].set_color('red')
            ax[ver_seg[0][i], ver_seg[1][i]].spines['right'].set_linewidth(10)
        
        for j in range(0,len(hor_seg[1])):
            ax[hor_seg[0][j], hor_seg[1][j]].spines['bottom'].set_color('red')
            ax[hor_seg[0][j], hor_seg[1][j]].spines['bottom'].set_linewidth(10)
        print ("Done.\n")
        # ------------------
        
        print ("Waiting on plot...")
        plt.xlim(peak-0.7,peak+0.7)

        if save == True:
            plt.savefig('{}/{}/{}_individual_pixel_LCs_{}.png'.format(indir, tic,tic, idx), format='png')

        if show == True:
            plt.show()
        else:
            plt.close()

# full light curve with the momentum dumps
def plot_full_md(tic, indir, alltime, allflux,allline,alltimebinned,allfluxbinned, peak_list, show = False, FFI = False):
    '''
    
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir : str
        path to where the files will be saved.
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
    peak_list  :  list
        list of all the marked transit-events

    Returns
    -------
        Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)). 
        The plot is always saved as this funtion is only called if the 'save' option was chosen
    '''
    gs = len(peak_list)

    plt.figure(figsize=(15,10))
    plt.tight_layout()

    if FFI == False:
        # rename things so that the code didn't have to be changed - not a very 'tidy' solution.
        time_dd = alltime
        flux_dd = allflux
        time_dd_binned = alltimebinned
        flux_dd_binned = allfluxbinned

    else:
        # --- do some sigma clipping to make the LC look better ---- 
        MAD = median_absolute_deviation(allflux)
        madrange = (5 * MAD * 1.456)
        ymask = (allflux < 1 + madrange) * (allflux > 1 - madrange) 

        time_dd = np.array(alltime)[ymask]
        flux_dd = np.array(allflux)[ymask]
        time_dd_binned = np.array(alltime)[ymask]
        flux_dd_binned = np.array(allflux)[ymask]
        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- 

    line_dd = allline

    minf_list = []
    maxf_list = []

    for peak in peak_list:
        mask_dd = (np.array(time_dd) < peak+0.5) & (np.array(time_dd) > peak-0.5)
        minf0 = np.nanmin(np.array(flux_dd)[mask_dd])
        maxf0 = np.nanmax(np.array(flux_dd)[mask_dd])
        minf_list.append(minf0)
        maxf_list.append(maxf0)

    minf = np.nanmin(minf_list)
    maxf = np.nanmin(maxf_list)
    height = maxf - minf

    for g,peak in enumerate(peak_list):

        mask_dd = (np.array(time_dd) < peak+0.5) & (np.array(time_dd) > peak-0.5)
        mask_dd_binned = (np.array(time_dd_binned) < peak+0.5) & (np.array(time_dd_binned) > peak-0.5)

        if np.sum(mask_dd) != 0:

            if gs == 1:
                plt.subplot(2,4,(6,7))

            else:
                plt.subplot(2,gs,(gs+g+1))

            plt.plot(np.array(time_dd)[mask_dd], np.array(flux_dd)[mask_dd], 'o', markersize = 4, color = 'orange', alpha = 0.8, label = "unbinned", markerfacecolor='white', zorder=1)
            plt.plot(np.array(time_dd_binned)[mask_dd_binned], np.array(flux_dd_binned)[mask_dd_binned], marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 5, label = 'binning = 7', markerfacecolor='k', zorder=2)

            plt.vlines(line_dd, minf,minf + height*0.25 , colors = 'r', label = "Momentum Dump", zorder=3)
            plt.vlines([peak], minf,minf + height*0.25 , linewidth = 3, colors = 'k', linestyle = '--', zorder=4)

            plt.xlim(peak-0.5, peak+0.5)
            #plt.axvline(peak, color = 'k')

            try:
                plt.ylim(minf,maxf)
            except:
                print ('axis limits error (momentun dumps)')

            plt.xlabel('BJD-2457000')
            plt.title('Transit {}'.format(g+1))

        else:
            print ("size 0 so skip")


    if gs == 1:
        plt.subplot(2,4,(1,4))
    else:
        plt.subplot(2,gs,(1,gs))

    if FFI == False:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = 'orange', alpha = 0.9, label = "unbinned", markerfacecolor='white', zorder=1)
        plt.plot(np.array(time_dd_binned), np.array(flux_dd_binned), marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 1, label = 'binning = 7', markerfacecolor='k', zorder=2)
    else:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = '#054950', alpha = 0.9, label = "unbinned", markerfacecolor='#054950', zorder=1)

    minf_full = np.nanmin(np.array(flux_dd))
    maxf_full = np.nanmax(np.array(flux_dd))

    height_full = maxf - minf

    plt.vlines(line_dd, minf_full,minf_full + height_full*0.3, colors = 'r', label = "Momentum Dump", zorder=3)
    plt.ylim(minf_full, maxf_full)
    plt.xlim(np.nanmin(np.array(time_dd)), np.nanmax(np.array(time_dd)))
    plt.xlabel('BJD-2457000')

    plt.vlines(peak_list, minf_full,minf_full + height*0.3 , colors = 'k', linestyle = '--', linewidth = 3, zorder=4)

    plt.savefig("{}/{}/{}_fullLC_md.png".format(indir,tic,tic), format='png', bbox_inches = 'tight')
    
    plt.close()


def plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model, results,period,duration,t0, in_transit = [0], save = False, show = False):

    if len(in_transit) == 1:  # conditions for the first 'round' of plotting

        color1 = '#DC143C'
        color2 = 'orange'
        title = 'Initial BLS'
        name = '{}_bls_first.png'.format(tic)

    else:  # conditions for the second 'round' of plotting once the first event has been removed
        color1 = 'deepskyblue'
        color2 = '#4682B4'
        title = 'Initial event removed'
        name = '{}_bls_second.png'.format(tic)
        
    fig, axes = plt.subplots(3, 1, figsize=(5, 7))
    
    # Highlight the harmonics of the peak period
    ax = axes[0]
    ax.axvline(period, alpha=0.4, lw=5, color = color1)
    for n in range(2, 15):
        ax.axvline(n*period, alpha=0.4, lw=2, linestyle="dashed", color = color2)
        ax.axvline(period / n, alpha=0.4, lw=2, linestyle="dashed", color = color2)
    
    # Plot the periodogram
    ax.plot(results.period, results.power, "k", lw=0.5, label = 'P = %.3f T0 = %.3f' % (period,t0))
    
    ax.set_title(title)
    ax.set_xlim(results.period.min(), results.period.max())
    ax.set_xlabel("period (days)")
    ax.set_ylabel("log likelihood")
    ax.legend(fontsize = 10, loc = 1)
    
    # Plot the light curve and best-fit model
    ax = axes[1]
    
    if len(in_transit) == 1:
        ax.plot(alltime, allflux, marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(alltimebinned, allfluxbinned, marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')
    else:
        ax.plot(alltime, allflux, marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(alltimebinned[~in_transit], allfluxbinned[~in_transit], marker ="o", alpha = 0.6, color = 'black',  MarkerFaceColor = 'none', ms=3, lw = 0)

    x = np.linspace(alltimebinned.min(), alltimebinned.max(), 3*len(alltimebinned))
    f = model.model(x, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(alltimebinned.min(), alltimebinned.max())
    ax.set_xlabel("time (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    
    ax = axes[2]
    if len(in_transit) == 1: 
        x_binned = (alltimebinned - t0 + 0.5*period) % period - 0.5*period
        x = (alltime - t0 + 0.5*period) % period - 0.5*period
    else:
        x_binned = (alltimebinned[~in_transit] - t0 + 0.5*period) % period - 0.5*period
        x = (alltime - t0 + 0.5*period) % period - 0.5*period
    
    m_binned = np.abs(x_binned) < 0.5 
    m = np.abs(x) < 0.5 
    
    if len(in_transit) == 1: 
        ax.plot(x[m], allflux[m],marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(x_binned[m_binned], allfluxbinned[m_binned], marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')
        
    else:
        ax.plot(x[m], allflux[m],marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(x_binned[m_binned], allfluxbinned[~in_transit][m_binned], marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')
        
    x = np.linspace(-0.5, 0.5, 1000)
    f = model.model(x + t0, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("time since transit (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    plt.tight_layout()

    if save == True:
        plt.savefig('{}/{}/{}'.format(indir, tic, name), format='png')
    
    if show == True:
        plt.show()
    else:
        plt.close()


def plot_bls_FFI(tic, indir, alltime, allflux, model, results,period,duration,t0, in_transit = [0], save = False, show = False):

    if len(in_transit) == 1:  # conditions for the first 'round' of plotting

        color1 = 'navy'
        color2 = 'orangered'

        title = 'Initial BLS'
        name = '{}_bls_first.png'.format(tic)

    else:  # conditions for the second 'round' of plotting once the first event has been removed
        color1 = 'deepskyblue'
        color2 = '#4682B4'
        title = 'Initial event removed'
        name = '{}_bls_second.png'.format(tic)
        
    fig, axes = plt.subplots(3, 1, figsize=(5, 7))
    
    # Highlight the harmonics of the peak period
    ax = axes[0]
    ax.axvline(period, alpha=0.4, lw=5, color = color1)
    for n in range(2, 15):
        ax.axvline(n*period, alpha=0.4, lw=2, linestyle="dashed", color = color2)
        ax.axvline(period / n, alpha=0.4, lw=2, linestyle="dashed", color = color2)
    
    # Plot the periodogram
    ax.plot(results.period, results.power, "k", lw=0.5, label = 'P = %.3f T0 = %.3f' % (period,t0))
    
    ax.set_title(title)
    ax.set_xlim(results.period.min(), results.period.max())
    ax.set_xlabel("period (days)")
    ax.set_ylabel("log likelihood")
    ax.legend(fontsize = 10, loc = 1)
    
    # Plot the light curve and best-fit model
    ax = axes[1]
    
    if len(in_transit) == 1:
        ax.plot(alltime, allflux, marker =".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
    else:
        ax.plot(alltime, allflux, marker =".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')

    x = np.linspace(alltime.min(), alltime.max(), 3*len(alltime))
    f = model.model(x, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(alltime.min(), alltime.max())
    ax.set_xlabel("time (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    
    ax = axes[2]
    if len(in_transit) == 1: 
        x = (alltime - t0 + 0.5*period) % period - 0.5*period
    else:
        x = (alltime - t0 + 0.5*period) % period - 0.5*period
    
    m = np.abs(x) < 0.5 
    
    if len(in_transit) == 1: 
        ax.plot(x[m], allflux[m],marker =".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')

    else:
        ax.plot(x[m], allflux[m],marker =".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')


    x = np.linspace(-0.5, 0.5, 1000)
    f = model.model(x + t0, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("time since transit (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    plt.tight_layout()

    if save == True:
        plt.savefig('{}/{}/{}'.format(indir, tic, name), format='png')
    
    if show == True:
        plt.show()
    else:
        plt.close()


def plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, save = False, show = False):
    

    plt.figure(figsize=(16,3.5*len(T0_list)))

    count = 0

    for idx, X4 in enumerate(X4_list):
        
        oot = oot_list[idx]
        intr = intr_list[idx]
        T0 = T0_list[idx]
        t = t_list[idx]
        tpf_filt  =  tpf_filt_list[idx]
        
        intr = abs(T0-t) < 0.25
        oot = (abs(T0-t) < 0.5) * (abs(T0-t) < 0.3)
        img_intr = tpf_filt[intr,:,:].sum(axis=0)/float(intr.sum())
        img_oot = tpf_filt[oot,:,:].sum(axis=0)/float(oot.sum())
        img_diff = img_oot-img_intr
        
        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_intr)
        plt.colorbar()
        plt.title("t = {} days \n In Transit Flux (e-/candence)".format(T0))

        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_oot)
        plt.colorbar()
        plt.title("Out of Transit Flux (e-/candence)")

        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_diff)
        plt.colorbar()
        plt.title("Difference Flux (e-/candence)")


    plt.tight_layout()

    if save == True:
        plt.savefig('{}/{}/{}_flux_comparison.png'.format(indir, tic, tic), format='png')
    
    if show == True:
        plt.show()
    else:
        plt.close()

# -----------------------
# other functions
# -----------------------
def rebin(arr,new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
        new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

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

# -----------------------
# only used for the Jupyter Notebook version
def transit_finder(transit, alltime, allline, allflux,alltimebinned, allfluxbinned, start_sec, end_sec, flux_min = None, flux_max = None):
    '''
    # only used for the Jupyter Notebook version
    '''
    # WHOLE
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(alltime, allflux, marker='o',lw = 0, markersize = 1, color = 'orange', alpha = 0.5, label = 'unbinned', markerfacecolor='white')
    
    ax.plot(alltimebinned, allfluxbinned, marker='o',color = 'k', alpha = 0.6, lw = 0, markersize = 1, label = 'binning = 7', markerfacecolor='k')
        
    #ax.set_ylim(np.nanmin(allflux),np.nanmax(allflux)-0.002)
        
    for p, peak in enumerate(start_sec):
        ax.axvline(end_sec[p], color = 'k', linewidth = 0.2, alpha =1, label= 'Sector limits')

    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf

 
    ax.tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax.tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax.set_xlabel("Time (BJD-2457000)", fontsize = 12)
    ax.set_ylabel("Normalised Flux", fontsize = 12)
    
    ax.axvline(transit-1, color = 'b', linewidth = 2)
    ax.axvline(transit+1, color = 'b', linewidth = 2)
    
    ax.axvline(transit, color = 'r', linewidth = 2)
    
    minorLocator = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(direction='in', which ='minor', colors='grey',length=3, labelsize=13)
    ax.vlines(allline, minf,minf + height*0.3 , colors = 'r', label = "Momentum Dump")

    minorLocator = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)
    
    if flux_min != None:
        ax.set_ylim(flux_min,flux_max)

    plt.show()
    
    # Zoomed in 

    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(alltime, allflux, marker='o',lw = 0, markersize = 4, color = 'orange', alpha = 0.8, label = 'unbinned', markerfacecolor='white')
    
    ax.plot(alltimebinned, allfluxbinned, marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, label = 'binning = 7', markerfacecolor='k')
    
    ax.vlines(allline, minf,minf + height*0.3 , colors = 'r', label = "Momentum Dump")

    if flux_min != None:
        ax.set_ylim(flux_min,flux_max)
    
    ax.axvline(transit, color = 'r', linewidth = 2)
        
    for p, peak in enumerate(start_sec):
        ax.axvline(end_sec[p], color = 'k', linewidth = 0.2, alpha =1, label= 'Sector limits')
    
    ax.tick_params(axis="y",direction="inout", labelsize = 12) #, pad= -20)
    ax.tick_params(axis="x",direction="inout", labelsize = 12) #, pad= -17)   
    ax.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax.set_xlabel("Time (BJD-2457000)", fontsize = 12)
    ax.set_ylabel("Normalised Flux", fontsize = 12)
    
    ax.set_xlim(transit-1, transit+1)
    mask = (np.array(alltimebinned) > transit-1) & (np.array(alltimebinned) < transit+1)
 

    minorLocator = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(direction='in', which ='minor', colors='grey',length=3, labelsize=13)
    
    minorLocator = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)
    
    print ("TRANSIT: {}".format(transit))
    plt.show()



