from __future__ import print_function, absolute_import, division

import os
import csv
import sys
import json
import astropy
import requests
import itertools
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sb
from glob import glob
import lightkurve as lk
from os.path import exists
import matplotlib.pyplot as plt

from astropy.wcs import WCS
import astropy.io.fits as pf
#from astropy.io.fits import pf.Undefined, pf.BinTableHDU
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.stats import median_absolute_deviation, BoxLeastSquares

from scipy.ndimage import label
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import find_optimal_celestial_wcs

from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons

import http.client as httplib
import mpl_toolkits.axes_grid1
from astroquery.mast import Catalogs
from sklearn.decomposition import PCA
from urllib.parse import quote as urlencode
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from astropy.utils.data import clear_download_cache
clear_download_cache()

# custom modules
from LATTE import filters
from LATTE import LATTEbrew as brew

main_pot_color = '#ffa31a'
main_plot_fontsize = 13
md_color = '#143085'
cutout_window = 0.7
individual_ylim = False

sb.set(style='ticks')
sb.set_color_codes('deep')

'''
Overview of LATTE scipts:

__main__.py      : Intitialises the parameters, what TIC ID, sector, checks for downloads of the data, FFI or not?

LATTEutils.py    : All the functions needed to download data and text files, runs the interactive gui, all of the plotting and data handling.

LATTE_DV.py      : Scipt to combine all of the results from LATTEbrew in order to generate a pdf data validation report.

LATTEbrew.py     : Calls the LATTEutils functions in turn in order to store the results and keeps track of what has been generated. Calls the LATTE_DV.py function to collate all the results.


The functions in this scipt are organised as follows:

1) 3 functions that call the interative GUI plotting for the short cadence and FFI data. There are two functions for the FFI, one that asks the user to
identify the aperture sizes manually and one that does it automatically. These functions are called from __main__.py and call the 'brew' scipt at the end.
The brew scipt then continues with collating the data by calling the remaining functions in this script.

2) 4 functions that download the text files needed for the LATTE to run. These are only exectuted the first time that LATTE is run, and when new TESS data is made available.
These functions are called from within __main__.py

3) 3 functions that determine the sectors in which the target was observed and TICs of nearest neighbours.

4) 7 functions needed to download and process the data, both from the FFIs and short cadence data.
These functions are called from within LATTEbrew.py and from withint the interact functions in this script.

5) 14 functions for plotting the results, both from the FFIs and short cadence data.
These functions are called from within LATTEbrew.py.

4) 4 short functions that are needed in various other functions within this scipt, such as one to bin the data.

'''

# check whether pyaneti has been sucessfully installed - if it has not been installed, don't give the option to model the data
# NOTE: in this version the Pyaneti modeling will not work (except on Nora's computer) - this is being handled in the next release.

try:
    from LATTE import pyaneti_LATTE
    pyaneti_installed = True
except:
    pyaneti_installed = False

# --------------------------------------------
# -----------------interact ------------------
# --------------------------------------------

# the main interactive tool used to identify the times of the transit-like events
# This function is called from __main__.py
def interact_LATTE(tic, indir, syspath, sectors_all, sectors, ra, dec, args):

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

    # function needed in the rebinning - changes the shape of the array
    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)

    # ---------------
    # this needs to be global for the sector slider
    global in_sec

    # call function to download the data from MAST
    print ("Start lightcurve data download.....", end =" ")
    alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir, sectors, tic)
    print ("done.\n")

    plt.close('all') # make sure that all figures are close.

    # ---------------

    # -------------------------
    # Plot the interactive plot - uses matplolib
    # -------------------------

    fig, ax = plt.subplots(2, 1, figsize=(11,7.5))
    plt.tight_layout()

    # Adjust the plots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.3)

    fluxmin = np.nanmin(allflux)
    fluxmax = np.nanmax(allflux)
    flux_diff = fluxmax - fluxmin
    # function to define the plotting area around the transit event.
    # this needs to be in a function as the area changes with the interactive slider.
    def cutout(transit):
        mask_binned = (np.array(alltimebinned) > transit-1) & (np.array(alltimebinned) < transit+1)
        mask = (np.array(alltime) > transit-1) & (np.array(alltime) < transit+1)

        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux), np.array(alltimebinned)[mask_binned], np.array(allfluxbinned)[mask_binned], np.array(alltimebinned), np.array(allfluxbinned)]

    # ---------------

    # function to define the bin factor
    # this needs to be in a function as binning can be changed with the interactive buttons.
    def binning(binfac):
        N      = len(alltime)
        n      = int(np.floor(N/binfac)*binfac)
        X      = np.zeros((2,n))
        X[0,:]  = alltime[:n]
        X[1,:]  = allflux[:n]
        Xb    = rebin(X, (2,int(n/binfac)))

        time_binned = Xb[0]
        flux_binned = Xb[1]

        return [time_binned, flux_binned]
    # ---------------

    # define the slider that let's you chose the sector to look at
    class PageSlider(matplotlib.widgets.Slider):

        def __init__(self, ax, label, numpages = 10, valinit=0, valfmt='%1d',
                     closedmin=True, closedmax=True,
                     dragging=True, **kwargs):

            self.facecolor=kwargs.get('facecolor',"w")
            self.activecolor = kwargs.pop('activecolor',"b")
            self.fontsize = kwargs.pop('fontsize', 10)
            self.numpages = numpages

            super(PageSlider, self).__init__(ax, label, 0, numpages,
                                valinit=valinit, valfmt=valfmt, **kwargs)

            # make it so that numpages is the index of the array of the
            global in_sec

            self.poly.set_visible(False)
            self.vline.set_visible(False)
            self.pageRects = []
            for i in range(numpages):

                facecolor = self.activecolor if i==valinit else self.facecolor
                r  = Rectangle((float(i)/numpages, 0), 1./numpages, 1,
                                    transform=ax.transAxes, facecolor=facecolor)
                ax.add_artist(r)
                self.pageRects.append(r)

                seclist_names = ['All'] + in_sec

                sec_i = seclist_names[i]

                ax.text(float(i)/numpages+0.5/numpages, 0.5, str(sec_i),
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=self.fontsize)
            self.valtext.set_visible(False)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            bax = divider.append_axes("right", size="5%", pad=0.05)
            fax = divider.append_axes("right", size="5%", pad=0.05)
            self.button_back = matplotlib.widgets.Button(bax, label=r'$\blacktriangleleft$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_forward = matplotlib.widgets.Button(fax, label=r'$\blacktriangleright$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_back.label.set_fontsize(self.fontsize)
            self.button_forward.label.set_fontsize(self.fontsize)
            self.button_back.on_clicked(self.backward)
            self.button_forward.on_clicked(self.forward)

        def _update(self, event):
            super(PageSlider, self)._update(event)
            i = int(self.val)
            if i >=self.valmax:
                return
            self._colorize(i)

        def _colorize(self, i):
            for j in range(self.numpages):
                self.pageRects[j].set_facecolor(self.facecolor)
            self.pageRects[i].set_facecolor(self.activecolor)

        def forward(self, event):
            current_i = int(self.val)
            i = current_i+1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

        def backward(self, event):
            current_i = int(self.val)
            i = current_i-1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

    # ---------------


    # define the initial conditions - this is what appears on the plot initially but can change later.
    transit = np.nanmean(alltimebinned)  # start the cut-out in the centre of the plot.
    binfac = 7

    # FIRST PLOT - FULL LC
    # plot the unbinned flux
    [line_full] = ax[0].plot(alltime, allflux , marker='o',lw = 0, markersize = 4, color = main_pot_color, alpha = 0.8,  markerfacecolor='white')
    # plot the binned flux - it is called through the above function so that the binning can be changed.
    [line_full_binned] = ax[0].plot(binning(binfac)[0], binning(binfac)[1],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3,  markerfacecolor='k')
    # ---------------

    # SECOND PLOT - CUT OUT LC
    # plot the cut out around the time of the transit - called with function as it can be changed with the slider.
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',lw = 0, markersize = 4, color = main_pot_color, alpha = 0.8, markerfacecolor='white')
    [line_binned] =  ax[1].plot(cutout(transit)[4], cutout(transit)[5],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='k')
    # ---------------

    global transit_slider_ax
    global transit_slider

    if len(in_sec) == 1:

        # Define two sliders to change the y axis scale
        scale_slider_upper_ax  = fig.add_axes([0.24, 0.19, 0.33, 0.03], facecolor='silver')
        scale_slider_upper = Slider(scale_slider_upper_ax, 'Flux limits', 1, fluxmax + (flux_diff*0.1), valinit=fluxmax, color='white')
        scale_slider_upper.valtext.set_visible(False)

        scale_slider_lower_ax  = fig.add_axes([0.57, 0.19, 0.33, 0.03], facecolor='white') # x, y, width, height
        scale_slider_lower = Slider(scale_slider_lower_ax, ' ', fluxmin - (flux_diff*0.1), 1,valinit=fluxmin, color='silver')
        scale_slider_lower.valtext.set_visible(False)

        # Define the slider to change the transit-event time (and cut out region)
        transit_slider_ax  = fig.add_axes([0.24, 0.15, 0.66, 0.03])  # location of the slider
        transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=transit, color='teal')

    # if there it's in more than one sector there is an extra slide bar to chose to look at individual sectors.
    else:
        # Define the slider to change the y axis scale
        scale_slider_upper_ax  = fig.add_axes([0.24, 0.2, 0.33, 0.02], facecolor='silver')
        scale_slider_upper = Slider(scale_slider_upper_ax, 'Flux limits', 1, fluxmax + (flux_diff*0.1), valinit=fluxmax, color='white')
        scale_slider_upper.valtext.set_visible(False)

        scale_slider_lower_ax  = fig.add_axes([0.57, 0.2, 0.33, 0.02], facecolor='white')
        scale_slider_lower = Slider(scale_slider_lower_ax, ' ', fluxmin - (flux_diff*0.1), 1,valinit=fluxmin, color='silver')
        scale_slider_lower.valtext.set_visible(False)

        # Define the slider to change the transit-event time (and cut out region)
        transit_slider_ax  = fig.add_axes([0.24, 0.17, 0.66, 0.02])  # location of the slider
        transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=transit, color='teal')

        # define the slider to change the plotting region
        sector_slider_ax  = fig.add_axes([0.24, 0.11, 0.66, 0.03])
        sector_slider = PageSlider(sector_slider_ax, 'Sector', len(in_sec) + 1, activecolor=main_pot_color)

    # ---------------

    # deifne the intial x and y axis limits - y limit can be changed with slider, x limit cannot.
    ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])

    ax[1].set_xlim([np.nanmean(alltime)-1, np.nanmean(alltime)+1])
    ax[1].set_ylim([fluxmin, fluxmax])

    # ---------------
    # Define an action for acessing the required cut-out data and drawing the line when the slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])

        line_binned.set_xdata(cutout(transit_slider.val)[4])
        line_binned.set_ydata(cutout(transit_slider.val)[5])

        fig.canvas.draw_idle()  # update figure

    # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)

    # ---------------
    # Define an action for modifying the plot region on the second plot
    def update_axis(val):
        # ---------------
        # only ever plot one line so get rid of the old ones before plotting the new ones...
        if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

        # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
        lver0 = ax[0].axvline(transit_slider.val, color = 'r', linewidth = 2)
        lver1 = ax[1].axvline(transit_slider.val, color = 'r', linewidth = 2)

        ax[1].set_xlim([transit_slider.val - 1,transit_slider.val + 1])
        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)

    def update_plotting_region(val):
        global in_sec
        seclist_names = ['All'] + in_sec

        index = int(val)

        if index == 0:
            ax[0].set_title("TIC {}     Sectors {}".format(tic,str(in_sec)[1:-1]))
            ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
        else:
            ax[0].set_title("TIC {}   Sector {}".format(tic,seclist_names[index]))
            ax[0].set_xlim([start_sec[index-1][0],end_sec[index-1][0]])


    def update_yaxis_upper(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])

    def update_yaxis_lower(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])


    # ------ ON CLICK -------
    # define the action that lets you click on the image in order to chose the location instead of using the slider - probably more useful
    def onclick(event):

        # just need to know the x (time) position of the click. The y (flux) position doesn't matter.
        val = event.xdata

        # check that the click was within the plotting region. If not ignore the click
        if (event.inaxes == ax[0]) or (event.inaxes == ax[1]):

            # SECOND PLOT - CUT OUT LC
            # plot the cut out around the time of the transit - changed by clicked on the image (either one of them)
            line.set_xdata(cutout(val)[0])
            line.set_ydata(cutout(val)[1])

            line_binned.set_xdata(cutout(val)[4])
            line_binned.set_ydata(cutout(val)[5])

            # ---------------
            # only ever plot one line so get rid of the old ones before plotting the new ones...
            if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

            # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
            lver0 = ax[0].axvline(val, color = 'r', linewidth = 2)
            lver1 = ax[1].axvline(val, color = 'r', linewidth = 2)

            # if clicked on the image to select the time to zoom in on then update the plot in this function.
            ax[1].set_xlim([val - 1,val + 1])
            lver0.set_xdata(val)
            lver1.set_xdata(val)

            # also update the location of the slider (aesthetic reasons only)
            # Define the slider to change the transit-event time (and cut out region) - this might slow the code down as it re-plots the slider everytime but I can't find an efficient way to do this differently for now.
            fig.canvas.draw_idle()  # update figure

            # these values need to be global as are needed but the clicking and the slider events
            global transit_slider_ax
            global transit_slider

            transit_slider_ax.remove() # remove the old bar and plot a new one which has the colour bar in the right place

            if len(in_sec) == 1:
                transit_slider_ax  = fig.add_axes([0.24, 0.15, 0.66, 0.03])  # location of the slider
                transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=val, color='teal')
            else:
                transit_slider_ax  = fig.add_axes([0.24, 0.17, 0.66, 0.02])  # location of the slider
                transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=val, color='teal')


        # if clicked outside of the region activate the slider bar
        else:
            # Link the functions (actions) to the corresponding sliders
            transit_slider.on_changed(update_axis)
            scale_slider_upper.on_changed(update_yaxis_upper)
            scale_slider_lower.on_changed(update_yaxis_lower)
            transit_slider.on_changed(sliders_on_changed)

    # alternatively one can click on the image in order to change the zoom in region of the plot.
    fig.canvas.mpl_connect('button_press_event', onclick)

    # ---------------------

    # Define buttons which allow the user to interactively choose options:
    # run simple code, BLS, model the data usign Pyaneti, save the data and figures, creata DV report

    # only give the model option if pyaneti has been sucessfully installed
    if pyaneti_installed == True:
        var_ax = fig.add_axes([0.025, 0.52, 0.119, 0.19]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'model', 'Save', 'Report'), (False, args.noshow, args.north, False, False, True, True))
    else:
        var_ax = fig.add_axes([0.025, 0.53, 0.119, 0.18]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'Save', 'Report'), (False, args.noshow, args.north, False, True, True))

    if len(in_sec) > 1:
        sector_slider.on_changed(update_plotting_region)


    aster_ax = fig.add_axes([0.025, 0.23, 0.119, 0.045], facecolor='gold') # x, y, width, height

    asteroseismic_button = Button(aster_ax, 'Seismology', color = 'white')

    asteroseismic_colors = itertools.cycle(['gold','white'])
    asteroseismic_status = itertools.cycle([True, False])

    global asteroseismology
    asteroseismology = False

    def change_color(event):
        global asteroseismology
        asteroseismic_button.color = next(asteroseismic_colors)

        asteroseismology = next(asteroseismic_status)

        # If you want the button's color to change as soon as it's clicked, you'll
        # need to set the hovercolor, as well, as the mouse is still over it
        asteroseismic_button.hovercolor = asteroseismic_button.color
        fig.canvas.draw()

    asteroseismic_button.on_clicked(change_color)

    # Initial values for each option
    simple = False
    hide = args.noshow
    north = args.north
    BLS = False
    model = False
    save = True
    DV = True

    # function to get the status of each button and to save it
    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        hide = not status[1]
        north = status[2]
        BLS = status[3]

        if pyaneti_installed == True:
            model = status[4]
            save = status[5]
            DV = status[6]
        else:
            model = False
            save = status[4]
            DV = status[5]

    # ---------------

    # Add a set of radio buttons for changing the binning of the data
    binning_ax = fig.add_axes([0.035, 0.77, 0.1, 0.15]) # x, y, width, height
    binning_radios = RadioButtons(binning_ax, ('2', '5', '7', '10'), active=0)

    # this function accesses the binning functino.
    def binning_button(label):
        line_full_binned.set_xdata(binning(int(label))[0])
        line_full_binned.set_ydata(binning(int(label))[1])
        fig.canvas.draw_idle()

    binning_radios.on_clicked(binning_button)
    # ---------------

    # define the paramaters of the plot
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf

    ax[0].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[0].vlines(all_md, minf-1,minf + height*0.3 , colors = md_color, label = "Momentum Dump")

    ax[1].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("Time (BJD-2457000)", fontsize = main_plot_fontsize)
    ax[1].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[1].vlines(all_md, minf-1,minf + height*0.3, lw =1,  colors = md_color, label = "Momentum Dump")

    # ---------------
    # Create a label for the 'binning' box to clarify what it does
    plt.text(0.08, 1.1, "Binning Factor", fontsize=10, verticalalignment='center')
    # ---------------
    # Create a label for the 'settings' box to clarify what it does
    plt.text(0.25, -0.28, "Settings", fontsize=10, verticalalignment='center')
    # ---------------
    # Create a label for the 'quick look options'
    plt.text(0.22, -1.85, "Quick look", fontsize=10, verticalalignment='center')
    # define a box to enter the nickname of the target - this is useful if you analyse a lot of different candidates and you need a way of identifying them easily

    # ---------------
    # define a 'text box' that lets you enter the transit times

    initial_text = ""  # no initial text (the text box is empty)
    nick_name = []  # list of the entered transit times

    # function to store the entererd transit times.
    def submit(text):
        nick_name.append(text)

    axbox = plt.axes([0.24, 0.01, 0.22, 0.04]) # x, y, width, height
    text_box = TextBox(axbox, 'Memorable name (optional): ', initial=initial_text)
    text_box.on_submit(submit)
    # ---------------


    # Create a button to close the figure and more onto the next stage of the code.
    ebx = plt.axes([0.77, 0.01, 0.13, 0.04])
    exit = Button(ebx, 'Done', color=main_pot_color)

    # make button to exit the plot and continue with the code.
    # pop up warning if no transit time is entered

    def close(event):

        if len(transit_times) > 0:
            plt.close('all')
        elif (len(transit_times) == 0) and (asteroseismology == True):
            print ("asteroseismic analysis only")
            plt.close('all')

        else:
            print ("Must enter at least one transit time!")
            global deltxt

            if len(in_sec) > 1:
                deltxt = ax[1].text(3.7,-10.5,"Please enter at least one transit time!", color = 'red', size=9)
            else:
                deltxt = ax[1].text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)

            plt.draw()


    exit.on_clicked(close)
    # ---------------

    global ttxt
    global deltxt

    if len(in_sec) > 1:
        deltxt = plt.text(3.7,-10.5," ", color = 'red', size=9)
    else:
        deltxt = plt.text(3.7,-10," ", color = 'red', size=9)


    # text that will be updated to list the 'wanted' transit event times
    if len(in_sec) > 1:
        ttxt = plt.text(-4.88,1.5, "Transit times: ", weight='bold')
    else:
        ttxt = plt.text(-4.88,2, "Transit times: ", weight='bold')

    # Create a button to store the 'current' value of the slider
    stx = plt.axes([0.53, 0.01, 0.11, 0.04])
    store_val = Button(stx, 'Add time', color='gold')

    transit_times = []

    # function to store a transit event time using a button
    def storeval(time):

        storetime = (transit_slider.val)
        transit_times.append(round(storetime,2))
        global ttxt
        global deltxt

        ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
        deltxt.set_text(" ")

        plt.draw()


    store_val.on_clicked(storeval)
    # ---------------

    # Create a button to delete the last entered value
    delx = plt.axes([0.65, 0.01, 0.11, 0.04])
    del_val = Button(delx, 'Remove time', color='grey')

    transit_times = []

    # function to delete the last entry if using a button
    def deleval(time):
        global ttxt
        if len (transit_times)>0:
            del transit_times[-1]

            ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
            plt.draw()


    del_val.on_clicked(deleval)

    # ---------------

    # Make some quick look buttons to look at certain plots without having to run the program/make a report

    plotbgkx = plt.axes([0.025, 0.45, 0.12, 0.025]) # x,y, xwidth, ywidth
    plotbgk = Button(plotbgkx, 'Background', color='white')
    plotbgk.label.set_fontsize(10)

    plotcentx = plt.axes([0.025, 0.415, 0.12, 0.025]) # x,y, xwidth, ywidth
    plotcent = Button(plotcentx, 'Centroids', color='white')
    plotcent.label.set_fontsize(10)

    ploteepx = plt.axes([0.025, 0.38, 0.12, 0.025]) # x,y, xwidth, ywidth
    ploteep = Button(ploteepx, 'Stellar tracks', color='white')
    ploteep.label.set_fontsize(10)

    plotperiodx = plt.axes([0.025, 0.345, 0.12, 0.025]) # x,y, xwidth, ywidth
    plotperiod = Button(plotperiodx, 'Periodogram', color='white')
    plotperiod.label.set_fontsize(10)

    # ------
    # these require a different data download of the TPFs

    #plotapx = plt.axes([0.025, 0.3, 0.12, 0.03]) # x,y, xwidth, ywidth
    #plotap = Button(plotapx, 'Ap. size', color='gold')

    #plotpixlevelx = plt.axes([0.025, 0.25, 0.12, 0.03]) # x,y, xwidth, ywidth
    #plotpixlevel = Button(plotpixlevelx, 'Pixel square', color='gold')

    # - - - - - - - - - -
    def quick_look_bkg(time):
        global deltxt
        if len (transit_times) > 0:
            plot_background(tic, indir, alltime, allfbkg, transit_times, args, ql = True)
        else:
            deltxt = plt.text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)
            print ("select at least one transit time to run this")

    plotbgk.on_clicked(quick_look_bkg)

    # - - - - - - - - - -
    def quick_look_centroid(time):
        global deltxt
        if len (transit_times) > 0:

            plot_centroid(tic, indir, alltime12, allx1, ally1, allx2, ally2, transit_times, args, ql = True)
        else:
            deltxt = plt.text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)
            print ("select at least one transit time to run this")


    plotcent.on_clicked(quick_look_centroid)

    # - - - - - - - - - -
    def quick_look_aperturesize(time):


        if (len (transit_times) > 0):

            print ("downloading the TPF data for this plot - this will take a moment.")

            transit_sectors = transit_sec(in_sec,start_sec, end_sec, transit_times)

            _, _, _, _, _, _, _, _, _, _,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, _ = download_tpf(indir, transit_sectors, transit_times, tic, url_list)

            plot_aperturesize(tic, indir, TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, transit_times, args, ql = True)

        else:
            print ("select at least one transit time to run this")

    #plotap.on_clicked(quick_look_aperturesize)


    # - - - - - - - - - -
    def quick_look_pixellevel(time):


        if (len (transit_times) > 0):

            print ("downloading the TPF data for this plot - this will take a moment.")

            transit_sectors = transit_sec(in_sec,start_sec, end_sec, transit_times)

            X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, _, _,_, _, _, _, _, _, tpf_list = download_tpf(indir, transit_sectors, transit_times, tic)

            plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, tpf_list, apmask_list, arrshape_list, t_list, transit_times, args, ql = True)

        else:
            print ("select at least one transit time to run this")


    #plotpixlevel.on_clicked(quick_look_pixellevel)

    # - - - - - - - - - -
    def quick_look_eep(time):

        eep_target(tic, indir, syspath,teff, srad, args, ql = True)

    ploteep.on_clicked(quick_look_eep)

    # - - - - - - - - - -

    def quick_look_periodogram(time):

        plot_periodogram(tic, indir, alltime, allflux, teff, 0, 0, args, ql = True)

    plotperiod.on_clicked(quick_look_periodogram)

    # - - - - - - - - - -


    plt.show()
    plt.close('all')

    # ------ END OF INTERACTIVE PART OF CODE ------

    # save the status of all the buttons
    end_status = save_var.get_status()

    if pyaneti_installed == True:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = end_status[4]
        save = end_status[5]
        DV = end_status[6]

    else:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = False
        save = end_status[4]
        DV = end_status[5]


    # ---------------

    # update the arguments if they are different from the ones entered in the command line:
    if args.noshow != hide:
        args.noshow = hide

    if args.north != north:
        args.north = north

    # if the save button was not selected, then change the iput argument to false
    if save == False:
        args.save = save

    # check whether the user identified a nickname, or memorable name, for the candidate
    if len(nick_name) != 0: # if the box was left untouched do nothing and continue without setting a new nickname - will be ignored
        if len(nick_name[-1]) > 0: # [-1] in case the user pressed enter and then deleted the nickname again
            args.nickname = str(nick_name[-1])

    # ---------------

    # get the entered transit times and make sure that all the values are floats
    transit_list = transit_times

    print ("Transit times : {} \n".format(str(transit_list)[1:-1]))

    # get the final upper and lower axis limit - used for later plotting
    upper_axis_lim_final = scale_slider_upper.val
    lower_axis_lim_final = scale_slider_lower.val

    if asteroseismology == True:
        interact_LATTE_asteroseismology(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, args)
    else:
        #  -----  BREW  ------
        brew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, -999, -999, [-111], args)


def interact_LATTE_asteroseismology(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, args):

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

    # -------------------------
    # Plot the interactive plot - uses matplolib
    # -------------------------


    fig, ax = plt.subplots(2, 1, figsize=(11,7.5))
    plt.tight_layout()

    # Adjust the plots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.3)


    alltime = np.array(alltime)
    allflux = np.array(allflux)
    # make sure that the data doesn't contain any nans
    finite_mask = np.isfinite(alltime) & np.isfinite(allflux)


    # use the lightkurve function to calculate the periodogram
    lc = lk.lightcurve.LightCurve(time = alltime[finite_mask], flux = allflux[finite_mask])

    # remove outliers from the periodogram
    lc = lc.remove_outliers().remove_nans()
    pg = lc.to_periodogram(method='lombscargle', normalization='psd') #, minimum_frequency=1, maximum_frequency=200)

    seis = pg.to_seismology()
    numax_w_units = seis.estimate_numax()
    numax = seis.estimate_numax().value

    #ls = lc.to_periodogram(normalization = 'psd')

    # create a smoothing filter - make so that this filter width can be changed?
    smooth = pg.smooth(method='boxkernel', filter_width=20.)


    freq = pg.frequency
    power = pg.power

    freq_smooth = smooth.frequency
    power_smooth = smooth.power

    freqmin = np.nanmin(freq.value)
    freqmax = np.nanmax(freq.value)
    freq_diff = freqmax - freqmin


    def change_smoothing_width(filter_width):

        smooth = pg.smooth(method='boxkernel', filter_width=filter_width)

        return [smooth.frequency, smooth.power]

    # - - - plot - - -

    #plt.plot(freq, power, color = 'k', lw = 0.5)
    #plt.plot(freq_smooth, power_smooth, color = '#bf4080', lw = 2)
    #
    #plt.xlabel(r"$Frequency~(\mu Hz)$", fontsize = main_plot_fontsize)
    #plt.ylabel(r"$Power~Spectral~Density~(A^{2} \mu Hz^{-1})$", fontsize = main_plot_fontsize)
    #
    #plt.loglog()

    # function to define the plotting area around the selected event.
    # this needs to be in a function as the area changes with the interactive slider.
    def cutout(numax):
        mask_smooth = (np.array(freq_smooth) > numax-freq_diff/5) & (np.array(freq_smooth) < numax+freq_diff/5)
        mask = (np.array(freq) > numax-freq_diff/5) & (np.array(freq) < numax+freq_diff/5)

        return [np.array(freq)[mask], np.array(power)[mask], np.array(freq), np.array(power), np.array(freq_smooth)[mask_smooth], np.array(power_smooth)[mask_smooth], np.array(freq_smooth), np.array(power_smooth)]

    # ---------------

    # define the slider that let's you chose the sector to look at
    class PageSlider(matplotlib.widgets.Slider):

        def __init__(self, ax, label, numpages = 10, valinit=0, valfmt='%1d',
                     closedmin=True, closedmax=True,
                     dragging=True, **kwargs):

            self.facecolor=kwargs.get('facecolor',"w")
            self.activecolor = kwargs.pop('activecolor',"b")
            self.fontsize = kwargs.pop('fontsize', 10)
            self.numpages = numpages

            super(PageSlider, self).__init__(ax, label, 0, numpages,
                                valinit=valinit, valfmt=valfmt, **kwargs)

            # make it so that numpages is the index of the array of the
            global in_sec

            self.poly.set_visible(False)
            self.vline.set_visible(False)
            self.pageRects = []
            for i in range(numpages):

                facecolor = self.activecolor if i==valinit else self.facecolor
                r  = Rectangle((float(i)/numpages, 0), 1./numpages, 1,
                                    transform=ax.transAxes, facecolor=facecolor)
                ax.add_artist(r)
                self.pageRects.append(r)

                seclist_names = ['All'] + in_sec

                sec_i = seclist_names[i]

                ax.text(float(i)/numpages+0.5/numpages, 0.5, str(sec_i),
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=self.fontsize)
            self.valtext.set_visible(False)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            bax = divider.append_axes("right", size="5%", pad=0.05)
            fax = divider.append_axes("right", size="5%", pad=0.05)
            self.button_back = matplotlib.widgets.Button(bax, label=r'$\blacktriangleleft$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_forward = matplotlib.widgets.Button(fax, label=r'$\blacktriangleright$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_back.label.set_fontsize(self.fontsize)
            self.button_forward.label.set_fontsize(self.fontsize)
            self.button_back.on_clicked(self.backward)
            self.button_forward.on_clicked(self.forward)

        def _update(self, event):
            super(PageSlider, self)._update(event)
            i = int(self.val)
            if i >=self.valmax:
                return
            self._colorize(i)

        def _colorize(self, i):
            for j in range(self.numpages):
                self.pageRects[j].set_facecolor(self.facecolor)
            self.pageRects[i].set_facecolor(self.activecolor)

        def forward(self, event):
            current_i = int(self.val)
            i = current_i+1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

        def backward(self, event):
            current_i = int(self.val)
            i = current_i-1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

    # ---------------

    # define the initial conditions - this is what appears on the plot initially but can change later.
    numax = numax  # start the cut-out in the centre of the plot.
    filter_width = 20

    # FIRST PLOT - loglog plot
    # plot the unbinned flux
    [line_full] = ax[0].plot(freq, power, lw = 0.5, color = 'k', alpha = 1)
    # plot the binned flux - it is called through the above function so that the binning can be changed.
    #[line_full_smooth] = ax[0].plot(freq_smooth, power_smooth, color = main_pot_color, lw = 2, alpha = 0.9)
    [line_full_smooth] = ax[0].plot(change_smoothing_width(filter_width)[0], change_smoothing_width(filter_width)[1], color = main_pot_color, lw = 2, alpha = 0.9)

    # ---------------
    # SECOND PLOT - zoom in in linear space
    # plot the cut out around the time of the numax - called with function as it can be changed with the slider.
    [line] =  ax[1].plot(cutout(numax)[0], cutout(numax)[1], lw = 1, color = 'k', alpha = 1)
    #[line_binned] =  ax[1].plot(cutout(numax)[4], cutout(numax)[5],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='k')
    # ---------------

    global numax_slider_ax
    global numax_slider
    global analysis_slider

    # Define the slider to change the y axis scalex
    scale_slider_ax  = fig.add_axes([0.24, 0.2, 0.66, 0.02], facecolor='silver')
    scale_slider = Slider(scale_slider_ax, 'display window', 20, freq_diff/5, valinit=200, color='white')
    scale_slider.valtext.set_visible(False)

    # Define the slider to change the numax-event time (and cut out region)
    analysis_slider_ax  = fig.add_axes([0.24, 0.17, 0.66, 0.02])  # location of the slider
    analysis_slider = Slider(analysis_slider_ax, 'analysis window', 1, 350, valinit=50, color='teal')

    # Define the slider to change the numax-event time (and cut out region)
    numax_slider_ax  = fig.add_axes([0.24, 0.14, 0.66, 0.02])  # location of the slider
    numax_slider = Slider(numax_slider_ax, 'nu max', np.nanmin(freq.value), np.nanmax(freq.value), valinit=numax, color='teal')

    # ---------------

    # deifne the intial x and y axis limits - y limit can be changed with slider, x limit cannot.
    ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
    ax[0].set_xlim([10, np.nanmax(freq.value)])

    ax[0].loglog()
    ax[1].set_xlim([numax -200, numax + 200])

    # ---------------
    # Define an action for acessing the required cut-out data and drawing the line when the slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(numax_slider.val)[0])
        line.set_ydata(cutout(numax_slider.val)[1])

        fig.canvas.draw_idle()  # update figure

    # draw a line at the time of the chosen numax-event (from slider) on the top and bottom plot
    lver0 = ax[0].axvline(numax, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(numax, color = 'r', linewidth = 2)

    # draw the lines for the limits - these are shown in blue and changed with a different slider (and not with the click)

    analysis_lower_lver0 = ax[0].axvline(numax - 20, color = md_color, linewidth = 1.5)
    analysis_lower_lver1 = ax[1].axvline(numax - 20, color = md_color, linewidth = 1.5)

    analysis_upper_lver0 = ax[0].axvline(numax + 20, color = md_color, linewidth = 1.5)
    analysis_upper_lver1 = ax[1].axvline(numax + 20, color = md_color, linewidth = 1.5)

    # keep track of the order that the lines were plotted in - start with numax, limit, limit

    # ---------------
    # Define an action for modifying the plot region on the second plot
    def update_axis(val):

        # this function
        # ---------------
        # only ever plot one line so get rid of the old ones before plotting the new ones...
        if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
            # have to do this three times because we need to remove all of the lines from the numax, and the two limits (and I don't know how to condense to make it one line of code so here we are)
            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()


        # draw a line at the time of the chosen numax-event (from slider) on the top and bottom plot
        lver0 = ax[0].axvline(numax_slider.val, color = 'r', linewidth = 2)
        lver1 = ax[1].axvline(numax_slider.val, color = 'r', linewidth = 2)

        # update the analysis window lines too
        analysis_lim      =  analysis_slider.val
        numax_value       =  numax_slider.val

        analysis_lower_lver0 = ax[0].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)
        analysis_lower_lver1 = ax[1].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)

        analysis_upper_lver0 = ax[0].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)
        analysis_upper_lver1 = ax[1].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)


        ax[1].set_xlim([numax_slider.val - scale_slider.val, numax_slider.val + scale_slider.val])
        lver0.set_xdata(numax_slider.val)
        lver1.set_xdata(numax_slider.val)

    # draw a line at the time of the chosen numax-event (from slider) on the top and bottom plot
    #lver0 = ax[0].axvline(numax, color = 'r', linewidth = 2)
    #lver1 = ax[1].axvline(numax, color = 'r', linewidth = 2)
    # starting positions for the sliders


    # ---------------
    # Define an action for modifying the plot region on the second plot
    def analysis_xaxis(val):
        # ---------------
        # only ever plot one line so get rid of the old ones before plotting the new ones...
        if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
            # have to do this three times because we need to remove all of the lines from the numax, and the two limits (and I don't know how to condense to make it one line of code so here we are)
            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()


        # draw a line at the time of the chosen numax-event (from slider) on the top and bottom plot
        lver0 = ax[0].axvline(numax_slider.val, color = 'r', linewidth = 2)
        lver1 = ax[1].axvline(numax_slider.val, color = 'r', linewidth = 2)

        # update the analysis window lines too
        analysis_lim      =  analysis_slider.val
        numax_value       =  numax_slider.val

        analysis_lower_lver0 = ax[0].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)
        analysis_lower_lver1 = ax[1].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)

        analysis_upper_lver0 = ax[0].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)
        analysis_upper_lver1 = ax[1].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)


        ax[1].set_xlim([numax_slider.val - scale_slider.val, numax_slider.val + scale_slider.val])
        lver0.set_xdata(numax_slider.val)
        lver1.set_xdata(numax_slider.val)


    def update_xaxis(val):

        axis_lim = scale_slider.val
        numax_value   = numax_slider.val

        ax[1].set_xlim([numax_value - axis_lim, numax_value + axis_lim])


    # ------ ON CLICK -------
    # define the action that lets you click on the image in order to chose the location instead of using the slider - probably more useful
    def onclick(event):

        # just need to know the x (time) position of the click. The y (flux) position doesn't matter.
        val = event.xdata

        # check that the click was within the plotting region. If not ignore the click
        if (event.inaxes == ax[0]) or (event.inaxes == ax[1]):

            # SECOND PLOT - CUT OUT LC
            # plot the cut out around the time of the numax - changed by clicked on the image (either one of them)
            line.set_xdata(cutout(val)[0])
            line.set_ydata(cutout(val)[1])


            if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
                # have to do this three times because we need to remove all of the lines from the numax, and the two limits (and I don't know how to condense to make it one line of code so here we are)
                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()


            # draw a line at the time of the chosen numax-event (from slider) on the top and bottom plot
            lver0 = ax[0].axvline(val, color = 'r', linewidth = 2)
            lver1 = ax[1].axvline(val, color = 'r', linewidth = 2)

            # update the analysis window lines too
            analysis_lim      =  analysis_slider.val
            numax_value       =  val

            analysis_lower_lver0 = ax[0].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)
            analysis_lower_lver1 = ax[1].axvline(numax_value - analysis_lim, color = md_color, linewidth = 1.5)

            analysis_upper_lver0 = ax[0].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)
            analysis_upper_lver1 = ax[1].axvline(numax_value + analysis_lim, color = md_color, linewidth = 1.5)


            ax[1].set_xlim([val - scale_slider.val, val + scale_slider.val])
            lver0.set_xdata(val)
            lver1.set_xdata(val)


            # also update the location of the slider (aesthetic reasons only)
            # Define the slider to change the numax-event time (and cut out region) - this might slow the code down as it re-plots the slider everytime but I can't find an efficient way to do this differently for now.
            fig.canvas.draw_idle()  # update figure

            # these values need to be global as are needed but the clicking and the slider events
            global numax_slider_ax
            global numax_slider

            numax_slider_ax.remove() # remove the old bar and plot a new one which has the colour bar in the right place

            numax_slider_ax  = fig.add_axes([0.24, 0.14, 0.66, 0.02])  # location of the slider
            numax_slider = Slider(numax_slider_ax, 'numax', np.nanmin(freq.value), np.nanmax(freq.value), valinit=val, color='teal')


        # if clicked outside of the region activate the slider bar
        else:
            # Link the functions (actions) to the corresponding sliders
            numax_slider.on_changed(update_axis)
            scale_slider.on_changed(update_xaxis)
            analysis_slider.on_changed(analysis_xaxis)
            numax_slider.on_changed(sliders_on_changed)

    # alternatively one can click on the image in order to change the zoom in region of the plot.
    fig.canvas.mpl_connect('button_press_event', onclick)

    # ---------------

    # Add a set of radio buttons for changing the binning of the data
    smoothing_length_ax = fig.add_axes([0.02, 0.77, 0.1, 0.15]) # x, y, width, height
    smoothing_length_radios = RadioButtons(smoothing_length_ax, ('10', '20', '40', '60'), active=0)

    # this function accesses the binning functino.
    def smoothing_button(label):
        line_full_smooth.set_xdata(change_smoothing_width(int(label))[0])
        line_full_smooth.set_ydata(change_smoothing_width(int(label))[1])
        fig.canvas.draw_idle()

    smoothing_length_radios.on_clicked(smoothing_button)


    # ---------------------

    # define the paramaters of the plot
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf


    #plt.text(4.1,-4.7, r"Change the value of $\nu_{max}$ and the analysis window to desired values", color = 'k', size=11)
    plt.text(2.57,-4.7, r"---------------------- Quick look ----------------------", color = 'k', size=11)

    ax[0].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')

    #ax[0].set_ylabel(r"$Power~Spectral~Density~(A^{2} \mu Hz^{-1})$", fontsize = main_plot_fontsize)

    ax[1].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_ylabel(r"$Power~Spectral~Density~(A^{2} \mu Hz^{-1})$", fontsize = main_plot_fontsize)
    ax[1].set_xlabel(r"$Frequency~(\mu Hz)$", fontsize = main_plot_fontsize)

    # Create a label for the 'binning' box to clarify what it does
    plt.text(0.12, 1.1, "Filter width", fontsize=10, verticalalignment='center')


    qlx = plt.axes([0.24, 0.01, 0.11, 0.04])
    ql = Button(qlx, 'Nu max', color='white')

    qlx2 = plt.axes([0.36, 0.01, 0.11, 0.04])
    ql2 = Button(qlx2, 'Delta nu', color='white')

    qlx3 = plt.axes([0.48, 0.01, 0.11, 0.04])
    ql3 = Button(qlx3, 'Echelle', color='white')


    global radtxt
    global masstxt
    global loggtxt
    global numaxtxt
    global deltanutxt

    mass_ast_text = ''
    radius_ast_text = ''
    logg_ast_text = ''
    numax_aast_text = ''
    deltanu_aast_text = ''

    tefftxt = plt.text(-4.15,17, r"$T_{{eff}}: {} K$".format(str(round(teff, 2))), weight='bold')
    radtxt_tic =  plt.text(-4.15,16, r"$Radius_{{tic}}: {}$".format(round(srad, 3)), weight='bold')
    radtxt =  plt.text(-4.15,15, r"$Radius_{{ast}}: $", weight='bold')
    masstxt = plt.text(-4.15,14, r"$Mass_{{ast}}: $", weight='bold')
    loggtxt = plt.text(-4.15,13, r"$logg: $", weight='bold')
    numaxtxt = plt.text(-4.15,12, r"$\nu_{\rm max}: $", weight='bold')
    deltanutxt = plt.text(-4.15,11, r"$\Delta \nu: $", weight='bold')

    tmagtxt = plt.text(-4.15,9, r"$Tmag: {}$".format(str(round(tessmag, 2))), weight='bold')
    ratxt = plt.text(-4.15,8, r"$RA: {}$".format(str(round(ra, 5))), weight='bold')
    dectxt = plt.text(-4.15,7, r"$Dec: {}$".format(str(round(dec, 5))), weight='bold')


    def quick_look(time):

        global radtxt
        global masstxt
        global loggtxt
        global numaxtxt
        global deltanutxt

        pg_ql = lc.to_periodogram(method='lombscargle', normalization='psd', minimum_frequency=numax_slider.val-analysis_slider.val, maximum_frequency=numax_slider.val+analysis_slider.val)
        snr = pg_ql.flatten()
        seis = snr.to_seismology()

        if len(ax[1].texts) > 0:
            for current_text in ax[1].texts:
                current_text.set_visible(False)

        # determine delta nu and numax - these are needed for the diagnostic plots
        # these are LK objects that have lots of information in them - like magic...
        try:
            numax = seis.estimate_numax()
            deltanu = seis.estimate_deltanu()

            mass_ast_text = str(round(seis.estimate_mass(teff).value,3))
            radius_ast_text = str(round(seis.estimate_radius(teff).value,3))
            logg_ast_text = str(round(seis.estimate_logg(teff).value,3))

            numax_ast_text = str(round(numax.value,3))
            deltanu_ast_text = str(round(deltanu.value,3))


            radtxt.set_text(r"Radius: {} $R_{{\odot}}$".format(radius_ast_text))
            masstxt.set_text(r"Mass: {} $M_{{\odot}}$".format(mass_ast_text))
            loggtxt.set_text(r"logg: {} $dex$".format(logg_ast_text))
            numaxtxt.set_text(r"$\nu_{{\rm max}}: {}$".format(numax_ast_text))
            deltanutxt.set_text(r"$\Delta \nu: {}$".format(deltanu_ast_text))

            plt.draw()
            plot_diagnostic_numax(tic, indir, numax, deltanu, pg_ql, args, ql = True)

        except:
            print ("Asteroseismic analysis did not work. Try increasing the analysis window.")
            ax[1].text(0.5,-0.6, "Asteroseismic analysis did not work. Try increasing the analysis window.", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            current_text = ax[1].texts
            current_text[0].set_visible(True)
            fig.canvas.draw()


    def quick_look2(time):

        global radtxt
        global masstxt
        global loggtxt
        global numaxtxt
        global deltanutxt

        pg_ql = lc.to_periodogram(method='lombscargle', normalization='psd', minimum_frequency=numax_slider.val-analysis_slider.val, maximum_frequency=numax_slider.val+analysis_slider.val)
        snr = pg_ql.flatten()
        seis = snr.to_seismology()

        if len(ax[1].texts) > 0:
            for current_text in ax[1].texts:
                current_text.set_visible(False)

        # determine delta nu and numax - these are needed for the diagnostic plots
        # these are LK objects that have lots of information in them - like magic...
        try:
            numax = seis.estimate_numax()
            deltanu = seis.estimate_deltanu()

            mass_ast_text = str(round(seis.estimate_mass(teff).value,3))
            radius_ast_text = str(round(seis.estimate_radius(teff).value,3))
            logg_ast_text = str(round(seis.estimate_logg(teff).value,3))

            numax_ast_text = str(round(numax.value,3))
            deltanu_ast_text = str(round(deltanu.value,3))

            radtxt.set_text(r"Radius: {} $R_{{\odot}}$".format(radius_ast_text))
            masstxt.set_text(r"Mass: {} $M_{{\odot}}$".format(mass_ast_text))
            loggtxt.set_text(r"logg: {} $dex$".format(logg_ast_text))
            numaxtxt.set_text(r"$\nu_{{\rm max}}: {}$".format(numax_ast_text))
            deltanutxt.set_text(r"$\Delta \nu: {}$".format(deltanu_ast_text))

            plt.draw()
            plot_scaled_correlation(tic, indir, deltanu, pg_ql, args, ql = True)


        except:
            print ("Asteroseismic analysis did not work. Try increasing the analysis window.")
            ax[1].text(0.5,-0.6, "Asteroseismic analysis did not work. Try increasing the analysis window.", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            current_text = ax[1].texts
            current_text[0].set_visible(True)
            fig.canvas.draw()


    def quick_look3(time):

        global radtxt
        global masstxt
        global loggtxt
        global numaxtxt
        global deltanutxt

        pg_ql = lc.to_periodogram(method='lombscargle', normalization='psd', minimum_frequency=numax_slider.val-analysis_slider.val, maximum_frequency=numax_slider.val+analysis_slider.val)
        snr = pg_ql.flatten()
        seis = snr.to_seismology()

        if len(ax[1].texts) > 0:
            for current_text in ax[1].texts:
                current_text.set_visible(False)

        # determine delta nu and numax - these are needed for the diagnostic plots
        # these are LK objects that have lots of information in them - like magic...
        try:
            numax = seis.estimate_numax()
            deltanu = seis.estimate_deltanu()

            mass_ast_text = str(round(seis.estimate_mass(teff).value,3))
            radius_ast_text = str(round(seis.estimate_radius(teff).value,3))
            logg_ast_text = str(round(seis.estimate_logg(teff).value,3))

            numax_ast_text = str(round(numax.value,3))
            deltanu_ast_text = str(round(deltanu.value,3))

            radtxt.set_text(r"Radius: {} $R_{{\odot}}$".format(radius_ast_text))
            masstxt.set_text(r"Mass: {} $M_{{\odot}}$".format(mass_ast_text))
            loggtxt.set_text(r"logg: {} $dex$".format(logg_ast_text))
            numaxtxt.set_text(r"$\nu_{{\rm max}}: {}$".format(numax_ast_text))
            deltanutxt.set_text(r"$\Delta \nu: {}$".format(deltanu_ast_text))

            plt.draw()
            plot_echelle(tic, indir, seis, args,ql = True)

        except:
            print ("Asteroseismic analysis did not work. Try increasing the analysis window.")
            ax[1].text(5,-0.6, "Asteroseismic analysis did not work. Try increasing the analysis window.", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            current_text = ax[1].texts
            current_text[0].set_visible(True)
            fig.canvas.draw()


    ql.on_clicked(quick_look)
    ql2.on_clicked(quick_look2)
    ql3.on_clicked(quick_look3)

    # Create a button to close the figure and more onto the next stage of the code.
    ebx_nosig = plt.axes([0.67, 0.01, 0.11, 0.04])
    exit_nosig = Button(ebx_nosig, 'No signal', color='grey')

     # Create a button to close the figure and more onto the next stage of the code.
    ebx = plt.axes([0.79, 0.01, 0.11, 0.04])
    exit = Button(ebx, 'Done', color=main_pot_color)


       # make button to exit the plot and continue with the code.
    # pop up warning if no numax time is entered

    global final_numax
    global final_analysis_window

    def close(event):
        plt.close('all')

        global final_numax
        global final_analysis_window

        final_numax = numax_slider.val
        final_analysis_window = analysis_slider.val

        print ("numax : {} ".format(str(round(final_numax,3))))
        print ("frequency window : {} \n".format(str(round(final_analysis_window,3))))



    def close_nosig(event):
        plt.close('all')

        global final_numax
        global final_analysis_window

        final_numax = -999
        final_analysis_window = -999

        print ("Continue without asteroseismic analysis.")

        if len (transit_list) == 0:
            print ("No analysis selected (please select at least one transit event or choose asteroseismic anlaysis). Quit prgram.")
            sys.exit('')


    exit.on_clicked(close)
    exit_nosig.on_clicked(close_nosig)

    # ---------------
    # ---------------
    plt.show()
    plt.close('all')

    # ------ END OF INTERACTIVE PART OF CODE ------
    # ---------------------------------------------

    # get the entered numax times and make sure that all the values are floats

    #  -----  BREW  ------
    brew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, final_numax, final_analysis_window,[-111], args)


# interactive tool to identify the aperture masks when run in the FFI mode
# This function is called from __main__.py
def interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, args):

    '''
    Function to run the Interactive LATTE code for the FFIs using the matplotlib interactive tool.
    This code allows the user to identify the aperture size for the small and large apertures which are then used for the LC extraction.
    Normalises the data and perfoms PCA analysis in order to correct for systematics.

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
    ra   :  float
        the right ascension of the target

    Returns
    -------
        saves the correct and uncorrected LCs (for comaprison to ensure that the detrending is oworking correctiyl) and the chosen apertures.


    alltime  :  list
        times
    allflux  :  list
        normalized flux extracted with the larger aperture (PCA corrected)
    allflux_small  : list
        normalized flux extracted with the smaller aperture (PCA corrected)
    allflux_flat   : list
        normalized detrended flux extracted with the larger aperture
    all_md  :  list
        times of the momentum dumps
    allfbkg  :  list
        background flux
    allfbkg_t  :  list
        times used to plot the background
    start_sec  :  list
        times of the start of the sector
    end_sec  :  list
        times of the end of the sector
    in_sec  :  list
        the sectors for which data was downloaded
    X1_list  :  list
        flux vs time for each pixel (for each sector)
    X4_list  :  list
        PCA corrected flux vs time for each pixel (for each sector)
    apmask_list  :  list
        aperture masks from the pipeline
    arrshape_list  :  list
        list of the shape of the array (for each sector)
    tpf_filt_list   : list
        list of the filtered (masked) corrected target pixel data - from X4. (for each sector)
    t_list  :  list
        list of the time arrays (for each sector)
    bkg_list  :  list
        the flux that was used to normalise each pixel - i.e. what is used to make the background plot colour for each pixel.
    tpf_list   : list
        list of the target pixel files (for each sector)
    '''
    # call function to download the data from MAST
    print ("Start TPF data download.....", end =" ")
    alltime_list, all_md, start_sec, end_sec, in_sec, X1_list, X1flux_list,  X4_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI_interact(indir, sectors, sectors_all, tic, save = False)
    print ("done.\n")

    # define the button that closes the event
    # pop up warning if no aperture is selected
    def close(event):
        if (np.sum(aperture) > 0 ) * (np.sum(aperture2) > 0 ):
            plt.close('all')
        else:
            print ("Must select at least one pixel per aperture!")
            ax[2].text(1.1,-0.29, "Select at least one pixel per aperture!", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            fig.canvas.draw()

    # ---------------
    # function to extract the LC using the defined large aperture
    def extract_LC(aperture):
        ax[0].cla()
        flux = X4[:,aperture.flatten()].sum(axis=1)
        m = np.nanmedian(flux)
        return t, flux/m

    # function to extract the LC using the defined small aperture
    def extract_LC2(aperture):
        flux = X4[:,aperture.flatten()].sum(axis=1)
        m = np.nanmedian(flux)

        return t, flux/m
    # ---------------

    allfbkg = []
    allfbkg_t = []

    allflux = []
    alltime = []
    allflux_flat = []
    allflux_small = []
    apmask_list = []

    for index, X4 in enumerate(X4_list):
        tpf = tpf_list[index]
        t = t_list[index]
        X1 = X1flux_list[index]
        sec = sectors[index]

        global mask
        global mask2
        global aperture
        global aperture2

        mask = []
        mask2 = []

        # Define the aperture arrays - initialise with all 'False' (i.e. no aperture)
        # This is only for the start - apertures are re-defined with every click
        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
        aperture2 = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)

        # place a random aperture in the centre to start off with - just as an example.
        for i in range(int(len(aperture)/2 ), int(len(aperture)/2 + 2)):
            for j in range(int(len(aperture)/2 ), int(len(aperture)/2 + 2)):
                mask.append((i,j))
        # ---------------

        # Define a function that registers where the user clicks on the aperture - what pixel you select.
        def onclick(event):
            global mask
            global mask2
            global aperture
            global aperture2

            # event.xdata and event.ydata are the locations of the clicks on the plot
            # need to add 0.5 as the pixels are shifted relative to the axes.
            events = ((int(event.xdata+0.5), int(event.ydata+0.5)))

            # check whether the click is on the defined axis.
            if event.inaxes in [ax[1]]:

                # the 'patches' are the highlighted select pixels
                # remove all the previous highlighted pixels and replot them all
                [p.remove() for p in reversed(ax[1].patches)]

                # if the pixel has already been selected, get rid of it from the list instead of adding it.
                # otherwise add it to the list of desired mask pixels.
                if (len(mask) > 0) and (events in list(mask)):
                    mask = [x for x in mask if x != events]
                else:
                    mask.append(events) # otherwise it's a new event and it should be added.

                # define a colour and alpha for the highlighted pixels.
                sqcol = '#ffffee'
                alpha = 0.5

                # highlight all the selected pixels
                for pixel in mask:
                    m = int(pixel[0])
                    n = int(pixel[1])
                    r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
                    ax[1].add_patch(r)

                # update figure.
                fig.canvas.draw()

            # ---------------
            # Repeat for the second mask
            if event.inaxes in [ax[2]]:
                [p.remove() for p in reversed(ax[2].patches)]

                if (len(mask2) > 0) and (events in list(mask2)):
                    mask2 = [x for x in mask2 if x != events]
                else:
                    mask2.append(events)

                for pixel2 in mask2:
                    m = int(pixel2[0])
                    n = int(pixel2[1])
                    r2 = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='#c33c7d', facecolor='#e7b1cb', alpha = 0.6)
                    ax[2].add_patch(r2)

                # update figure
                fig.canvas.draw()
            # ---------------

            # redefine an aperture with every click - all False
            aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
            aperture2 = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)

            # fill in with True depending on the selected pixels
            for coord in mask:
                aperture[coord] = True

            for coord2 in mask2:
                aperture2[coord2] = True
            # ---------------

            # plot the LC to the right of the apertures plots with the selected (highlighted) apertures.

            ax[0].plot(extract_LC(aperture)[0], extract_LC(aperture)[1],marker='o',color = '#255c74', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#054950')
            ax[0].plot(extract_LC2(aperture2)[0], extract_LC2(aperture2)[1],marker='x',color = '#cc0066', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#c94f8a')

            ax[0].set_xlabel("Time")
            ax[0].set_ylabel("Normalized Flux")

            # update plot
            fig.canvas.draw_idle()
            plt.draw()

        # ---------------
        # define the plotting parameters
        fig, ax = plt.subplots(1,3, figsize=(10,3), gridspec_kw={'width_ratios': [3, 1, 1]})

        sqcol = '#ffffee' # this is for the initial plot
        alpha = 0.5

        # ---------------
        # start off with the plotting of the larger aperture as defined above.
        for pixel in mask:
            m = int(pixel[0])
            n = int(pixel[1])
            r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
            ax[1].add_patch(r)

        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)

        for coord in mask:
            aperture[coord] = True

        # this is the initial LC plot - ipdates with every click.
        ax[0].plot(extract_LC(aperture)[0], extract_LC(aperture)[1],marker='o',color = '#054950', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#054950')

        # ---------------

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
        exit = Button(ebx, 'Done', color=main_pot_color)
        exit.on_clicked(close)

        plt.show()

        # ---- END OF INTERACTIVE APERTURE LOT ----

        # redefine the large aperture in case it wasn't changed from the 'inital' startig one - otherwise will get an error.
        aperture = np.array(np.zeros_like(X1.sum(axis=0)), dtype=bool)
        for coord in mask[:-1]:
            aperture[coord] = True

        # save the aperture mask
        apmask_list.append(aperture)

        # if they didn't mark a small aperture, just move on an make the small aperture the same size as the large aperture
        if (np.sum(aperture2) == 0):
            target_mask_small = aperture
        else:
            target_mask_small = aperture2
        # ---------------

        # rename it and plot and save the chosen apertures
        target_mask = aperture

        # ----------
        # plot the mean image and plot the extraction apertures on top of it so that one can verify that the used apertures make sense
        im = np.nanmean(tpf.flux, axis = 0)
        # set up the plot - these are stored and one of the images saved in the report
        fig, ax = plt.subplots(1,2, figsize=(10,5), subplot_kw={'xticks': [], 'yticks': []})
        kwargs = {'interpolation': 'none', 'vmin': im.min(), 'vmax': im.max()}
        color = ['red', 'deepskyblue']
        label = ['small (~60 %)', 'pipeline (100 %)']

        # plot both the large and the small apertures
        for idx,arr in enumerate([target_mask_small,target_mask]):

            mask = np.zeros(shape=(arr.shape[0], arr.shape[1]))
            mask= arr

            f = lambda x,y: mask[int(y),int(x)]
            g = np.vectorize(f)

            x = np.linspace(0,mask.shape[1], mask.shape[1]*100)
            y = np.linspace(0,mask.shape[0], mask.shape[0]*100)
            X, Y= np.meshgrid(x[:-1],y[:-1])
            Z = g(X[:-1],Y[:-1])

            ax[idx].set_title('Aperture: {}'.format(label[idx]), fontsize = 18)
            ax[idx].imshow(im, cmap=plt.cm.viridis, **kwargs, origin = 'lower')
            ax[idx].contour(Z, [0.5], colors=color[idx], linewidths=[4],
                        extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])

        # save the figure
        plt.savefig('{}/{}/{}_apertures_{}.png'.format(indir, tic, tic, index), format='png', bbox_inches = 'tight')
        plt.clf()
        plt.close('all')

        # ---------------

        # extract the flux lightcurves using the above defined apertures and
        flux = X4[:,target_mask.flatten()].sum(axis=1)
        flux_small = X4[:,target_mask_small.flatten()].sum(axis=1)

        m = np.nanmedian(flux)
        m_small = np.nanmedian(flux_small)

        # normalize the flux by dividing by the median value
        flux = flux/m
        flux_small = flux_small/m_small

        print ("done.\n")

        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---
        # ---------- flatten the large aperture lighcurve ----------

        print ("Flatten LC...", end =" ")

        fr_inj = flux

        T_dur = 0.7  #The transit duration - may need to change but this is a good average duration

        nmed = int(48*3*T_dur)  # 48 because the data has a cadence of 30 minutes and 24 hours / 30 mins = 48 data points
        nmed = 2*int(nmed/2)+1 # make it an odd number
        ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
        # ^ first number (nmed) should be roughly three time transit durations, the second quite small (10,20)

        # mask to only use the finite values - NaNs disrupt the code
        l = np.isfinite(ff)

        g = interp1d(t[l],ff[l],bounds_error=False,fill_value=np.nan)
        ff = g(t)

        fr = fr_inj / ff

        # --- do some sigma clipping to make the LC look better ----
        MAD = median_absolute_deviation(fr,ignore_nan = True)
        madrange = (100 * MAD)
        ymask = (fr < 1 + madrange) * (fr > 1 - madrange)
        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---

        # plot corrected and uncorrected LC for comparison - allows the user to check whether the detrending worked correctly.
        fix, ax = plt.subplots(3,1,figsize=(16,12))

        ax[0].plot(t,fr_inj, '.',label = 'Uncorrected')
        ax[0].plot(t,ff,'.',label = 'Model fit')
        ax[0].legend(fontsize = 16, loc = 1)

        ax[1].plot(t[~ymask],fr[~ymask], 'k.', markersize = 10, label = 'Clipped')
        ax[1].plot(t[ymask],fr[ymask], 'r.',label = 'Corrected')
        ax[1].legend(fontsize = 16, loc = 1)

        ax[2].plot(t[ymask],fr[ymask], '.', color = 'navy', label = 'Clipped + corrected')
        ax[2].legend(fontsize = 16, loc = 1)

        ax[2].set_xlabel("Time", fontsize = 16)
        ax[0].set_ylabel("Flux", fontsize = 16)
        ax[1].set_ylabel("Flux", fontsize = 16)
        ax[2].set_ylabel("Flux", fontsize = 16)

        plt.savefig('{}/{}/{}_fit_test.png'.format(indir, tic, tic), format='png')
        plt.clf()
        plt.close()

        # ---------------------------------------------
        # -------- extract the backrgound flux --------
        print ("done.\n")

        print ("Extract background...", end =" ")

        # the background flux has a low flux threshhold - pixels with low flux are assumed to be background and not stars.
        background_mask = ~tpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
        n_background_pixels = background_mask.sum()

        background_lc_per_pixel = tpf.to_lightcurve(aperture_mask=background_mask) / n_background_pixels
        n_target_pixels = target_mask.sum()
        background_estimate_lc = background_lc_per_pixel * n_target_pixels

        allfbkg.append(background_estimate_lc.flux) # uncorrected background
        allfbkg_t.append(background_estimate_lc.time)

        print ("done.\n")
         # ------------------------------------------------

        # add all the information to lists that are returned at the endof this function
        allflux_flat.append(list(fr))
        allflux.append(list(flux))
        allflux_small.append(list(flux_small))
        alltime.append(list(t))

    # flatten the lists that are lists of lists.
    allflux_flat =  list(np.hstack(allflux_flat) )
    allflux_small = list(np.hstack(allflux_small))
    alltime =       list(np.hstack(alltime)      )
    allflux =       list(np.hstack(allflux)      )
    allfbkg =       list(np.hstack(allfbkg)      )
    allfbkg_t =     list(np.hstack(allfbkg_t)    )

    # we're not back to where the same place as with interact_LATTE_FFI

    return alltime, allflux, allflux_small, allflux_flat, all_md, allfbkg,allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list


def interact_LATTE_FFI_OLD(tic, indir, syspath, sectors_all, sectors, ra, dec, args):
    '''
    Function to run the Interactive LATTE code for the FFI using the matplotlib interactive tool.
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
    ra  : float
        right ascension of the target
    dec   :
        declinatino of the target.

    *args
        if auto = False, the aperture will be chosen automatically using an intensity threshhold.

    '''

    '''
    # ----- download the data ------
    # by default you have to manually choose the different aperture sizes

    # CHECK what ffi data products are available
    data_urls, args = check_FFI_data_products(tic,sector, args)

    # see whether there is SPOC/QLP data
    if args.FFI == 'SPOC':
        alltime, allflux, allflux_err, all_md, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, url_list, args = download_data_FFI_SPOC(indir, sector, syspath, sectors_all, tic, args, save = False)

    elif args.FFI == 'QLP':
        alltime, allflux, allflux_small, allflux_err, all_md, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, X1_list, X4_list, t_list, tpf_list, tpf_filt_list,arrshape_list = download_data_FFI_QLP(indir, sector, syspath, sectors_all, tic, args, data_urls, save = False)
        alltime0 = alltime
        allflux_list = allflux
        bkg_list = [-111]
        apmask_list = [-111]
        allfbkg_t = alltime

        alltime, allflux_list, allflux_small, allflux, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, args
    # the above function returns [-111] if there is no SPOC or QLP data
    elif (args.FFI == False) and (args.auto == False):
        alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, args)

    # if the auto mode is selected in the command line, the aperture sizes are chosen automatically based on an intensity threshhold.
    elif (args.FFI == False) and (args.auto == True):
        print ("Start TPF data download.....", end =" ")
        alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t ,start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI(indir, sectors, syspath, sectors_all, tic, args)
        print ("done.\n")

    # make sure all the plots are closed before starting the next section
    plt.close('all')


    # --------------------------------------------
    # Plot the interactive plot - uses matplolib -
    # --------------------------------------------

    # allflux is detrended

    # ------- median absolute deviation in order to determine the clipping
    if args.FFI == False:
        MAD = median_absolute_deviation(allflux0, ignore_nan = True)

        madrange = (100 * MAD)

        # make a mask for the points we want to get rid of
        ymask = (allflux0 < 1 + madrange) * (allflux0 > 1 - madrange)


        alltime = np.array(alltime0)[ymask]
        allflux = np.array(allflux0)[ymask]

        # -------------------------------------

    fig, ax = plt.subplots(2, 1, figsize=(10,7))
    plt.tight_layout()

    # Adjust tbplots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.3)

    fluxmin = np.nanmin(allflux)
    fluxmax = np.nanmax(allflux)
    flux_diff = fluxmax - fluxmin
    # function to define the plotting area around the transit event.
    # this needs to be in a function as the area changes with the interactive slider.
    def cutout(transit):
        mask = (np.array(alltime) > transit-2) & (np.array(alltime) < transit+2)

        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux), np.array(allflux)]
    # ---------------

    transit = np.nanmean(alltime)
    # FIRST PLOT - FULL LC
    # plot the unbinned flux
    [line_full] = ax[0].plot(alltime, allflux , marker='o',lw = 0, markersize = 4, color = '#255c74', alpha = 0.8, label = 'unbinned', markerfacecolor='#003941')

    # SECOND PLOT - CUT OUT LC
    # plot the cut out around the time of the transit - called with function as it can be changed with the slider.
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',lw = 0, markersize = 4, color = '#255c74', alpha = 0.8, label = 'unbinned', markerfacecolor='#003941')
    # ---------------
    global transit_slider_ax
    global transit_slider

    # Define the slider to change the transit-event time (and cut out region)
    transit_slider_ax  = fig.add_axes([0.25, 0.14, 0.65, 0.03])
    transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=transit, color='#255c74')

    # Define the slider to change the y axis scale
    scale_slider_upper_ax  = fig.add_axes([0.25, 0.19, 0.325, 0.03], facecolor = 'silver') # x, y, width, height
    scale_slider_upper = Slider(scale_slider_upper_ax, 'Flux limits', 1, fluxmax + (flux_diff*0.1), valinit=fluxmax, color='white')
    scale_slider_upper.valtext.set_visible(False)

    # Define the slider to change the y axis scale
    scale_slider_lower_ax  = fig.add_axes([0.575, 0.19, 0.325, 0.03], facecolor = 'white')  # x, y, width, height
    scale_slider_lower = Slider(scale_slider_lower_ax, ' ', fluxmin - (flux_diff*0.1), 1,valinit=fluxmin, color='silver')
    scale_slider_lower.valtext.set_visible(False)

    # define the intial x and y axis limits - y limit can be changed with slider, x limit cannot.
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])

    ax[1].set_xlim([np.nanmean(alltime)-2, np.nanmean(alltime)+2])
    ax[1].set_ylim([fluxmin, fluxmax])

    # ---------------
    # Define an action for acessing the required cut-out data and drawing the line when the slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])

        fig.canvas.draw_idle()

    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)

    # ---------------
    # Define an action for modifying the plot region on the second plot
    def update_axis(val):
        ax[1].set_xlim([transit_slider.val - 2,transit_slider.val + 2])

        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)

    def update_yaxis_upper(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])

    def update_yaxis_lower(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])

    # ---------------

    # Link the functions (actions) to the corresponding sliders
    transit_slider.on_changed(update_axis)

    scale_slider_upper.on_changed(update_yaxis_upper)
    scale_slider_lower.on_changed(update_yaxis_lower)

    transit_slider.on_changed(sliders_on_changed)



    # ------ ON CLICK -------
    # define the action that lets you click on the image in order to chose the location instead of using the slider - probably more useful
    def onclick(event):

        # just need to know the x (time) position of the click. The y (flux) position doesn't matter.
        val = event.xdata

        # check that the click was within the plotting region. If not ignore the click
        if (event.inaxes == ax[0]) or (event.inaxes == ax[1]):

            # SECOND PLOT - CUT OUT LC
            # plot the cut out around the time of the transit - changed by clicked on the image (either one of them)
            line.set_xdata(cutout(val)[0])
            line.set_ydata(cutout(val)[1])
            # ---------------
            # only ever plot one line so get rid of the old ones before plotting the new ones...
            if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

            # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
            lver0 = ax[0].axvline(val, color = 'r', linewidth = 2)
            lver1 = ax[1].axvline(val, color = 'r', linewidth = 2)

            # if clicked on the image to select the time to zoom in on then update the plot in this function.
            ax[1].set_xlim([val - 1,val + 1])
            lver0.set_xdata(val)
            lver1.set_xdata(val)

            # also update the location of the slider (aesthetic reasons only)
            # Define the slider to change the transit-event time (and cut out region) - this might slow the code down as it re-plots the slider everytime but I can't find an efficient way to do this differently for now.
            fig.canvas.draw_idle()  # update figure

            # these values need to be global as are needed but the clicking and the slider events
            global transit_slider_ax
            global transit_slider
            transit_slider_ax.remove() # remove the old bar and plot a new one which has the colour bar in the right place
            transit_slider_ax  = fig.add_axes([0.25, 0.14, 0.65, 0.03])  # location of the slider
            transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=val, color='teal')

        # if clicked outside of the region activate the slider bar
        else:
            transit_slider.on_changed(update_axis)
            scale_slider_upper.on_changed(update_yaxis_upper)
            scale_slider_lower.on_changed(update_yaxis_lower)
            transit_slider.on_changed(sliders_on_changed)

    # alternatively one can click on the image in order to change the zoom in region of the plot.
    fig.canvas.mpl_connect('button_press_event', onclick)

    # ---------------------

    # Define buttons which allow the user to interactively choose options:
    # run simple code, BLS, model the data usign Pyaneti, save the data and figures, creata DV report


    # only give the model option if pyaneti has been sucessfully installed
    if pyaneti_installed == True:
        var_ax = fig.add_axes([0.025, 0.44, 0.119, 0.21]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'model', 'Save', 'Report'), (False, args.noshow, args.north, False, False, True, True))
    else:
        var_ax = fig.add_axes([0.025, 0.45, 0.119, 0.2]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'Save', 'Report'), (False, args.noshow, args.north, False, True, True))

    # Initial values for each option
    simple = False
    hide = args.noshow
    north = args.north
    BLS = False
    model = False
    save = True
    DV = True

    # function to get the status of each button and to save it
    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        hide = not status[1]
        north = status[2]
        BLS = status[3]

        if pyaneti_installed == True:
            model = status[4]
            save = status[5]
            DV = status[6]
        else:
            model = False
            save = status[4]
            DV = status[5]

    # ---------------

    # define the paramaters of the plot
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf

    # -------------------------------------
    # specify plotting params
    ax[0].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[0].vlines(all_md, minf-1,minf + height*0.3 , colors = md_color, label = "Momentum Dump")

    ax[1].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("Time (BJD-2457000)", fontsize = main_plot_fontsize)
    ax[1].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[1].vlines(all_md, minf-1,minf + height*0.3, lw  = 1, colors = md_color, label = "Momentum Dump")
    # -------------------------------------

    # Create a button to close the figure and more onto the next stage of the code.
    ebx = plt.axes([0.77, 0.03, 0.13, 0.04])
    exit = Button(ebx, 'Done', color=main_pot_color)

    # define a box to enter the nickname of the target - this is useful if you analyse a lot of different candidates and you need a way of identifying them easily
    # ---------------
    # define a 'text box' that lets you enter the transit times

    initial_text = "FFI"  # no initial text (the text box is empty)
    nick_name = []  # list of the entered transit times

    # function to store the entererd transit times.
    def submit(text):
        nick_name.append(text)

    axbox = plt.axes([0.25, 0.03, 0.22, 0.04]) # x, y, width, height
    text_box = TextBox(axbox, 'Memorable name (optional): ', initial=initial_text)
    text_box.on_submit(submit)

    # ---------------
    # make button to exit the plot and continue with the code.
    # pop up warning if no transit time is entered
    def close(event):
        if len(transit_times) > 0:
            plt.close('all')
        else:
            print ("Must enter at least one transit time!")
            ax[1].text(0.4,-0.95, "Please enter at least one transit time!", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            fig.canvas.draw()


    exit.on_clicked(close)
    # ---------------

    global ttxt
    # text that will be updated to list the 'wanted' transit event times
    ttxt = plt.text(0,1.5, "Transit times: ", weight='bold')

    # Create a button to store the 'current' value of the slider
    stx = plt.axes([0.53, 0.03, 0.11, 0.04])
    store_val = Button(stx, 'Add time', color='gold')

    transit_times = []

    # function to store a transit event time using a button
    def storeval(time):

        storetime = (transit_slider.val)
        transit_times.append(round(storetime,2))
        global ttxt

        ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))

        plt.draw()


    store_val.on_clicked(storeval)
    # ---------------

    # Create a button to delete the last entered value
    delx = plt.axes([0.65, 0.03, 0.11, 0.04])
    del_val = Button(delx, 'Remove time', color='grey')

    transit_times = []

    # function to delete the last entry if using a button
    def deleval(time):
        global ttxt
        if len (transit_times)>0:
            del transit_times[-1]

            ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
            plt.draw()


    del_val.on_clicked(deleval)

    # ---------------

    plt.show()

    exit.on_clicked(close)


    # ---------------
    plt.show()

    # ------ END OF INTERACTIVE PART OF CODE ------

    # save the status of all the buttons
    end_status = save_var.get_status()

    if pyaneti_installed == True:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = end_status[4]
        save = end_status[5]
        DV = end_status[6]

    else:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = False
        save = end_status[4]
        DV = end_status[5]

    # ---------------

    # update the arguments if they are different from the ones entered in the command line:
    if args.noshow != hide:
        args.noshow = hide

    if args.north != north:
        args.north = north

    # if the save button was not selected, then change the iput argument to false
    if save == False:
        args.save = save

    # check whether the user identified a nickname, or memorable name, for the candidate

    if len(nick_name) == 0: # if the box was left untouched
        args.nickname = "FFI"

    elif len(nick_name[-1]) > 0: # [-1] in case the user pressed enter and then deleted the nickname again
        args.nickname = str(nick_name[-1])

    # ---------------


    # get the entered transit times and turn them into a list
    transit_list = transit_times # redefine the to note that this is a list

    print ("Transit times : {}".format(str(transit_list)[1:-1]))

    # if the save button was not selected, then change the iput argument to false
    if save == False:
        args.save = save



    # ----------------------------------------
    #               START BREWING ....
    # ----------------------------------------
    if args.FFI == 'SPOC':
        print ("20 minute cadence SPOC data available")
        # if the SPOC data is available, you can use run the normal report
        brew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, args)

    else:
        brew.brew_LATTE_FFI(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux_list, allflux_small, allflux, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, args)

'''

# the main interactive tool used to identify the times of the transit-like events when run in FFI mode
def interact_LATTE_FFI(tic, indir, syspath, sectors_all, sectors, ra, dec, args):
    '''
    Function to run the Interactive LATTE code for the FFI using the matplotlib interactive tool.
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
    ra  : float
        right ascension of the target
    dec   :
        declinatino of the target.

    *args
        if auto = False, the aperture will be chosen automatically using an intensity threshhold.

    '''

    # ----- download the data ------
    # by default you have to manually choose the different aperture sizes

    # see whether there is SPOC/QLP data
    global in_sec

    # ----- download the data ------
    # by default you have to manually choose the different aperture sizes

    # CHECK what ffi data products are available
    data_urls, args = check_FFI_data_products(tic,sectors, args)

    # see whether there is SPOC/QLP data
    if args.FFI == 'SPOC':
        alltime, allflux, allflux_err, all_md, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, url_list, args = download_data_FFI_SPOC(indir, sectors, syspath, sectors_all, tic, args, data_urls, save = False)

    elif args.FFI == 'QLP':

        alltime0, allflux0, allflux_small, allflux_err, all_md, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, X1_list, X4_list, bkg_list, t_list, tpf_list, tpf_filt_list,arrshape_list = download_data_FFI_QLP(indir, sectors, syspath, sectors_all, tic, args, data_urls, save = False)

        allflux_list = allflux0
        bkg_list = bkg_list
        apmask_list = [-111]
        allfbkg_t = alltime0

    # the above function returns [-111] if there is no SPOC or QLP data
    elif (args.FFI == True) and (args.auto == False):
        alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = interact_LATTE_FFI_aperture(tic, indir, sectors_all, sectors, ra, dec, args)

    # if the auto mode is selected in the command line, the aperture sizes are chosen automatically based on an intensity threshhold.
    elif (args.FFI == True) and (args.auto == True):
        print ("Start TPF data download.....", end =" ")
        alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t ,start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = download_data_FFI(indir, sectors, syspath, sectors_all, tic, args)
        print ("done.\n")


    # --------------------------------------------
    # Plot the interactive plot - uses matplolib -
    # --------------------------------------------

    # allflux is detrended
    # ------- median absolute deviation in order to determine the clipping
    print (args.FFI)
    print ("- - - - - - - - - - - - - - -")

    if (args.FFI != 'SPOC') and (args.FFI != 'QLP'):
        MAD = median_absolute_deviation(allflux0, ignore_nan = True)

        madrange = (50 * MAD)

        # make a mask for the points we want to get rid of
        ymask = (allflux0 < 1 + madrange) * (allflux0 > 1 - madrange)


        alltime = np.array(alltime0)[ymask]
        allflux = np.array(allflux0)[ymask]
        # -------------------------------------

        fluxmin = np.nanmin(allflux)
        fluxmax = np.nanmax(allflux)
        flux_diff = fluxmax - fluxmin
    

    if args.FFI == 'QLP':
        alltime = np.array(alltime0)
        allflux = np.array(allflux0)

        fluxmin = np.nanmin(allflux)
        fluxmax = np.nanmax(allflux)
        flux_diff = fluxmax - fluxmin

    # function to define the plotting area around the transit event.
    # this needs to be in a function as the area changes with the interactive slider.
    fig, ax = plt.subplots(2, 1, figsize=(11,7.5))
    plt.tight_layout()
    
    # Adjust the plots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.3)

    fluxmin = np.nanmin(allflux)
    fluxmax = np.nanmax(allflux)
    flux_diff = fluxmax - fluxmin
    # function to define the plotting area around the transit event.
    # this needs to be in a function as the area changes with the interactive slider.
    def cutout(transit):
        mask = (np.array(alltime) > transit-1) & (np.array(alltime) < transit+1)

        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux)]

    # ---------------

    # define the slider that let's you chose the sector to look at
    class PageSlider(matplotlib.widgets.Slider):

        def __init__(self, ax, label, numpages = 10, valinit=0, valfmt='%1d',
                     closedmin=True, closedmax=True,
                     dragging=True, **kwargs):

            self.facecolor=kwargs.get('facecolor',"w")
            self.activecolor = kwargs.pop('activecolor',"b")
            self.fontsize = kwargs.pop('fontsize', 10)
            self.numpages = numpages

            super(PageSlider, self).__init__(ax, label, 0, numpages,
                                valinit=valinit, valfmt=valfmt, **kwargs)

            # make it so that numpages is the index of the array of the
            global in_sec

            self.poly.set_visible(False)
            self.vline.set_visible(False)
            self.pageRects = []
            for i in range(numpages):

                facecolor = self.activecolor if i==valinit else self.facecolor
                r  = Rectangle((float(i)/numpages, 0), 1./numpages, 1,
                                    transform=ax.transAxes, facecolor=facecolor)
                ax.add_artist(r)
                self.pageRects.append(r)

                seclist_names = ['All'] + in_sec

                sec_i = seclist_names[i]

                ax.text(float(i)/numpages+0.5/numpages, 0.5, str(sec_i),
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=self.fontsize)
            self.valtext.set_visible(False)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            bax = divider.append_axes("right", size="5%", pad=0.05)
            fax = divider.append_axes("right", size="5%", pad=0.05)
            self.button_back = matplotlib.widgets.Button(bax, label=r'$\blacktriangleleft$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_forward = matplotlib.widgets.Button(fax, label=r'$\blacktriangleright$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_back.label.set_fontsize(self.fontsize)
            self.button_forward.label.set_fontsize(self.fontsize)
            self.button_back.on_clicked(self.backward)
            self.button_forward.on_clicked(self.forward)

        def _update(self, event):
            super(PageSlider, self)._update(event)
            i = int(self.val)
            if i >=self.valmax:
                return
            self._colorize(i)

        def _colorize(self, i):
            for j in range(self.numpages):
                self.pageRects[j].set_facecolor(self.facecolor)
            self.pageRects[i].set_facecolor(self.activecolor)

        def forward(self, event):
            current_i = int(self.val)
            i = current_i+1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

        def backward(self, event):
            current_i = int(self.val)
            i = current_i-1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

    # ---------------

    # define the initial conditions - this is what appears on the plot initially but can change later.
    transit = np.nanmean(alltime)  # start the cut-out in the centre of the plot.
    binfac = 7

    # FIRST PLOT - FULL LC
    # plot the unbinned flux
    [line_full] = ax[0].plot(alltime, allflux , marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3,  markerfacecolor='k')
    # ---------------

    # SECOND PLOT - CUT OUT LC
    # plot the cut out around the time of the transit - called with function as it can be changed with the slider.
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3,  markerfacecolor='k')
    # ---------------

    global transit_slider_ax
    global transit_slider

    if len(in_sec) == 1:

        # Define two sliders to change the y axis scale
        scale_slider_upper_ax  = fig.add_axes([0.24, 0.19, 0.33, 0.03], facecolor='silver')
        scale_slider_upper = Slider(scale_slider_upper_ax, 'Flux limits', 1, fluxmax + (flux_diff*0.1), valinit=fluxmax, color='white')
        scale_slider_upper.valtext.set_visible(False)

        scale_slider_lower_ax  = fig.add_axes([0.57, 0.19, 0.33, 0.03], facecolor='white') # x, y, width, height
        scale_slider_lower = Slider(scale_slider_lower_ax, ' ', fluxmin - (flux_diff*0.1), 1,valinit=fluxmin, color='silver')
        scale_slider_lower.valtext.set_visible(False)

        # Define the slider to change the transit-event time (and cut out region)
        transit_slider_ax  = fig.add_axes([0.24, 0.15, 0.66, 0.03])  # location of the slider
        transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=transit, color='teal')

    # if there it's in more than one sector there is an extra slide bar to chose to look at individual sectors.
    else:
        # Define the slider to change the y axis scale
        scale_slider_upper_ax  = fig.add_axes([0.24, 0.2, 0.33, 0.02], facecolor='silver')
        scale_slider_upper = Slider(scale_slider_upper_ax, 'Flux limits', 1, fluxmax + (flux_diff*0.1), valinit=fluxmax, color='white')
        scale_slider_upper.valtext.set_visible(False)

        scale_slider_lower_ax  = fig.add_axes([0.57, 0.2, 0.33, 0.02], facecolor='white')
        scale_slider_lower = Slider(scale_slider_lower_ax, ' ', fluxmin - (flux_diff*0.1), 1,valinit=fluxmin, color='silver')
        scale_slider_lower.valtext.set_visible(False)

        # Define the slider to change the transit-event time (and cut out region)
        transit_slider_ax  = fig.add_axes([0.24, 0.17, 0.66, 0.02])  # location of the slider
        transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=transit, color='teal')

        # define the slider to change the plotting region
        sector_slider_ax  = fig.add_axes([0.24, 0.11, 0.66, 0.03])
        sector_slider = PageSlider(sector_slider_ax, 'Sector', len(in_sec) + 1, activecolor=main_pot_color)

    # ---------------

    # deifne the intial x and y axis limits - y limit can be changed with slider, x limit cannot.
    ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])

    ax[1].set_xlim([np.nanmean(alltime)-1, np.nanmean(alltime)+1])
    ax[1].set_ylim([fluxmin, fluxmax])

    # ---------------
    # Define an action for acessing the required cut-out data and drawing the line when the slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])

        fig.canvas.draw_idle()  # update figure

    # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)

    # ---------------
    # Define an action for modifying the plot region on the second plot
    def update_axis(val):
        # ---------------
        # only ever plot one line so get rid of the old ones before plotting the new ones...
        if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

        # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
        lver0 = ax[0].axvline(transit_slider.val, color = 'r', linewidth = 2)
        lver1 = ax[1].axvline(transit_slider.val, color = 'r', linewidth = 2)

        ax[1].set_xlim([transit_slider.val - 1,transit_slider.val + 1])
        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)

    def update_plotting_region(val):
        global in_sec
        seclist_names = ['All'] + in_sec

        index = int(val)

        if index == 0:
            ax[0].set_title("TIC {}     Sectors {}".format(tic,str(in_sec)[1:-1]))
            ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
        else:
            ax[0].set_title("TIC {}   Sector {}".format(tic,seclist_names[index]))
            ax[0].set_xlim([start_sec[index-1][0],end_sec[index-1][0]])


    def update_yaxis_upper(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])

    def update_yaxis_lower(val):

        upper_axis_lim = scale_slider_upper.val
        lower_axis_lim = scale_slider_lower.val

        ax[0].set_ylim([lower_axis_lim, upper_axis_lim])
        ax[1].set_ylim([lower_axis_lim, upper_axis_lim])


    # ------ ON CLICK -------
    # define the action that lets you click on the image in order to chose the location instead of using the slider - probably more useful
    def onclick(event):

        # just need to know the x (time) position of the click. The y (flux) position doesn't matter.
        val = event.xdata

        # check that the click was within the plotting region. If not ignore the click
        if (event.inaxes == ax[0]) or (event.inaxes == ax[1]):

            # SECOND PLOT - CUT OUT LC
            # plot the cut out around the time of the transit - changed by clicked on the image (either one of them)
            line.set_xdata(cutout(val)[0])
            line.set_ydata(cutout(val)[1])

            # ---------------
            # only ever plot one line so get rid of the old ones before plotting the new ones...
            if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

            # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
            lver0 = ax[0].axvline(val, color = 'r', linewidth = 2)
            lver1 = ax[1].axvline(val, color = 'r', linewidth = 2)

            # if clicked on the image to select the time to zoom in on then update the plot in this function.
            ax[1].set_xlim([val - 1,val + 1])
            lver0.set_xdata(val)
            lver1.set_xdata(val)

            # also update the location of the slider (aesthetic reasons only)
            # Define the slider to change the transit-event time (and cut out region) - this might slow the code down as it re-plots the slider everytime but I can't find an efficient way to do this differently for now.
            fig.canvas.draw_idle()  # update figure

            # these values need to be global as are needed but the clicking and the slider events
            global transit_slider_ax
            global transit_slider

            transit_slider_ax.remove() # remove the old bar and plot a new one which has the colour bar in the right place

            if len(in_sec) == 1:
                transit_slider_ax  = fig.add_axes([0.24, 0.15, 0.66, 0.03])  # location of the slider
                transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=val, color='teal')
            else:
                transit_slider_ax  = fig.add_axes([0.24, 0.17, 0.66, 0.02])  # location of the slider
                transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltime), np.nanmax(alltime), valinit=val, color='teal')


        # if clicked outside of the region activate the slider bar
        else:
            # Link the functions (actions) to the corresponding sliders
            transit_slider.on_changed(update_axis)
            scale_slider_upper.on_changed(update_yaxis_upper)
            scale_slider_lower.on_changed(update_yaxis_lower)
            transit_slider.on_changed(sliders_on_changed)

    # alternatively one can click on the image in order to change the zoom in region of the plot.
    fig.canvas.mpl_connect('button_press_event', onclick)

    # ---------------------

    # Define buttons which allow the user to interactively choose options:
    # run simple code, BLS, model the data usign Pyaneti, save the data and figures, creata DV report

    # only give the model option if pyaneti has been sucessfully installed
    if pyaneti_installed == True:
        var_ax = fig.add_axes([0.025, 0.52, 0.119, 0.19]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'model', 'Save', 'Report'), (False, args.noshow, args.north, False, False, True, True))
    else:
        var_ax = fig.add_axes([0.025, 0.53, 0.119, 0.18]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Show plots', 'North', 'BLS', 'Save', 'Report'), (False, args.noshow, args.north, False, True, True))

    if len(in_sec) > 1:
        sector_slider.on_changed(update_plotting_region)


    '''
    aster_ax = fig.add_axes([0.025, 0.23, 0.119, 0.045], facecolor='gold') # x, y, width, height

    asteroseismic_button = Button(aster_ax, 'Seismology', color = 'white')

    asteroseismic_colors = itertools.cycle(['gold','white'])
    asteroseismic_status = itertools.cycle([True, False])

    global asteroseismology
    asteroseismology = False

    def change_color(event):
        global asteroseismology
        asteroseismic_button.color = next(asteroseismic_colors)

        asteroseismology = next(asteroseismic_status)

        # If you want the button's color to change as soon as it's clicked, you'll
        # need to set the hovercolor, as well, as the mouse is still over it
        asteroseismic_button.hovercolor = asteroseismic_button.color
        fig.canvas.draw()

    asteroseismic_button.on_clicked(change_color)
    '''

    # Initial values for each option
    simple = False
    hide = args.noshow
    north = args.north
    BLS = False
    model = False
    save = True
    DV = True

    # function to get the status of each button and to save it
    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        hide = not status[1]
        north = status[2]
        BLS = status[3]

        if pyaneti_installed == True:
            model = status[4]
            save = status[5]
            DV = status[6]
        else:
            model = False
            save = status[4]
            DV = status[5]

    # ---------------

    # ---------------

    # define the paramaters of the plot
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf

    ax[0].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[0].vlines(all_md, minf-1,minf + height*0.3 , colors = md_color, label = "Momentum Dump")

    ax[1].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("Time (BJD-2457000)", fontsize = main_plot_fontsize)
    ax[1].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[1].vlines(all_md, minf-1,minf + height*0.3, lw =1,  colors = md_color, label = "Momentum Dump")


    # ---------------
    # Create a label for the 'settings' box to clarify what it does

    plt.text(0.25, 1.07, "Settings", fontsize=10, verticalalignment='center')

    # ---------------
    # Create a label for the 'quick look options'

    plt.text(0.22, -0.35, "Quick look", fontsize=10, verticalalignment='center')

    # define a box to enter the nickname of the target - this is useful if you analyse a lot of different candidates and you need a way of identifying them easily

    # ---------------
    # define a 'text box' that lets you enter the transit times

    initial_text = ""  # no initial text (the text box is empty)
    nick_name = []  # list of the entered transit times

    # function to store the entererd transit times.
    def submit(text):
        nick_name.append(text)

    axbox = plt.axes([0.24, 0.01, 0.22, 0.04]) # x, y, width, height
    text_box = TextBox(axbox, 'Memorable name (optional): ', initial=initial_text)
    text_box.on_submit(submit)
    # ---------------


    # Create a button to close the figure and more onto the next stage of the code.
    ebx = plt.axes([0.77, 0.01, 0.13, 0.04])
    exit = Button(ebx, 'Done', color=main_pot_color)

    # make button to exit the plot and continue with the code.
    # pop up warning if no transit time is entered

    def close(event):

        if len(transit_times) > 0:
            plt.close('all')

        else:
            print ("Must enter at least one transit time!")
            global deltxt

            if len(in_sec) > 1:
                deltxt = ax[1].text(3.7,-10.5,"Please enter at least one transit time!", color = 'red', size=9)
            else:
                deltxt = ax[1].text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)

            plt.draw()


    exit.on_clicked(close)
    # ---------------

    global ttxt
    global deltxt

    if len(in_sec) > 1:
        deltxt = plt.text(3.7,-10.5," ", color = 'red', size=9)
    else:
        deltxt = plt.text(3.7,-10," ", color = 'red', size=9)


    # text that will be updated to list the 'wanted' transit event times
    if len(in_sec) > 1:
        ttxt = plt.text(-4.88,1.5, "Transit times: ", weight='bold')
    else:
        ttxt = plt.text(-4.88,2, "Transit times: ", weight='bold')

    # Create a button to store the 'current' value of the slider
    stx = plt.axes([0.53, 0.01, 0.11, 0.04])
    store_val = Button(stx, 'Add time', color='gold')

    transit_times = []

    # function to store a transit event time using a button
    def storeval(time):

        storetime = (transit_slider.val)
        transit_times.append(round(storetime,2))
        global ttxt
        global deltxt

        ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
        deltxt.set_text(" ")

        plt.draw()


    store_val.on_clicked(storeval)
    # ---------------

    # Create a button to delete the last entered value
    delx = plt.axes([0.65, 0.01, 0.11, 0.04])
    del_val = Button(delx, 'Remove time', color='grey')

    transit_times = []

    # function to delete the last entry if using a button
    def deleval(time):
        global ttxt
        if len (transit_times)>0:
            del transit_times[-1]

            ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
            plt.draw()


    del_val.on_clicked(deleval)

    # ---------------

    # Make some quick look buttons to look at certain plots without having to run the program/make a report

    if args.FFI == 'SPOC':
        plotcentx = plt.axes([0.025, 0.45, 0.12, 0.025]) # x,y, xwidth, ywidth
        plotcent = Button(plotcentx, 'Centroids', color='white')
        plotcent.label.set_fontsize(10)

    plotbgkx = plt.axes([0.025, 0.415, 0.12, 0.025]) # x,y, xwidth, ywidth
    plotbgk = Button(plotbgkx, 'Background', color='white')
    plotbgk.label.set_fontsize(10)

    ploteepx = plt.axes([0.025, 0.38, 0.12, 0.025]) # x,y, xwidth, ywidth
    ploteep = Button(ploteepx, 'Stellar tracks', color='white')
    ploteep.label.set_fontsize(10)

    plotperiodx = plt.axes([0.025, 0.345, 0.12, 0.025]) # x,y, xwidth, ywidth
    plotperiod = Button(plotperiodx, 'Periodogram', color='white')
    plotperiod.label.set_fontsize(10)

    # ------
    # these require a different data download of the TPFs

    #plotapx = plt.axes([0.025, 0.3, 0.12, 0.03]) # x,y, xwidth, ywidth
    #plotap = Button(plotapx, 'Ap. size', color='gold')

    #plotpixlevelx = plt.axes([0.025, 0.25, 0.12, 0.03]) # x,y, xwidth, ywidth
    #plotpixlevel = Button(plotpixlevelx, 'Pixel square', color='gold')

    # - - - - - - - - - -
    def quick_look_bkg(time):
        global deltxt
        if len (transit_times) > 0:

            plot_background(tic, indir, alltime, allfbkg, transit_times, args, ql = True)
        else:
            deltxt = plt.text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)
            print ("select at least one transit time to run this")

    plotbgk.on_clicked(quick_look_bkg)

    # - - - - - - - - - -
    def quick_look_centroid(time):
        global deltxt
        if len (transit_times) > 0:

            plot_centroid(tic, indir, alltime12, allx1, ally1, allx2, ally2, transit_times, args, ql = True)
        else:
            deltxt = plt.text(3.7,-10,"Please enter at least one transit time!", color = 'red', size=9)
            print ("select at least one transit time to run this")

    if args.FFI == 'SPOC':
        plotcent.on_clicked(quick_look_centroid)

    # - - - - - - - - - -
    def quick_look_eep(time):

        eep_target(tic, indir, syspath,teff, srad, args, ql = True)

    ploteep.on_clicked(quick_look_eep)

    # - - - - - - - - - -

    def quick_look_periodogram(time):

        plot_periodogram(tic, indir, alltime, allflux, teff, 0, 0, args, ql = True)

    plotperiod.on_clicked(quick_look_periodogram)

    # - - - - - - - - - -

    plt.show()
    plt.close('all')

    # ------ END OF INTERACTIVE PART OF CODE ------

    # save the status of all the buttons
    end_status = save_var.get_status()

    if pyaneti_installed == True:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = end_status[4]
        save = end_status[5]
        DV = end_status[6]

    else:
        simple = end_status[0]
        hide = not end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = False
        save = end_status[4]
        DV = end_status[5]


    # ---------------

    # update the arguments if they are different from the ones entered in the command line:
    if args.noshow != hide:
        args.noshow = hide

    if args.north != north:
        args.north = north

    # if the save button was not selected, then change the iput argument to false
    if save == False:
        args.save = save

    # check whether the user identified a nickname, or memorable name, for the candidate
    if len(nick_name) != 0: # if the box was left untouched do nothing and continue without setting a new nickname - will be ignored
        if len(nick_name[-1]) > 0: # [-1] in case the user pressed enter and then deleted the nickname again
            args.nickname = str(nick_name[-1])

    # ---------------

    # get the entered transit times and make sure that all the values are floats
    transit_list = transit_times

    print ("Transit times : {} \n".format(str(transit_list)[1:-1]))

    # get the final upper and lower axis limit - used for later plotting
    upper_axis_lim_final = scale_slider_upper.val
    lower_axis_lim_final = scale_slider_lower.val

    #if asteroseismology == True:
    #    interact_LATTE_asteroseismology(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, args)
    #else:


    # ----------------------------------------
    #               START BREWING ....
    # ----------------------------------------

    if args.FFI == 'SPOC':
        print ("30 minute cadence SPOC data available")
        # if the SPOC data is available we can run the 'normal' LATTE mode as they have all the same data products as the 2 minute SPOC data products
        brew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltime, allflux, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, upper_axis_lim_final, lower_axis_lim_final, tessmag, teff, srad, ra, dec, -999, -999, url_list, args)

    else:
        if args.FFI == 'QLP':
            # if QLP data is available we have to extract the TPF data now, so that we can run the FFI
            print ("30 minute cadence QLP data available")


        brew.brew_LATTE_FFI(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, upper_axis_lim_final, lower_axis_lim_final, args)


# --------------------------------------------
#         Download the data access files      #
# --------------------------------------------

# The Functions needed to get the files that know how to acess the data
# These text files contain the curl scripts needed to download the TESS data (both FFIs and SC data)
# The files (~ 325M) are stored on your computer in order to make the download of the data faster when LATET is run (1 bulk download vs downloading these curl scipts every time)

def data_files(indir):
    '''
    NO LONGER USED IN NEW VERSION (1.0) OF LATTE
    Function to download all of the text files that we need to the local computer.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    '''

    if not os.path.exists("{}/data/tesscurl_sector_all_lc.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_lc.sh".format(indir),'w') as f:
            f.write("#all LC file links")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re-run to download new data)
        print ("Download all required text files for available sectors starting with sector 1")

    else: # if the file already exists check whether there is something in the file
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_lc.sh > {0}/data/temp.txt'.format(indir))
        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]

        if string == "#all LC file links": # if this is the last (and only) line
            first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re-run to download new data)
            print ("Download all required text files for available sectors starting with sector 1")
        else:
            first_sec = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here



    if not os.path.exists("{}/data/tesscurl_sector_all_tp.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_tp.sh".format(indir),'w') as f:
            f.write("#all LC file links") # give some information at the top of the saved file for later reference (the # means this line will be inored when read back later)
        first_sec_tp = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re-run to download new data)
        print ("Download all required text files for available sectors starting with sector 1")

    else:
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_tp.sh > {0}/data/temp_tp.txt'.format(indir))

        with open("{}/data/temp_tp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]

        if string == "#all LC file links": # if this is the last (and only) line
            first_sec_tp = 0
            print ("Download all required text files for available sectors starting with sector 1")
        else:
            first_sec_tp = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here


    for sec in range(first_sec+1,500): # 26 because that's how many TESS sectors there will be in total

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
                print("Finished adding text files for sector {}".format(sec))
                # write the contents of the response (r.content)
                # to a new file in binary mode.

        with open("{}/data/tesscurl_sector_{}_lc.sh".format(indir,sec), 'wb') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_LC.content)


    for sec in range(first_sec_tp+1,500): # 26 because that's how many TESS sectors there will be in total

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

def tp_files(indir):
    '''
    Function to download all of the TPF data that we want to the local computer.
    These are needed for the nearest neighbour analysis.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    '''

    if not os.path.exists("{}/data/all_targets_list.txt".format(indir)):

        empty_dataframe = pd.DataFrame(columns=['TICID', 'Camera', 'CCD', 'Tmag', 'RA', 'Dec', 'sec'])

        empty_dataframe.to_csv("{}/data/all_targets_list.txt".format(indir), index=False)

        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("We also need to have a record of all of the TIC ID available in each Sectors. \
                This will only take a minute.  Starting with Sector 1")

    else:
        files = np.sort(glob('{}/data/all_targets_S*'.format(indir)))

        exist = []
        for f in files:
            exist.append(int(f[-9:-7]))  # get the first sector number (last that has already been downloaded)

        first_sec = (np.max(exist))


    for sec in range(first_sec+1,500):  # how many tess sectors will there be in total?

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

        # need to add what sector it's from
        # open with pandas and save it again with the sectro column - there must be a better way to do this but I can't think of it now.
        this_sector_file = pd.read_csv("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), comment = '#', delimiter = '\t', names = ['TICID', 'Camera', 'CCD', 'Tmag', 'RA', 'Dec'])
        this_sector_file['sec'] = sec
        this_sector_file.to_csv("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), index=False)

        # also open the latest version of the 'all' file so that we can add the new part to it... again not the best way
        all_sectors_file = pd.read_csv("{}/data/all_targets_list.txt".format(indir, download_sector), comment = '#', delimiter = ',')

        # merge the dataframes
        full_data_frame_allsectors = pd.concat([all_sectors_file, this_sector_file])

        # save it again
        full_data_frame_allsectors.to_csv("{}/data/all_targets_list.txt".format(indir, download_sector), index=False)

        print("Finished adding target list text file for sector {}".format(sec))

def TOI_TCE_files(indir):
    '''
    Function to download the files that list all the known TOI's and TCEs.
    This is useful to display on the DV report.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    '''

    # ------ TOIs ------
    TOI_url = "https://tev.mit.edu/data/collection/193/csv/5/"

    r_TOI = requests.get(TOI_url) # create HTTP response object

    if r_TOI.status_code == 404:
        print ("Can't download the TOI list at the moment. Has the URL changed? If this persists contact Nora.")

    with open("{}/data/TOI_list.txt".format(indir),'wb') as f:
           f.write(r_TOI.content)

    # this TOI file is needed for the TOI stellar paramaters
    TOI_star_params_rls = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv'

    r_TOI_star_params = requests.get(TOI_star_params_rls) # create HTTP response object

    if r_TOI_star_params.status_code == 404:
        print ("Can't download the TOI list at the moment. Has the URL changed? If this persists contact Nora.")

    with open("{}/data/TOI_list_star_params.txt".format(indir),'wb') as f:
           f.write(r_TOI_star_params.content)

    # ------ TCE ------
    if not os.path.exists("{}/data/tesscurl_sector_all_dv.sh".format(indir)):
        with open("{}/data/tesscurl_sector_all_dv.sh".format(indir),'w') as f:
            f.write("#all LC file links")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ('- - - - -')
        print ("Finally we need to download the text files that tell us which targets already have DV (and where we can find them).\
            This will only take a minute. Startign with Sector 1")

    else:
        os.system('tail -n 1 {0}/data/tesscurl_sector_all_dv.sh > {0}/data/temp.txt'.format(indir))

        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]

        first_sec = int(string.split('-')[5][2:]) # this is the last imported sector 0 start from here


    for sec in range(first_sec+1,500): # 52 because that's how many TESS sectors there will be in the first 4 years

        TCE_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_dv.sh".format(sec)
        r_TCE = requests.get(TCE_url) # create HTTP response object

        if r_TCE.status_code == 404:
            print ("You're all caught up on your list of TCE and TOIS up to Sector {} -- try downloading more data later".format(sec-1))
            break

        with open("{}/data/tesscurl_sector_all_dv.sh".format(indir), 'ab') as f:
                '''
                Saving recieved content as a png file in binary format
                '''
                f.write(r_TCE.content)
                print("Finished adding DV links for sector {}".format(sec))

def get_data_codes(indir):

    if not os.path.exists("{}/data/sector_download_codes.txt".format(indir)):
        with open("{}/data/sector_download_codes.txt".format(indir),'w') as f:
            f.write("#all LC file links \n")
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re-run to download new data)
        print ("We need some 'codes' specific to each sector to download the data later. \
            Getting these codes and saving them to a text file will only take a minute. Starting with Sector 1")

    else: # if the file already exists check whether there is something in the file
        os.system('tail -n 1 {0}/data/sector_download_codes.txt > {0}/data/temp.txt'.format(indir))
        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()[-1]

        if string == "#all LC file links": # if this is the last (and only) line
            first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re-run to download new data)
            print ("We need some 'codes' specific to each sector to download the data later. \
                Getting these codes and saving them to a text file will only take a minute. Starting with Sector 1")
        else:

            first_sec = int(str(string).split(' ')[0]) # this is the last imported sector 0 start from here


    for sec in range(first_sec+1,500): # 200 - just a high number of (how many sectors will there be in TESS?)

        LC_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh".format(sec)
        r_LC = requests.get(LC_url) # create HTTP response object

        if r_LC.status_code == 404:
            print ("You're all caught up with your codes up to Sector {} -- try downloading more data later".format(sec-1))
            break


        with open("{}/data/sector_download_codes.txt".format(indir), 'a') as f:
                '''
                open the first line of the file in order to get the codes that are requred to download data from that sector
                '''
                code1 = (str(r_LC.content[0:200]).split('-')[4][6:])
                code2 = (str(r_LC.content[0:200]).split('-')[7])

                string = "{} {} {}\n".format(str(sec), code1, code2 )

                # write the contents of the response (r.content)
                # to a new file in binary mode.

                f.write(string)
                print("Finished adding sector codes for sector {}".format(sec))

def momentum_dumps_info(indir):
    '''
    function to the create a list of all of the momentum dump times for each sector - only needed in the FFIs.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    '''
    print (" - - - - - - - -")
    print ("We need to store the times of the momentum dumps for each sector. \
        this will quickly open one file from each sector (no files are stored). \
        This will take a minute.")
    print()

    if not os.path.exists("{}/data/tess_mom_dumps.txt".format(indir)):
        with open("{}/data/tess_mom_dumps.txt".format(indir),'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['sec', 'time'])
        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("Will determine the times of the momentum dumps for all sectors and store them... ")

    else:
        os.system('tail -n 1 {0}/data/tess_mom_dumps.txt > {0}/data/temp.txt'.format(indir))

        with open("{}/data/temp.txt".format(indir), 'r') as f:
            string = f.readlines()#[-1]
            if string == ['sec\ttime\n']:
                first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
            else:
                first_sec = int(string[0][0:2])


    sectorfile = "{}/data/sector_download_codes.txt".format(indir)
    infile = pd.read_csv(sectorfile, delimiter = ' ', usecols=[0,1,2],names = ['sec', 'first', 'second'], comment = '#')

    for sec in range(first_sec+1,500): # a large number - how many TESS sectors will there be?

        try:
            if sec < 10:
                download_sector = "00{}".format(sec)
            else:
                download_sector = "0{}".format(sec)

            # load the data file of the sector that the first marked transit appears in
            # sort the list to be ordered by camera, RA and then Dec
            tic_list = pd.read_csv("{}/data/all_targets_S{}_v1.txt".format(indir,download_sector))

            # select the first tic listed for that sector - the momentum dumps are the same for all tic ids - they are sector dependent.
            tic =  list(tic_list['TICID'])[0]

            this_sector_code = infile[infile.sec == int(sec)]

            sector = str(sec)

            download_url = (
                "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
                + str(this_sector_code['first'].values[0]).rjust(13, "0")
                + "-s"
                + sector.rjust(4, "0")
                + "-"
                + str(tic).rjust(16, "0")
                + "-"
                + str(this_sector_code['second'].values[0]).rjust(4, "0")
                + "-s_lc.fits")

            response = requests.get(download_url)
            lchdu  = pf.open(response.url) # This needs to be a URL - not a file

            # Open and view columns in lightcurve extension
            lcdata = lchdu[1].data
            quality = lcdata['QUALITY']
            time    = lcdata['TIME']

            lchdu.close()

            mom_dump_mask = np.bitwise_and(quality, 2**5) >= 1

            momdump = (list(time[mom_dump_mask]))
            sector = list([sec] * len(momdump))

            with open("{}/data/tess_mom_dumps.txt".format(indir), 'a') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(zip(sector,momdump))

        except:
            print ("Done")
            break

# -----------------------

def tess_point(indir,tic):
    '''
    THIS FUNCTION IS NO LONGER USED IN LATTE VERSION 1.0 AND BEYOND

    Use tess-point to find out in what sectors the dubject appears in.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    tic : str
        TIC (Tess Input Catalog) ID of the target

    Returns
    -------
    outSec  : list
        list of all the sectors in which this subejct has been /  will be observed

    '''

    if not exists('{}/tesspoint'.format(indir)):
        os.makedirs('{}/tesspoint'.format(indir))

    os.system('python3 -m tess_stars2px -t {} > {}/tesspoint/{}_tesspoint.txt'.format(tic,indir,tic))

    df = pd.read_csv('{}/tesspoint/{}_tesspoint.txt'.format(indir,tic), comment = '#', delimiter = '|', names = ['TIC','RA','Dec','EclipticLong','EclipticLat','Sector','Camera','Ccd','ColPix', 'RowPix'])

    return list(df['Sector']), float(df['RA'][0]), float(df['Dec'][0])


def transit_sec(in_sec,start_sec, end_sec, transit_list):
    '''
    determine in what sectors the marked transit are.

    Parameters
    ----------
    in_sec : list
        the sectors that we are looking at.
    start_sec : list
        the start time of these sectors
    end_sec : list
        the end time of these sectors
    transit_list : list
        list of the marked peaks

    Returns
    -------
    transit_sec  : list
        list of the sectors in which the peaks appear.

    '''

    transit_sec = []
    for peak in transit_list:
        for n,start in enumerate(start_sec):

            if peak >= start[0] and peak <= end_sec[n][0]:

                transit_sec.append(in_sec[n])

    return list(set(transit_sec))


def nn_ticids(indir, transit_sec, tic):
    '''
    function to find the TIC IDs of the 6 nearest neighbour stars that were observed by TESS (short cadence).

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_sec  : list
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

    neighbours_sector = transit_sec[0]

    if neighbours_sector < 10:
        download_sector = "00{}".format(neighbours_sector)
    else:
        download_sector = "0{}".format(neighbours_sector)

    # load the data file of the sector that the first marked transit appears in
    # sort the list to be ordered by camera, RA and then Dec
    tic_list = pd.read_csv("{}/data/all_targets_S{}_v1.txt".format(indir,download_sector)).sort_values(['Camera', 'RA', 'Dec']).reset_index()

    # acess the information for the target stars
    target = tic_list.loc[tic_list['TICID'] == float(tic)]
    tic_idx = target.index[0]  # get the index of the target in the list

    # get the RA and Dec of the target star - this will be returned but the function
    target_ra = float(target['RA'])
    target_dec = float(target['Dec'])

    # ------------

    # make a list of the closest stars to the target (only stars that are TESS target stars are considered)
    # ensure that the TIC ID is not at the end or begining of the tic list - otherwise we can't take 100 from wither side.
    if tic_idx < 101:
        tic_list_close = tic_list[0:tic_idx + 101]
    elif tic_idx > (len(tic_list) + 1):
        tic_list_close = tic_list[tic_idx - 100 :len(tic_list)]
    else:
        tic_list_close = tic_list[tic_idx - 100:tic_idx + 101]
    # ------------

    # function to alculated the angular separation from the target to the stars to find the nearest neighbours
    def star_sep(row, ra, dec):

        ra2 = float(row['RA'])
        dec2 = float(row['Dec'])

        c1 = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
        c2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')
        sep = c1.separation(c2)

        return sep.arcminute

    # make a new column in the pandas dataframe of the separation
    tic_list_close['dist'] = tic_list_close.apply(star_sep, args = (target_ra,target_dec), axis=1)

    # Select the 6 nearest neightbours and creare a list of their tic IDs and distance to them.
    closest_tic = tic_list_close.sort_values('dist')[0:6]
    ticids = closest_tic['TICID'].tolist()
    distance = closest_tic['dist'].tolist()

    return ticids, distance, target_ra, target_dec


# --------------------------------------------
#              download the data             #
# --------------------------------------------

# The functions to download the actual data
def download_data(indir,sectors, tic, binfac = 5, test = 'no'):

    '''
    Download the LCs for the target star for all the indicated sectors

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    sector  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    tic : str
        TIC (Tess Input Catalog) ID of the target
    binfac  :  int
        The factor by which the data should be binned. Default = 5 (which is what is shown on PHT)

    test   :   str
        in order to test the function with unittests we want to run it with an input file (string to input file)

    Returns
    -------
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    allflux_err  :  list
        normalized flux errors (not binned)
    all_md  :  list
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
    alltime12  :  list
        time used for the x and y centroid position plottin
    allfbkg  :  list
        background flux
    start_sec  :  list
        times of the start of the sector
    end_sec  :  list
        times of the end of the sector
    in_sec  :  list
        the sectors for which data was downloaded
    tessmag  :  float
        TESS magnitude of the target star
    teff  :  float
        effective temperature of the tagret star (K)
    srad  :  float
        radius of the target star (solar radii)

    '''

    # plot using the seaborn library
    sb.set(style='ticks')
    sb.set_color_codes('deep')

    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)


    # -!-!-!-!-!-!-!-
    if test != 'no':
        dwload_link = [test]  # in the unittest define a link to a pre downloaded file
    # -!-!-!-!-!-!-!-

    # if not a test find the link needed to access the data
    else:

        dwload_link = []

        sectorfile = "{}/data/sector_download_codes.txt".format(indir)

        infile = pd.read_csv(sectorfile, delimiter = ' ', names = ['sec', 'first', 'second'], comment = '#')


        lastavailablesector = list(infile['sec'])[-1]
        sectors = np.array(sectors)[np.array(sectors) <= int(lastavailablesector)]

        # search the LC download file for the URL to download the LC from MAST for the given target star.
        # if we are looking at all of the sectors search, the file that has all of the LC URLs in it.
        if sectors == 'all':
            for sector in sectors:
                try:
                    this_sector_code = infile[infile.sec == int(sector)]

                    sector = str(sector)

                    download_url = (
                        "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
                        + str(this_sector_code['first'].values[0]).rjust(13, "0")
                        + "-s"
                        + sector.rjust(4, "0")
                        + "-"
                        + str(tic).rjust(16, "0")
                        + "-"
                        + str(this_sector_code['second'].values[0]).rjust(4, "0")
                        + "-s_lc.fits")

                    dwload_link.append(download_url)
                except:
                    future_sectors.append(sector)
        # otherwise only load the files specific for each sector that we are looking at - this is to avoid loading the larger file.
        else:
            future_sectors = []
            # locate the file string
            # for s in sector:

            for sector in sectors:
                sector = str(sector)
                try:
                    this_sector_code = infile[infile.sec == int(sector)]

                    sector = str(sector)

                    download_url = (
                        "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
                        + str(this_sector_code['first'].values[0]).rjust(13, "0")
                        + "-s"
                        + sector.rjust(4, "0")
                        + "-"
                        + str(tic).rjust(16, "0")
                        + "-"
                        + str(this_sector_code['second'].values[0]).rjust(4, "0")
                        + "-s_lc.fits")

                    dwload_link.append(download_url)

                except:
                    future_sectors.append(sector)


            if len(future_sectors):
                print ("TIC {} has not (yet) been observed in sectors {} as part of the SC campaign (check the FFIs).".format(tic, future_sectors))

        if len(dwload_link) == 0:
            print ("TIC {} was not observed in Sector(s):   {} as part of the SC data. Try again with different sectors or check whetehr it's a n FFI instead.".format(tic, sector))

            print ("\n (Also check that the data files have been downloaded - run code with '--new-data' in the command line)")
            raise SystemExit

    # define all the empty lists to append to in order to return the data that will be requrides later on in the script

    alltimebinned = []
    allfluxbinned = []

    allx1 = []
    allx2 = []
    ally1 = []
    ally2 = []
    alltime12 = []
    allfbkg = []

    start_sec = []
    end_sec = []
    in_sec = []

    alltime = []
    allflux = []
    allflux_err = []
    all_md = []

    failed_sectors = []
    # loop through all the download links - all the data that we want to access
    print ("download_link")

    for lcfile in dwload_link:

        # !-!-!-!-!-!-!-
        # if this a test run, download the file already on the system
        if test != 'no':
            lchdu  = pf.open(lcfile)
        # !-!-!-!-!-!-!-

        else:
            try:
                # use the downlload link to download the file from the server - need an internet connection for this to work
                response = requests.get(lcfile)
                # open the file using the response url
                lchdu  = pf.open(response.url) # this needs to be a URL - not a file
            except:
                failed_sectors.append(int(lcfile.split('-')[1][1:]))
                continue

        # open and view columns in lightcurve extension
        lcdata = lchdu[1].data


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

        lchdu.close()

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
        alltime12.append(list(time[l2]))

        allfbkg.append(fbkg)

        start_sec.append([time[0]])
        end_sec.append([time[-1]])

    # flatten lists of list

    alltime =       list(np.hstack(alltime)      )
    allflux =       list(np.hstack(allflux)      )
    allflux_err =   list(np.hstack(allflux_err)  )
    all_md =        list(np.hstack(all_md)       )

    alltimebinned = list(np.hstack(alltimebinned))
    allfluxbinned = list(np.hstack(allfluxbinned))

    allx1 =         list(np.hstack(allx1)     )
    allx2 =         list(np.hstack(allx2)     )
    ally1 =         list(np.hstack(ally1)     )
    ally2 =         list(np.hstack(ally2)     )
    alltime12 =     list(np.hstack(alltime12) )
    allfbkg =       list(np.hstack(allfbkg)   )

    if len (failed_sectors) != 0:
        print ("Failed sectors: {}    (These may be future sectors of may have only been observed in the FFIs".format(failed_sectors))

    return alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad


def download_data_FFI_interact(indir,sector, sectors_all, tic, save = False):

    '''
    Download the LCs for the target star for all the indicated sectors from the FFIs. This uses the lighkurve package.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
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
    all_md  :  list
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
    alltime12  :  list
        time used for the x and y centroid position plottin
    allfbkg  :  list
        background flux
    start_sec  :  list
        times of the start of the sector
    end_sec  :  list
        times of the end of the sector
    in_sec  :  list
        the sectors for which data was downloaded
    tessmag  :  float
        TESS magnitude of the target star
    teff  :  float
        effective temperature of the tagret star (K)
    srad  :  float
        radius of the target star (solar radii)

    '''

    # plot using the seaborn library
    sb.set(style='ticks')
    sb.set_color_codes('deep')

    future_sectors = list(set(sectors_all) - set(sector))

    # dowload the data

    searchtic = 'TIC' + str(tic)

    # empty lists to append the data that will be returned and used later
    start_sec = []
    end_sec = []
    in_sec = []

    alltime_list = []
    all_md = []

    X1_list = []
    X1flux_list = []
    X4_list = []
    apmask_list = []
    arrshape_list = []
    tpf_filt_list = []
    t_list = []
    bkg_list = []
    tpf_list = []

    # print a warning if more than one sector is chosen - as this will use more internet data and will take longer to run
    if len(sector) > 1:
        print ("Warning: Downloading data from multiple FFI sectors may take a while to run.")

    # open file which shows the momentum dumps
    momentum_dumps_list = "{}/data/tess_mom_dumps.txt".format(indir)
    mom_df = pd.read_csv(momentum_dumps_list, comment = '#', delimiter = '\t')

    for sec in sector:
        print (sec)
        # import the data
        print ("Importing FFI data sector {} ...".format(sec))

        # download data using lightkurve
        search_result = lk.search_tesscut(searchtic, sector=sec)

        # sometimes this failes so try multiple times... (faisl due to large data file)
        count = 0
        it_workd = False
        tpf = search_result.download(cutout_size=11)

        while not it_workd and count < 5:  # try up to 4 times to download it
            try:
                tpf = search_result.download(cutout_size=11)
                it_workd = True
            except:
                count += 1
                continue

        if count == 5: # if it fails 4 times, exit the program with an error message
            print ("ERROR: Could not download the tpf file at this time. Please try again.")
            sys.exit('')
        # ---------

        # get rid of the 'bad' quality data - use the data flags and only take data where the quality = 0.
        try:
            quality = tpf.quality
            tpf = tpf[quality == 0]
        except:
            print ('This file has no quality flag.')

        # save the raw tpf (target pixel file) data which will be used later - need to store it in this format as it will be used for the reprojection into the right coordinates.
        tpf_list.append(tpf)

        print ("done.\n")

        # extract the information and perform PCA
        print ("Start PCA analysis...", end =" ")

        # acess the flux information
        X1 = tpf.flux
        X1flux_list.append(X1)  # store tghe flux array in the list

        arrshape_list.append(X1.shape)  # store the shape in the list - not always the same shape but mostly 11 x 11 pixels.

        # calculate the average pixel value across all frames.
        # this will be used as the background color of the pixel level LCs plot.

        bkg = X1
        bkg = np.nanmean(bkg,axis = 0)

        bkg_list.append(bkg) # store the information in a list

        # reshape the array in order to perform PCA on it
        # reshape into 2 D (from 3D)
        s = X1.shape
        X1 = X1.reshape(s[0],s[1]*s[2])

        # maker sure that there are no NaN values
        lkeep = np.isfinite(X1.sum(axis=1)) * (X1.sum(axis=1)>0)
        X1 = X1[lkeep,:]

        # make 'similar' arrays filled with zeros
        X2 = np.zeros_like(X1)
        M = np.zeros(X1.shape[1])
        S = np.zeros(X1.shape[1])

        # normalize the data
        for n in range(X1.shape[1]):
            a = X1[:,n]
            x, m, s = norm(a)
            X2[:,n]=x
            M[n] = m
            S[n] = s

        ncomp = 2 # number of PCA components to consider
        pca = PCA(n_components=ncomp) # PCA
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

        all_md.append(list((mom_df.loc[mom_df['sec'] == int(sec)])['time']))

        start_sec.append([alltime[0]])
        end_sec.append([alltime[-1]])

        X1_list.append(X1) #  not corrected
        X4_list.append(X4) #  PCA corrected
        t_list.append(np.array(alltime))  # add it here because this list isn't flattened but the other one is



    alltime_list = list(np.hstack(alltime_list))
    all_md =       list(np.hstack(all_md) )

    del mom_df

    return alltime_list, all_md, start_sec, end_sec, in_sec, X1_list, X1flux_list,  X4_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list

def download_data_FFI(indir, sector, syspath, sectors_all, tic, save = False):
    '''
    Download the LCs for the target star for all the indicated sectors from the FFIs. This uses the lighkurve package.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
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
    all_md  :  list
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
    alltime12  :  list
        time used for the x and y centroid position plottin
    allfbkg  :  list
        background flux
    start_sec  :  list
        times of the start of the sector
    end_sec  :  list
        times of the end of the sector
    in_sec  :  list
        the sectors for which data was downloaded
    tessmag  :  float
        TESS magnitude of the target star
    teff  :  float
        effective temperature of the tagret star (K)
    srad  :  float
        radius of the target star (solar radii)

    '''
    # plot using the seaborn library
    sb.set(style='ticks')
    sb.set_color_codes('deep')


    future_sectors = list(set(sectors_all) - set(sector))

    # dowload the data

    searchtic = 'TIC' + str(tic)

    allfbkg = []
    allfbkg_t = []

    start_sec = []
    end_sec = []
    in_sec = []

    alltime_list = []
    allflux = []
    allflux_flat = []
    allflux_small = []
    all_md = []

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

    for idx, sec in enumerate(sector):

        # import the data
        print ("Importing FFI data sector {}...".format(sec), end =" ")

        searchtic = 'TIC' + str(tic)
        search_result = lk.search_tesscut(searchtic, sector=sec)

        count = 0
        it_workd = False

        while not it_workd and count < 5:  # try up to 4 times to download it
            try:
                tpf = search_result.download(cutout_size=15)
                it_workd = True
            except:
                count += 1
                continue

        if count == 5: # if it fails 4 times, exit the program with an error message
            print ("ERROR: Could not download the tpf file at this time. Please try again.")
            sys.exit('')
        # ---------

        # get rid of the 'bad' quality data - use the data flags and only take data where the quality = 0.
        try:
            quality = tpf.quality
            tpf = tpf[quality == 0]
        except:
            print ('This file has no quality flag.')

        tpf_list.append(tpf)

        print ("done.\n")

        # extract the information and perform PCA
        print ("Start PCA analysis...", end =" ")
        X1 = tpf.flux
        arrshape_list.append(X1.shape)

        # identify the target mask - this will need to be imporved - a bit 'hacky' at the moment.
        # the exctractiom mask size is based on the based on the average aperture size of the TESS pipeline apertures as afunction of radius.
        # this can be loosely quanified as:
        # magnitude = m*apsize + c   --> apsize = (magnitude - c) / m

        m = -0.415 # determined empirically
        c = 15.311

        # get the magnitude of this target
        magnitude = tpf.header['TESSMAG']
        #magnitude = 10

        # find the optimum aperture size
        large_ap_count = (magnitude - c) / m

        # aim to make the big aperture around 40 % smaller than the large aperture
        small_ap_count = round(large_ap_count * 0.60) # 60% of pipeline aperture

        # we are using the lighkurve optimization to extract the aperture.
        # this places an aperture on the central pixel and selects the brightest surrounding ones based on a threshhold value
        # determine this threshhold value based on the number of pixels that we want using scipy minimizaton and the 'find aperure' function as defined below under 'other functions'
        large_ap_thresh_val = minimize_scalar(find_aperture_FFI, bounds = [0,10], args = (large_ap_count, tpf)).x
        small_ap_thresh_val = minimize_scalar(find_aperture_FFI, bounds = [0,10], args = (small_ap_count, tpf)).x

        # using the optimal threshhold values, determine the masks
        target_mask = tpf.create_threshold_mask(threshold=large_ap_thresh_val, reference_pixel='center')
        target_mask_small = tpf.create_threshold_mask(threshold=small_ap_thresh_val, reference_pixel='center')


        # if this is run with an input file then the needed folder might not exist yet...
        if not os.path.exists("{}/{}/".format(indir, tic)): # if this folder doesn't already exist, make it
            os.makedirs("{}/{}/".format(indir, tic))

        # ---------

        # ----------
        # plot the mean image and plot the extraction apertures on top of it so that one can verify that the used apertures make sense
        im = np.nanmean(tpf.flux, axis = 0)
        # set up the plot - these are stored and one of the images saved in the report
        fig, ax = plt.subplots(1,2, figsize=(10,5), subplot_kw={'xticks': [], 'yticks': []})
        kwargs = {'interpolation': 'none', 'vmin': im.min(), 'vmax': im.max()}
        color = ['red', 'deepskyblue']
        label = ['small (~60 %)', 'pipeline (100 %)']

        # plot both the small and the
        for i, arr in enumerate([target_mask_small,target_mask]):

            mask = np.zeros(shape=(arr.shape[0], arr.shape[1]))
            mask = arr

            f = lambda x,y: mask[int(y),int(x)]
            g = np.vectorize(f)

            x = np.linspace(0,mask.shape[1], mask.shape[1]*100)
            y = np.linspace(0,mask.shape[0], mask.shape[0]*100)
            X, Y= np.meshgrid(x[:-1],y[:-1])
            Z = g(X[:-1],Y[:-1])

            ax[i].set_title('Aperture: {}'.format(label[i]), fontsize = 18)
            ax[i].imshow(im, cmap=plt.cm.viridis, **kwargs, origin = 'lower')
            ax[i].contour(Z, [0.5], colors=color[i], linewidths=[4],
                        extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])

        # save the figure
        plt.savefig('{}/{}/{}_apertures_{}.png'.format(indir, tic, tic, idx), format='png', bbox_inches = 'tight')
        plt.clf()
        plt.close()

        # ---------------

        bkg = X1
        bkg = np.nanmean(bkg,axis = 0)
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

        ncomp = 2
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

        print ("done.\n")

        # -------- flatten the normal lighcurve --------
        print ("Flatten LC...", end =" ")

        l = np.isfinite(flux)

        fr_inj = flux
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
        MAD = median_absolute_deviation(fr,ignore_nan = True)
        madrange = (5 * MAD * 100)

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

        print ("done.\n")

        print ("Extract background...", end =" ")
        # -------- extract the backrgound flux -----------
        background_mask = ~tpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
        n_background_pixels = background_mask.sum()

        background_lc_per_pixel = tpf.to_lightcurve(aperture_mask=background_mask) / n_background_pixels
        n_target_pixels = target_mask.sum()
        background_estimate_lc = background_lc_per_pixel * n_target_pixels

        allfbkg.append(background_estimate_lc.flux) # uncorrected background
        allfbkg_t.append(background_estimate_lc.time)

        print ("done.\n")
         # ------------------------------------------------

        # add all the information to lists that are returned at the end.
        tpf_filt_list.append(X4.reshape(tpf.flux[lkeep,:,:].shape))

        in_sec.append(sec)

        alltime_list.append(list(alltime))
        allflux_flat.append(list(fr))
        allflux.append(list(flux))
        allflux_small.append(list(flux_small))

        all_md.append(list((mom_df.loc[mom_df['sec'] == int(sec)])['time']))

        start_sec.append([alltime[0]])
        end_sec.append([alltime[-1]])

        X1_list.append(X1) #  not corrected
        X4_list.append(X4) #  PCA corrected
        t_list.append(np.array(alltime))  # add it here because this list isn't flattened but the other one is

    alltime_list =  list(np.hstack(alltime_list) )
    allflux_flat =  list(np.hstack(allflux_flat) )
    allflux_small = list(np.hstack(allflux_small))
    allflux =       list(np.hstack(allflux)      )
    all_md =        list(np.hstack(all_md)       )
    allfbkg =       list(np.hstack(allfbkg)      )
    allfbkg_t =     list(np.hstack(allfbkg_t)    )

    del mom_df

    return alltime_list, allflux, allflux_small, allflux_flat, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list


def download_data_neighbours(indir, sector, tics, distance, binfac = 5):
    ''''
    Dowloand the data for the 6 nearest neightbour target pixel LC - in the future want to
    take this so it used the FFI data as this will be closer targets and a better diagnostic.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    sector  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    distance:  list
        list of the distance to the nearest neighbouring tics from the target (in arcminutes)
    tics : str
        TIC (Tess Input Catalog) ID of the target
    binfac  :  int
        The factor by which the data should be binned. Default = 5 (which is what is shown on PHT)

    Returns
    ------
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    all_md  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    transit_list  :  list
        list of all the marked transit-events
    outtics  :  list
        the tic ids of the nearby targets
    tessmag_list  :  list
        the magnitudes of the nearby targets
    distance  :  list
        the distances to the nearby targets
    '''

    sb.set(style='ticks')
    sb.set_color_codes('deep')

    #try:
    #    lc_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_lc.sh'.format(indir, str(sector)), dtype = str)
    #except:
    #    print ("Sector {}  has not yet been osberved - come back to this target later.".format(sector))

    alltime_neighbours = []
    allflux_neighbours = []
    all_md_neighbours = []
    alltimebinned_neighbours = []
    allfluxbinned_neighbours = []
    outtics = []

    dwload_link = []

    infile = pd.read_csv("{}/data/sector_download_codes.txt".format(indir), delimiter = ' ', names = ['sec', 'first', 'second'], comment = '#')

    for tic in tics:

        this_sector_code = infile[infile.sec == int(sector)]

        sector = str(sector)

        download_url = (
            "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
            + str(this_sector_code['first'].values[0]).rjust(13, "0")
            + "-s"
            + sector.rjust(4, "0")
            + "-"
            + str(tic).rjust(16, "0")
            + "-"
            + str(this_sector_code['second'].values[0]).rjust(4, "0")
            + "-s_lc.fits")

        dwload_link.append(download_url)


    alltimebinned = []
    allfluxbinned = []

    start_sec = []
    end_sec = []

    alltime = []
    allflux = []
    all_md = []

    tessmag_list = []

    for num,lcfile in enumerate(dwload_link):

        print ("Downloading nearest neighbour   {}   of   {}....".format(num + 1, len(dwload_link), end ='' ))
        try:
            response = requests.get(lcfile)

            # open the file using the response url
            lchdu  = pf.open(response.url) # this needs to be a URL - not a file
            outtics.append(int(lchdu[0].header['TICID']))

            #Open and view columns in lightcurve extension
            lcdata = lchdu[1].data
            lchdu[1].columns

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

            lchdu.close()

            # binned data
            N       = len(time)
            n       = int(np.floor(N/binfac)*binfac)
            X       = np.zeros((2,n))
            X[0,:]  = time[:n]
            X[1,:]  = f1[:n]
            Xb      = rebin(X, (2,int(n/binfac)))

            time_binned    = Xb[0]
            flux_binned    = Xb[1]

            mom_dump = np.bitwise_and(quality, 2**5) >= 1

            alltime.append(list(time)) #[~bad_data]
            allflux.append(list(f1)) #[~bad_data])
            all_md.append(list(time[mom_dump]))

            alltimebinned.append(list(time_binned))
            allfluxbinned.append(list(flux_binned))

            start_sec.append([time[0]])
            end_sec.append([time[-1]])
            tessmag_list.append(tessmag)

            print("done.")
        except:
            continue

    return alltime, allflux, all_md, alltimebinned, allfluxbinned, outtics, tessmag_list, distance


def check_FFI_data_products(tic,sector, args):
    '''
    Functino to check what data products are available for the FFIS.
    '''

    # query MAST to find out what data products are available
    objectOfInterest = 'TIC {}'.format(tic)

    resolverRequest = {'service':'Mast.Name.Lookup',
                        'params':{'input':objectOfInterest,
                                    'format':'json'}}

    headers,resolvedObjectString = mastQuery(resolverRequest)

    resolvedObject = json.loads(resolvedObjectString)

    objRa = resolvedObject['resolvedCoordinate'][0]['ra']
    objDec = resolvedObject['resolvedCoordinate'][0]['decl']


    mastRequest = {'service':'Mast.Caom.Cone',
                    'params':{'ra':objRa,
                            'dec':objDec,
                            'radius':0.01},
                    'format':'json',
                    'pagesize':2000,
                    'page':1,
                    'removenullcolumns':True,
                    'cachebreaker':"2020-02-12T0749"}

    headers,mastDataString = mastQuery(mastRequest)

    mastData = json.loads(mastDataString)

    # convert to a table that makes sense
    mastDataTable = Table()


    for col,atype in [(x['name'],x['type']) for x in mastData['fields']]:
        if atype=="string":
            atype="str"
        if atype=="boolean":
            atype="bool"
        mastDataTable[col] = np.array([x.get(col,None) for x in mastData['data']],dtype=atype)

    mastDataTable = mastDataTable.to_pandas()
    # filter for the sectors that we want to look at
    mastDataTable = mastDataTable[mastDataTable['sequence_number'].isin(sector)]


    # filter for TESS-SPOC data products
    TESSSPOC_mastDataTable  = mastDataTable[mastDataTable['provenance_name'] == 'TESS-SPOC']
    QLP_mastDataTable  = mastDataTable[mastDataTable['provenance_name'] == 'QLP']

    base_url = 'https://mast.stsci.edu/api/v0.1/Download/file/?uri='

    if len(TESSSPOC_mastDataTable) > 0:
        # this is the best option as the data are reduced using the PSOC pipeline so we have all of the same data products as for the 2 minute cadence.
        args.FFI = 'SPOC'

        # these are the data urls that we need to download the data reduced using the spoc tess pipeline (similar to 2 min candence data)
        data_urls = list(TESSSPOC_mastDataTable['dataURL'])

    elif len(QLP_mastDataTable) > 0:
        # this is the best option as the data are reduced using the PSOC pipeline so we have all of the same data products as for the 2 minute cadence.
        args.FFI = 'QLP'

        if len(mastDataTable) > len(sector):
            QLP_mastDataTable_df = pd.DataFrame(QLP_mastDataTable)
            QLP_mastDataTable = QLP_mastDataTable_df.drop_duplicates(subset ='sequence_number')
            
        # these are the data urls that we need to download the data reduced using the spoc tess pipeline (similar to 2 min candence data)
        data_urls = QLP_mastDataTable['dataURL']

    else:
        args.FFI = True
        data_urls = [-111]

    return data_urls, args


def download_data_FFI_SPOC(indir, sector, syspath, sectors_all, tic, args, data_urls, save = False):

    allx1 = []
    allx2 = []
    ally1 = []
    ally2 = []
    alltime12 = []
    allfbkg = []

    start_sec = []
    end_sec = []
    in_sec = []

    alltime = []
    allflux = []
    allflux_err = []
    all_md = []

    failed_sectors = []
    # loop through all the download links - all the data that we want to access
    base_url = 'https://mast.stsci.edu/api/v0.1/Download/file/?uri='

    for lcfile in data_urls:

        # !-!-!-!-!-!-!-
        # if this a test run, download the file already on the system
        #if test != 'no':
        #    lchdu  = pf.open(lcfile)
        ## !-!-!-!-!-!-!-

        try:
            # use the downlload link to download the file from the server - need an internet connection for this to work
            response = requests.get("{}{}".format(base_url, lcfile))
            # open the file using the response url
            lchdu  = pf.open(response.url) # this needs to be a URL - not a file
        except:
            continue

        # open and view columns in lightcurve extension
        lcdata = lchdu[1].data

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


        allx1.append(list(x1[l2]))
        allx2.append(list(x2[l2]))
        ally1.append(list(y1[l2]))
        ally2.append(list(y2[l2]))

        lchdu.close()

        # store the sector we are looking at
        in_sec.append(sec)

        # the time of the momentum dumps are indicated by the quality flag
        mom_dump = np.bitwise_and(quality, 2**5) >= 1

        # store the relevant information in the given list
        alltime.append(list(time))
        allflux.append(list(f1))
        allflux_err.append(list(f1_err))
        all_md.append(list(time[mom_dump]))

        alltime12.append(list(time[l2]))

        allfbkg.append(fbkg)

        start_sec.append([time[0]])
        end_sec.append([time[-1]])

    # flatten lists of list
    alltime =       list(np.hstack(alltime)      )
    allflux =       list(np.hstack(allflux)      )
    allflux_err =   list(np.hstack(allflux_err)  )
    all_md =        list(np.hstack(all_md)       )

    allx1 =         list(np.hstack(allx1)        )
    allx2 =         list(np.hstack(allx2)        )
    ally1 =         list(np.hstack(ally1)        )
    ally2 =         list(np.hstack(ally2)        )
    alltime12 =     list(np.hstack(alltime12)    )
    allfbkg =       list(np.hstack(allfbkg)      )

    if len (failed_sectors) != 0:
        print ("Failed sectors: {}    (These may be future sectors of may have only been observed in the FFIs".format(failed_sectors))

    return alltime, allflux, allflux_err, all_md, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, data_urls, args



def download_data_FFI_QLP(indir, sector, syspath, sectors_all, tic, args, data_urls, save = False):


    alltime12 = []
    allfbkg = []

    start_sec = []
    end_sec = []
    in_sec = []

    alltime = []
    allflux = []
    allflux_err = []
    allflux_small = []
    allflux_large = []
    all_md = []

    X1_list = []
    X4_list = []
    arrshape_list = []
    tpf_filt_list = []
    t_list = []
    tpf_list = []
    bkg_list = []

    failed_sectors = []
    # loop through all the download links - all the data that we want to access

    base_url = 'https://mast.stsci.edu/api/v0.1/Download/file/?uri='
    searchtic = 'TIC' + str(tic)

    for lcfile in data_urls:

        # !-!-!-!-!-!-!-
        # if this a test run, download the file already on the system
        #if test != 'no':
        #    lchdu  = pf.open(lcfile)
        ## !-!-!-!-!-!-!-

        #try:
            # use the downlload link to download the file from the server - need an internet connection for this to work
        response = requests.get("{}{}".format(base_url, lcfile))
        # open the file using the response url
        lchdu  = pf.open(response.url) # this needs to be a URL - not a file
        #except:
        #    continue

        # open and view columns in lightcurve extension
        lcdata = lchdu[1].data

        f02 = lcdata['KSPSAP_FLUX'] # Presearch Data Conditioning
        f02_small = lcdata['KSPSAP_FLUX_SML'] # Presearch Data Conditioning
        f02_large = lcdata['KSPSAP_FLUX_LAG'] # Presearch Data Conditioning

        f02_err = lcdata['KSPSAP_FLUX_ERR']

        quality = lcdata['QUALITY']  # quality flags as determined by the SPOC pipeline
        time    = lcdata['TIME']

        f0     = lcdata['SAP_FLUX']  #
        fbkg     = lcdata['SAP_BKG']  # background flux

        f1 = f02/np.nanmedian(f02)  # normalize by dividing byt the median flux
        f1_small = f02_small/np.nanmedian(f02_small)  # normalize by dividing byt the median flux
        f1_large = f02_large/np.nanmedian(f02_large)  # normalize by dividing byt the median flux
        f1_err = f02_err/np.nanmedian(f02)  # normalise the errors on the flux


        l       = (quality>0)   # good quality data
        l2      = (quality<=16)  # bad quality data
        sec     = int(lchdu[0].header['SECTOR'])  # the TESS observational sector
        tessmag = lchdu[0].header['TESSMAG']  # magnitude in the FITS header
        teff    = lchdu[0].header['TEFF']     # effective temperature in the FITS header (kelvin)
        srad    = lchdu[0].header['RADIUS']   # stellar radius in the FITS header (solar radii)

        lchdu.close()

        # store the sector we are looking at
        in_sec.append(sec)

        # the time of the momentum dumps are indicated by the quality flag
        mom_dump = np.bitwise_and(quality, 2**5) >= 1

        # store the relevant information in the given list
        alltime.append(list(time))
        allflux.append(list(f1))
        allflux_small.append(list(f1_small))
        allflux_large.append(list(f1_large))

        allflux_err.append(list(f1_err))
        all_md.append(list(time[mom_dump]))

        alltime12.append(list(time[l2]))

        allfbkg.append(fbkg)

        start_sec.append([time[0]])
        end_sec.append([time[-1]])

    # flatten lists of list
    alltime =              list(np.hstack(alltime)      )
    allflux =              list(np.hstack(allflux)      )
    allflux_small =        list(np.hstack(allflux_small)      )
    allflux_large =        list(np.hstack(allflux_large)      )
    allflux_err =          list(np.hstack(allflux_err)  )
    all_md =               list(np.hstack(all_md)       )
    allfbkg =              list(np.hstack(allfbkg)      )

    # get the TP data (using lighktkurve)
    for idx, sec in enumerate(sector):

        # import the data
        print ("Importing FFI data sector {}...".format(sec), end =" ")
        search_result = lk.search_tesscut(searchtic, sector=sec)

        count = 0
        it_workd = False

        while not it_workd and count < 5:  # try up to 4 times to download it
            try:
                tpf = search_result.download(cutout_size=15)
                it_workd = True
            except:
                count += 1
                continue

        if count == 5: # if it fails 4 times, exit the program with an error message
            print ("ERROR: Could not download the tpf file at this time. Please try again.")
            sys.exit('')

        # get rid of the 'bad' quality data - use the data flags and only take data where the quality = 0.
        try:
            quality = tpf.quality
            tpf = tpf[quality == 0]
        except:
            print ('This file has no quality flag.')



        X1 = tpf.flux
        arrshape_list.append(X1.shape)

        if not os.path.exists("{}/{}/".format(indir, tic)): # if this folder doesn't already exist, make it
            os.makedirs("{}/{}/".format(indir, tic))

        bkg = X1
        bkg = np.nanmean(bkg,axis = 0)
        bkg_list.append(bkg) # store the information in a list

        # reshape the array in order to perform PCA on it.
        s = X1.shape
        X1 = X1.reshape(s[0],s[1]*s[2])

        lkeep = np.isfinite(X1.sum(axis=1)) * (X1.sum(axis=1)>0)
        X1 = X1[lkeep,:]
        t=tpf.time[lkeep].value

        #X2 = np.zeros_like(X1)
        #M = np.zeros(X1.shape[1])
        #S = np.zeros(X1.shape[1])

        #for n in range(X1.shape[1]):
        #    a = X1[:,n]
        #    x, m, s = norm(a)
        #    X2[:,n]=x
        #    M[n] = m
        #    S[n] = s

        #ncomp = 2
        #pca = PCA(n_components=ncomp)
        #trends = pca.fit_transform(X2)
        #weights = pca.components_

        #X3 = np.copy(X2)
        #for n in range(X2.shape[1]):
        #    for m in range(ncomp):
        #        X3[:,n] -= trends[:,m] * weights[m,n]

        #X4 = np.zeros_like(X3)
        #for n in range(X2.shape[1]):
        #    x = X3[:,n]
        #    a = unnorm(x, M[n], S[n])
        #    X4[:,n] = a


        X1_list.append(X1) #  not corrected
        X4_list.append(X1) #  PCA corrected
        t_list.append(np.array(t))
        tpf_list.append(tpf)
        tpf_filt_list.append(X1.reshape(tpf.flux[lkeep,:,:].shape))


    if len (failed_sectors) != 0:
        print ("Failed sectors: {}    (These may be future sectors of may have only been observed in the FFIs".format(failed_sectors))

    return alltime, allflux, allflux_small, allflux_err, all_md, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, X1_list, X4_list, bkg_list, t_list, tpf_list, tpf_filt_list, arrshape_list



def download_tpf(indir, transit_sec, transit_list, tic, url_list, test = 'no'):

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
        list of the filtered tpfs

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
    tpf_list   :  list
        list of the tpf files - these are needed in order to access the projection for the rotation of images later.

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

    TESS_unbinned_t_l = []
    TESS_binned_t_l = []
    small_binned_t_l= []
    TESS_unbinned_l= []
    TESS_binned_l= []
    small_binned_l= []
    tpf_list = []

    # -!-!-!-!-!-!-!-
    # if this is a test run, download the file that's already on the system as to not rely on the server response
    if test != 'no':
        dwload_link_tp = [test]  # in the unittest define a link to a pre downloaded file
    # -!-!-!-!-!-!-!-


    elif url_list != [-111]:
        dwload_link_tp = []

        for url in url_list:
            for sector in transit_sec:
                if "{}/target".format(sector) in (url):
                    dwload_link_tp.append('https://mast.stsci.edu/api/v0.1/Download/file/?uri=' + url.replace('lc.fits', 'tp.fits'))

    # if not a test find the link needed to access the data
    else:
        dwload_link_tp = []

        sectorfile = "{}/data/sector_download_codes.txt".format(indir)
        infile = pd.read_csv(sectorfile, delimiter = ' ', names = ['sec', 'first', 'second'], comment = '#')

        # search the LC download file for the URL to download the LC from MAST for the given target star.
        # if we are looking at all of the sectors search, the file that has all of the LC URLs in it.


        for sector in transit_sec:
            sector = str(sector)
            try:
                this_sector_code = infile[infile.sec == int(sector)]

                sector = str(sector)

                download_url = (
                    "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
                    + str(this_sector_code['first'].values[0]).rjust(13, "0")
                    + "-s"
                    + sector.rjust(4, "0")
                    + "-"
                    + str(tic).rjust(16, "0")
                    + "-"
                    + str(this_sector_code['second'].values[0]).rjust(4, "0")
                    + "-s_tp.fits")
                dwload_link_tp.append(download_url)
            except:
                print ("file can't be downlaoded at this moment")

    # download each file (i.e. each sector)

    for idx, file in enumerate(dwload_link_tp):
        try:
            response = requests.get(file)
            tpf = pf.open(response.url)   # open the file

        except:
            print ("\n !!! This target pixel file is corrupt and cannot be downloaded at this time. Please try again later or a different file. \n")
            return [-111], [-111], [-111], [-111], [-111], [-111], [-111],[-111], [-111], [-111], [-111], [-111], [-111], [-111], [-111], [-111], [-111] # flag an error message


        for T0 in transit_list:

            tpf_list.append(tpf)

            X1_original = tpf[1].data['FLUX']
            arrshape = X1_original.shape

            Q = tpf[1].data['QUALITY']
            t = tpf[1].data['TIME']

            aperture = tpf[2].data
            apmask = aperture >= np.nanmax(aperture)


            # check whether any of the transits are in that sector - if not skip that sector.
            if (T0 > np.nanmin(t)) and (T0 < np.nanmax(t)):

                lkeep = (Q==0)
                bkg = X1_original[lkeep,:]
                bkg = np.nanmean(bkg,axis = 0)

                bkg_list.append(bkg)
                apmask_list.append(apmask)
                arrshape_list.append(arrshape)

                s = X1_original.shape
                X1 = X1_original.reshape(s[0],s[1]*s[2])

                lkeep = np.isfinite(X1.sum(axis=1)) * (Q==0) * (X1.sum(axis=1)>0)
                X1 = X1[lkeep,:]
                t=t[lkeep]

                # - - - - - - - - - - - - -
                # exctract LCs using different aperture sizes....

                # tess mask
                pipeline_mask  = tpf[2].data & 2 > 0  # the tess mask has values greater than two
                tess_ap_size = np.sum(pipeline_mask)  # size of the tess mask - needed to later find a smaller mask

                binfac = 7

                # extract the images given the apertures
                # pipeline apertures
                flux = X1[:,pipeline_mask.flatten()].sum(axis=1)
                flux /= np.nanmedian(flux)

                # bin the data
                N       = len(t)
                n       = int(np.floor(N/binfac)*binfac)
                X       = np.zeros((2,n))
                X[0,:]  = t[:n]
                X[1,:]  = flux[:n]
                Xb      = rebin(X, (2,int(n/binfac)))

                time_binned    = Xb[0]
                flux_binned    = Xb[1]

                # aim to make the big aperture around 40 % smaller than the TESS pipeline aperture
                small_ap_count = round(tess_ap_size * 0.60) # 60% of pipeline aperture

                # we are using a lighkurve function optimization to extract the aperture.
                # this places an aperture on the central pixel and selects the brightest surrounding ones based on a threshhold value
                # determine this threshhold value based on the number of pixels that we want using scipy minimizaton and the 'find aperure' function as defined below under 'other functions'
                small_ap_thresh_val = minimize_scalar(find_aperture, bounds = [0,10], args = (small_ap_count,X1_original)).x

                # using the optimal threshhold values, determine the masks
                smaller_mask = extract_aperture(X1_original, small_ap_thresh_val)

                flux_small = X1[:,smaller_mask.flatten()].sum(axis=1)
                flux_small /= np.nanmedian(flux_small)

                # bin the flux from the smaller aperture
                N       = len(t)
                n       = int(np.floor(N/binfac)*binfac)
                X       = np.zeros((2,n))
                X[0,:]  = t[:n]
                X[1,:]  = flux_small[:n]
                Xb      = rebin(X, (2,int(n/binfac)))

                time_binned_small    = Xb[0]
                flux_binned_small    = Xb[1]

                # ----------
                # time
                TESS_unbinned_t_l.append(t)
                TESS_binned_t_l.append(time_binned)
                small_binned_t_l.append(time_binned_small)

                # flux
                TESS_unbinned_l.append(flux)
                TESS_binned_l.append(flux_binned)
                small_binned_l.append(flux_binned_small)


                # - - - - - - - - - - - - - - - - - - - - -
                # do PCA analysis to try to detrend the LC.

                X2 = np.zeros_like(X1)
                M = np.zeros(X1.shape[1])
                S = np.zeros(X1.shape[1])
                for n in range(X1.shape[1]):
                    a = X1[:,n]
                    x, m, s = norm(a)
                    X2[:,n]=x
                    M[n] = m
                    S[n] = s

                ncomp = 2
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


                tpf_filt = X4.reshape(tpf[1].data['FLUX'][lkeep,:,:].shape)

                oot = (abs(T0-t) < 0.3) * (abs(T0-t) < 0.3)
                intr = abs(T0-t) < 0.1

                X1_list.append(X1) # not corrected
                X4_list.append(X4) #  PCA corrected
                oot_list.append(oot)  # out of transit filter
                intr_list.append(intr)  # in transit filter
                t_list.append(t)
                T0_list.append(T0)
                tpf_filt_list.append(tpf_filt)



                # ----------
                # plot the mean image and plot the extraction apertures on top of it so that one can verify that the used apertures make sense
                im = np.nanmean(X1_original, axis = 0)
                # set up the plot - these are stored and one of the images saved in the report
                fig, ax = plt.subplots(1,2, figsize=(10,5), subplot_kw={'xticks': [], 'yticks': []})
                kwargs = {'interpolation': 'none', 'vmin': im.min(), 'vmax': im.max()}
                color = ['red', 'deepskyblue']
                label = ['small (~60 %)', 'pipeline (100 %)']


        for i,arr in enumerate([smaller_mask,pipeline_mask]):

            mask = np.zeros(shape=(arr.shape[0], arr.shape[1]))
            mask= arr

            f = lambda x,y: mask[int(y),int(x)]
            g = np.vectorize(f)

            x = np.linspace(0,mask.shape[1], mask.shape[1]*100)
            y = np.linspace(0,mask.shape[0], mask.shape[0]*100)
            X, Y= np.meshgrid(x[:-1],y[:-1])
            Z = g(X[:-1],Y[:-1])

             # get the pixels for the columns and rows
            start_column = tpf[1].header['1CRV5P']
            start_row = tpf[1].header['2CRV5P']

            # only show the start and end pixel
            y_label_list = ['{}'.format(start_row),'{}'.format(start_row + im.shape[0])] # make a pixel label for the x and y axis
            ax[i].set_yticks([0, y[:-1].max()-1])
            ax[i].set_yticklabels(y_label_list)

            x_label_list = ['{}'.format(start_column),'{}'.format(start_column + im.shape[1])]
            ax[i].set_xticks([0, x[:-1].max()-1])
            ax[i].set_xticklabels(x_label_list)

            ax[i].set_title('Aperture: {}'.format(label[i]), fontsize = 18)
            ax[i].imshow(im, cmap=plt.cm.viridis, **kwargs, origin = 'lower')
            ax[i].contour(Z, [0.5], colors=color[i], linewidths=[4],
                        extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])

            ax[i].set_xlabel("column (pixels)", labelpad=-9.5)
            ax[i].set_ylabel("row (pixels)", labelpad=-20)

        # save the figure (unless this is in TEST mode)
        if test == 'no':
            plt.savefig('{}/{}/{}_apertures_{}.png'.format(indir, tic, tic, idx), format='png', bbox_inches = 'tight')

        plt.clf()
        plt.close()



    TESS_unbinned_t_l = list(np.hstack(TESS_unbinned_t_l))
    TESS_binned_t_l   = list(np.hstack(TESS_binned_t_l))
    small_binned_t_l =  list(np.hstack(small_binned_t_l))

    TESS_unbinned_l =   list(np.hstack(TESS_unbinned_l))
    TESS_binned_l =     list(np.hstack(TESS_binned_l))
    small_binned_l =    list(np.hstack(small_binned_l))


    return X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list


def download_tpf_lightkurve(indir, transit_list, sector, tic, test = 'no'):
    '''

    NO LONGER USED - WILL BE DELETED

    function to download the tpf data using LightKurve. This is used in order to extract the LC in different aperture sizes.

    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_list  : list
        list of the marked transit-events
    sector  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    tic : str
        TIC of the target

    Returns
    ------
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
    tpf_list   :  list
        list of the tpd files - these are needed in order to access the projection for the rotation of images later.

    '''

    # -!-!-!-!-!-!-!-
    # if this is a test run, download the file that's already on the system as to not rely on the server response
    if test != 'no':
        dwload_link_tp = [test]  # in the unittest define a link to a pre downloaded file
    # -!-!-!-!-!-!-!-

    # if not test, find the link needed to downlaod the data

    else:
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
            print ("\n (Also check that the data files have been downloaded - run code with '--new-data' in the command line)")

        dwload_link_tp = dwload_link


    TESS_unbinned_t_l = []
    TESS_binned_t_l = []
    small_binned_t_l= []
    TESS_unbinned_l= []
    TESS_binned_l= []
    small_binned_l= []
    tpf_list = []

    for idx,file in enumerate(dwload_link_tp):

        try:
            tpf = lk.TessTargetPixelFile(file) # dowload the Target Pixel File
        except:
            # on very occasioal files the target pixel file is corrupt and cannot download with the TypeError: Buffer is too small.
            print ("\n !!! This target pixel file is corrupt and cannot be downloaded at this time. Please try again later or a different file. \n")
            return -111, -111, -111, -111, -111, -111, -111 # flag an error message

        # if the transit is within that sector then add it to the list to be used later
        for T0 in transit_list:
            if (T0 > np.nanmin(tpf.time)) and (T0 < np.nanmax(tpf.time)):
                tpf_list.append(tpf)

        try:
            # calculate the masks that we want to look at - three masks (small, pipeline mask and large)
            # determine the size of the tess pipeline mask - this correlates with the magnitude of the target star
            tess_ap_size = np.sum(tpf.pipeline_mask)

            # aim to make the big aperture around 40 % smaller than the TESS pipeline aperture
            small_ap_count = round(tess_ap_size * 0.60) # 60% of pipeline aperture

            # we are using the lighkurve optimization to extract the aperture.
            # this places an aperture on the central pixel and selects the brightest surrounding ones based on a threshhold value
            # determine this threshhold value based on the number of pixels that we want using scipy minimizaton and the 'find aperure' function as defined below under 'other functions'
            small_ap_thresh_val = minimize_scalar(find_aperture, bounds = [0,10], args = (small_ap_count, tpf)).x

            # using the optimal threshhold values, determine the masks
            smaller_mask = tpf.create_threshold_mask(threshold=small_ap_thresh_val, reference_pixel='center')

            # TESS binned
            TESS_unbinned = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask).flatten(window_length=100001)
            TESS_unbinned = TESS_unbinned.remove_outliers(6)

            # TESS binned
            TESS_binned = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask).flatten(window_length=100001)
            TESS_binned = TESS_binned.remove_outliers(6).bin(7)

            # Use a custom aperture binned
            small_binned = tpf.to_lightcurve(aperture_mask=smaller_mask).flatten(window_length=100001)
            small_binned = small_binned.remove_outliers(6).bin(7)

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


            # ----------
            # plot the mean image and plot the extraction apertures on top of it so that one can verify that the used apertures make sense
            im = np.nanmean(tpf.flux, axis = 0)
            # set up the plot - these are stored and one of the images saved in the report
            fig, ax = plt.subplots(1,2, figsize=(10,5), subplot_kw={'xticks': [], 'yticks': []})
            kwargs = {'interpolation': 'none', 'vmin': im.min(), 'vmax': im.max()}
            color = ['red', 'deepskyblue']
            label = ['small (~60 %)', 'pipeline (100 %)']


            for i,arr in enumerate([smaller_mask,tpf.pipeline_mask]):

                mask = np.zeros(shape=(arr.shape[0], arr.shape[1]))
                mask= arr

                f = lambda x,y: mask[int(y),int(x)]
                g = np.vectorize(f)

                x = np.linspace(0,mask.shape[1], mask.shape[1]*100)
                y = np.linspace(0,mask.shape[0], mask.shape[0]*100)
                X, Y= np.meshgrid(x[:-1],y[:-1])
                Z = g(X[:-1],Y[:-1])

                # get the pixels for the columns and rows
                start_column = tpf[1].header['1CRV5P']
                start_row = tpf[1].header['2CRV5P']

                # only show the start and end pixel
                y_label_list = ['{}'.format(start_row),'{}'.format(start_row + im.shape[0])] # make a pixel label for the x and y axis
                ax[i].set_yticks([0, y[:-1].max()-1])
                ax[i].set_yticklabels(y_label_list)

                x_label_list = ['{}'.format(start_column),'{}'.format(start_column + im.shape[1])]
                ax[i].set_xticks([0, x[:-1].max()-1])
                ax[i].set_xticklabels(x_label_list)

                ax[i].set_title('Aperture: {}'.format(label[i]), fontsize = 18)
                ax[i].imshow(im, cmap=plt.cm.viridis, **kwargs, origin = 'lower')
                ax[i].contour(Z, [0.5], colors=color[i], linewidths=[4],
                            extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])

                ax[i].set_xlabel("column (pixels)", labelpad=-9.5)
                ax[i].set_ylabel("row (pixels)", labelpad=-20)

            # save the figure
            plt.savefig('{}/{}/{}_apertures_{}.png'.format(indir, tic, tic, idx), format='png', bbox_inches = 'tight')
            plt.clf()
            plt.close()

        except:
            continue


    TESS_unbinned_t_l = list(np.hstack(TESS_unbinned_t_l))
    TESS_binned_t_l   = list(np.hstack(TESS_binned_t_l)  )
    small_binned_t_l =  list(np.hstack(small_binned_t_l) )

    TESS_unbinned_l =   list(np.hstack(TESS_unbinned_l)  )
    TESS_binned_l =     list(np.hstack(TESS_binned_l)    )
    small_binned_l =    list(np.hstack(small_binned_l)   )


    return TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list

def download_tpf_mast(indir, transit_sec, transit_list, tic, test = 'no'):
    '''

    NO LONGER USED - WILL BE DELETED

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


    # -!-!-!-!-!-!-!-
    # if this is a test run, download the file that's already on the system as to not rely on the server response
    if test != 'no':
        dwload_link_tp = [test]  # in the unittest define a link to a pre downloaded file
    # -!-!-!-!-!-!-!-


    # if not test, find the link needed to downlaod the data
    else:

        dwload_link_tp = []

        for sec in transit_sec: #the sector that this image is in

            tpf_all = np.genfromtxt('{}/data/tesscurl_sector_{}_tp.sh'.format(indir,sec), dtype = str)

            for i in tpf_all:

                if str(tic) in str(i[6]):
                    dwload_link_tp.append(i[6])

    # download each file (i.e. each sector)
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
                bkg = np.nanmean(bkg,axis = 0)

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

                ncomp = 2
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


                oot = (abs(T0-t) < 1.5) * (abs(T0-t) < 1.5)
                intr = abs(T0-t) < 0.1

                X1_list.append(X1) # not corrected
                X4_list.append(X4) #  PCA corrected
                oot_list.append(oot)  # out of transit filter
                intr_list.append(intr)  # in transit filter
                t_list.append(t)
                T0_list.append(T0)
                tpf_filt_list.append(tpf_filt)


    return X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list

def data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, args):
    '''
    function that runs the BLS routine and plots the results. The BLS is run twice and in the second
    the most significant result found in the first run is removed.
    Prior to running the BLS the data is detrended.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times for the LC of the target
    allflux  :  list
        normalized flux for the LC of the target (not binned)
    alltimebinned  :  list
        binned times for the target
    allfluxbinned  :  list
        normalized binned flux for the LC of the target

    Returns
    -------
        two lists of the statistics of the to BLS runs. Each list contains:
    stats_period
    stats_t0
    stats_duration
    stats_depth
    stats_depth_phased
    stats_depth_half
    stats_depth_odd
    stats_depth_even

    '''
    # make sure that there are no nan value sin the data - they cause everything to crash
    mask_binned = np.isfinite(alltimebinned) * np.isfinite(allfluxbinned)
    mask = np.isfinite(alltime) * np.isfinite(allflux)

    alltimebinned = np.array(alltimebinned)[mask_binned]
    allfluxbinned = np.array(allfluxbinned)[mask_binned]
    alltime = np.array(alltime)[mask]
    allflux = np.array(allflux)[mask]

    # -------------------

    # detrend the data before running the BLS
    print ("Detrend data prior to BLS analysis...", end =" ")

    fr_inj = allfluxbinned

    T_dur = 0.7  #The transit duration - may need to change!! Larger means it corrects less which is probably better.

    nmed = int(144*3*T_dur) # 144 because the data is binned to 10 minutes and 24 hours / 10 mins = 144 data points
    nmed = 2*int(nmed/2)+1 # make it an odd number
    ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
    # first number (nmed) is three time transit durations, the second quite small (10,20 )

    l = np.isfinite(ff)
    g = interp1d(alltimebinned[l],ff[l],bounds_error=False,fill_value=np.nan)
    ff = g(alltimebinned)

    allfluxbinned = fr_inj / ff

    print ("done.")

    # -------------------
    if args.save == True:
        fix, ax = plt.subplots(2,1,figsize=(16,12))

        ax[0].plot(alltimebinned,fr_inj, '.',label = 'Uncorrected')
        ax[0].plot(alltimebinned,ff,'.',label = 'Model fit')
        ax[0].legend(fontsize = 16, loc = 1)

        ax[1].plot(alltimebinned,allfluxbinned, '.', color = 'navy', label = 'Clipped + corrected')
        ax[1].legend(fontsize = 16, loc = 1)

        ax[1].set_xlabel("Time", fontsize = 16)
        ax[0].set_ylabel("Flux", fontsize = 16)
        ax[1].set_ylabel("Flux", fontsize = 16)
        ax[1].set_ylabel("Flux", fontsize = 16)

        plt.savefig('{}/{}/{}_fit_test.png'.format(indir, tic, tic), format='png', bbox_inches = 'tight')
        plt.clf()
        plt.close()


    # get rid of nans again...
    mask_binned = np.isfinite(alltimebinned) * np.isfinite(allfluxbinned)
    alltimebinned = np.array(alltimebinned)[mask_binned]
    allfluxbinned = np.array(allfluxbinned)[mask_binned]
    # -----------------------

    durations = np.linspace(0.05, 0.5, 15) # ????? CHECK THESE
    periods = np.arange(0.51, (np.nanmax(alltimebinned) - np.nanmin(alltimebinned)), 0.01)

    model = BoxLeastSquares(alltimebinned, allfluxbinned)
    results = model.power(periods, durations)


    index = np.argmax(results.power)
    period = results.period[index]
    t0 = results.transit_time[index]
    duration = results.duration[index]
    mid_transit_t0 = model.compute_stats(period, duration, t0)['transit_times'][0]

    # call the first round of plotting


    plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model, results, period, duration, t0, mid_transit_t0, args)

    stats_period = period
    stats_t0 = mid_transit_t0
    stats_duration = duration
    stats_depth = model.compute_stats(period, duration, t0)['depth']
    stats_depth_phased = model.compute_stats(period, duration, t0)['depth_phased']
    stats_depth_half = model.compute_stats(period, duration, t0)['depth_half']
    stats_depth_odd = model.compute_stats(period, duration, t0)['depth_odd']
    stats_depth_even = model.compute_stats(period, duration, t0)['depth_even']


    if (2*duration) >= period: # if the 'found' events are very short period, don't rn the BLS twice as the code would crash.
        return [stats_period, stats_t0, stats_duration, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [-999]

    # Find the in-transit points using a longer duration as a buffer to avoid ingress and egress
    in_transit = model.transit_mask(alltimebinned, period, 2*duration, t0)
    in_transit_notbinned = model.transit_mask(alltime, period, 2*duration, t0)

    # Re-run the algorithm, and plot the results
    model2 = BoxLeastSquares(alltimebinned[~in_transit], allfluxbinned[~in_transit])
    results2 = model2.power(periods, durations)

    # Extract the parameters of the best-fit model
    index = np.argmax(results2.power)
    period2 = results2.period[index]
    t02 = results2.transit_time[index]
    duration2 = results2.duration[index]
    mid_transit_t02 = model.compute_stats(period2, duration2, t02)['transit_times'][0]

    # call the second round of plotting - once the intitial transit has been removed
    plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model2, results2,period2,duration2,t02, mid_transit_t02, args, in_transit = in_transit, in_transit_notbinned = in_transit_notbinned)

    stats2_period = period2
    stats2_t0 = mid_transit_t02
    stats2_duration = duration2
    stats2_depth = model2.compute_stats(period2, duration2, t0)['depth']
    stats2_depth_phased = model2.compute_stats(period2, duration2, t0)['depth_phased']
    stats2_depth_half = model2.compute_stats(period2, duration2, t0)['depth_half']
    stats2_depth_odd = model2.compute_stats(period2, duration2, t0)['depth_odd']
    stats2_depth_even = model2.compute_stats(period2, duration2, t0)['depth_even']

    return [stats_period, stats_t0, stats_duration, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [stats2_period, stats2_t0, stats2_duration, stats2_depth, stats2_depth_phased, stats2_depth_half, stats2_depth_odd, stats2_depth_even]

def data_bls_FFI(tic, indir, alltime, allflux, args):
    '''
    function that rund the BLS routine for the FFI images. This code is differs from data_bls (which is used for the atrget pixel files) in as it doesn't have binned data.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times for the LC of the target
    allflux  :  list
        normalized flux for the LC of the target (not binned)

    Returns
    -------
        two lists of the statistics of the to BLS runs. Each list contains:
    stats_period
    stats_t0
    stats_duration
    stats_depth
    stats_depth_phased
    stats_depth_half
    stats_depth_odd
    stats_depth_even
    '''
    # make sure that there are no nan value sin the data - they cause everything to crash

    mask = np.isfinite(alltime) * np.isfinite(allflux)

    alltime = np.array(alltime)[mask]
    allflux = np.array(allflux)[mask]

    #durations = np.linspace(0.05, 0.2, 10)
    #model = BoxLeastSquares(alltime, allflux)
    #results = model.autopower(durations, frequency_factor=5.0)

    # define the durations and the period to test. The period ranges from 0.5 days to the length of the available data
    durations = np.linspace(0.05, 0.2, 10)
    periods = np.arange(0.5, (np.nanmax(alltime) - np.nanmin(alltime)), 0.01)

    model = BoxLeastSquares(alltime, allflux)
    #results = model.autopower(durations, frequency_factor=5.0)
    results = model.power(periods, durations)


    index = np.argmax(results.power)
    period = results.period[index]
    t0 = results.transit_time[index]
    duration = results.duration[index]
    mid_transit_t0 = model.compute_stats(period, duration, t0)['transit_times'][0]

    # call the first round of plotting
    plot_bls_FFI(tic, indir, alltime, allflux, model, results, period, duration, t0, mid_transit_t0, args)

    stats_period = period
    stats_t0 = t0
    stats_duration = duration
    stats_depth = model.compute_stats(period, duration, t0)['depth']
    stats_depth_phased = model.compute_stats(period, duration, t0)['depth_phased']
    stats_depth_half = model.compute_stats(period, duration, t0)['depth_half']
    stats_depth_odd = model.compute_stats(period, duration, t0)['depth_odd']
    stats_depth_even = model.compute_stats(period, duration, t0)['depth_even']


    if (2*duration) >= period: # if the 'found' events are very short period, don't rn the BLS twice as the code would crash.
        return [stats_period, stats_t0, stats_duration, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [-999]


    # Find the in-transit points using a longer duration as a buffer to avoid ingress and egress
    in_transit = model.transit_mask(alltime, period, 2*duration, t0)
    # Re-run the algorithm, and plot the results
    model2 = BoxLeastSquares(alltime[~in_transit], allflux[~in_transit])
    results2 = model2.power(periods, durations)

    # Extract the parameters of the best-fit model
    index = np.argmax(results2.power)
    period2 = results2.period[index]
    t02 = results2.transit_time[index]
    duration2 = results2.duration[index]
    mid_transit_t02 = model.compute_stats(period2, duration2, t02)['transit_times'][0]

    # call the second round of plotting - once the intitial transit has been removed
    plot_bls_FFI(tic, indir, alltime, allflux, model2, results2,period2,duration2,t02, mid_transit_t02, args, in_transit = in_transit)

    stats2_period = period2
    stats2_t0 = t02
    stats2_duration = duration2
    stats2_depth = model2.compute_stats(period2, duration2, t0)['depth']
    stats2_depth_phased = model2.compute_stats(period2, duration2, t0)['depth_phased']
    stats2_depth_half = model2.compute_stats(period2, duration2, t0)['depth_half']
    stats2_depth_odd = model2.compute_stats(period2, duration2, t0)['depth_odd']
    stats2_depth_even = model2.compute_stats(period2, duration2, t0)['depth_even']

    return [stats_period, stats_t0, stats_duration, stats_depth, stats_depth_phased, stats_depth_half, stats_depth_odd, stats_depth_even], [stats2_period, stats2_t0, stats2_duration, stats2_depth, stats2_depth_phased, stats2_depth_half, stats2_depth_odd, stats2_depth_even]

# --------------------------------------------
#                    plots                   #
# --------------------------------------------

# plot nearest neighbour LCs
def plot_nn(tic, indir, alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, transit_list, outtics, tessmag_list, distance, args):

    '''
    Plot the lighcurves of the 6 nearest neighbours to the target.

    Distance to the targets is calculated using the RA and DEC of the targte and adding them in quadrature

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime_nn  :  list
        times for all the nearest neighbours (not binned)
    allflux_nn  :  list
        normalized flux times for all the nearest neighbours (not binned)
    alltimebinned_nn  :  list
        binned time for all the nearest neighbours
    allfluxbinned_nn  :  list
        normalized binned flux for all the nearest neighbours
    transit_list  :  list
        list of all the marked transits
    outtics  :  list
        the tic IDs of the 6 nearest neighbours

    Returns
    -------
        Plot of the 6 nearest neighbours. The vertical line indicated the location(s) of the marked transit.

    '''

    fig, ax = plt.subplots(len(alltime_nn), 1, figsize=(13,8), sharex=True, gridspec_kw={'hspace': 0})
    plt.tight_layout()

    colors = ['r', 'darkorange', 'gold', 'seagreen', 'royalblue', 'navy','magenta' ]
    colors2 = ['k', 'k', 'k', 'k', 'k', 'grey','k' ]

    for i in range(0,len(alltime_nn)):

        for line in (transit_list):
            ax[i].axvline(line, color = 'k', linewidth = 2.2, alpha = 1, linestyle = '-')
        if str(outtics[i]) == str(tic):
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "*** {}  Tmag = {:3f} ***".format(tic, tessmag_list[i]), marker = '.', ms = 2, linewidth = 0)

        else:
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "{}  Tmag = {:3f}   d = {:3f} arcmins".format(outtics[i], tessmag_list[i], distance[i]), marker = '.', ms = 2, linewidth = 0)

        ax[i].plot(alltimebinned_nn[i], np.array(allfluxbinned_nn[i]), color = colors2[i], marker = '.', ms = 1, linewidth = 0)

        ax[i].legend(loc = 1)

    ax[0].set_title("LCs of Nearby Stars", fontsize = 13)
    ax[0].set_title("LCs of Nearby Stars", fontsize = 13)
    ax[len(alltime_nn) - 1].set_xlabel("Time (BJD-2457000)", fontsize = 13)
    ax[int(len(alltime_nn)/2)].set_ylabel("Normalised Flux", fontsize = 13)
    ax[0].set_xlim(np.nanmin(alltime_nn), np.nanmax(alltime_nn))

    plt.xlim(np.nanmin(alltime_nn), np.nanmax(alltime_nn))

    if args.save == True:
        plt.savefig('{}/{}/{}_nearest_neighbours.png'.format(indir, tic, tic), format='png', bbox_inches='tight')


    if args.noshow == False:
        plt.show()
    else:
        plt.close()

# plot the image cut-out of the nearby stars
def plot_cutout(image):
    """
    Plot image cut out of the target.
    """
    plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r,
           vmax = np.percentile(image, 92),
           vmin = np.percentile(image, 5))

    plt.grid(axis = 'both',color = 'white', ls = 'solid')

# plot the centroid positions and whether thet move in transit (only zoomed in on transit events and not for the whole LC)
def plot_centroid(tic, indir, alltime12, allx1, ally1, allx2, ally2, transit_list, args, ql = False):
    '''
    Plot the x and y centroids around the time(s) of the marked transit.

    Download the LCs for the target star for all the indicated sectors

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime12  :  list
        time used for the x and y centroid position plotting
    allx1  :  list
        CCD column position of target’s flux-weighted centroid. In x direction
    allx2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In x direction
    ally1  :  list
        CCD column position of target’s flux-weighted centroid. In y direction
    ally2  :  list
        The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects. In y direction
    transit_list   :  list
        list of the marked transit events

    Returns
    -------
        Plot of the centroid positions in the x and y direction.
        The black points are the centroid position - moves if this is a blend.

    '''


    gs = len(transit_list) # the grid size

    centroid1 = [allx1, ally1]  # CCD column position of target’s flux-weighted centroid.
    centroid2 = [allx2, ally2]  # The CCD column local motion differential velocity aberration (DVA), pointing drift, and thermal effects.

    if gs == 1:
        plt.figure(figsize=(7, 7))
    else:
        plt.figure(figsize=(5 * gs, 7))

    for g,peak in enumerate(transit_list):

        for i,cen1 in enumerate(centroid1):

            mask = (np.array(alltime12) < peak+2) & (np.array(alltime12) > peak-2)

            # ------ ALL X --------
            ax1 = plt.subplot(2,gs,g+1) # define the plotting area

            minf = np.nanmin(np.concatenate((np.array(centroid1[0])[mask],np.array(centroid2[0])[mask]), axis=0) )
            maxf = np.nanmax(np.concatenate((np.array(centroid1[0])[mask],np.array(centroid2[0])[mask]), axis=0) )

            #maxf = np.nanmax(np.array(centroid1[0])[mask])
            height = maxf - minf
            if i == 0:
                ax1.plot(np.array(alltime12)[mask],np.array(centroid2[0])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='Local Motion')
                ax1.plot(np.array(alltime12)[mask],np.array(centroid1[0])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='Flux-weighted Centroid')
            else:
                ax1.plot(np.array(alltime12)[mask],np.array(centroid2[0])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='_nolegend_')
                ax1.plot(np.array(alltime12)[mask],np.array(centroid1[0])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='_nolegend_')

            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax1.get_yticklabels(), visible=False)


            plt.tight_layout()
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = main_pot_color)
            plt.ylim(minf,maxf)
            if g == 0:
                plt.legend()
                plt.ylabel('x centroid shift', fontsize = main_plot_fontsize)
            plt.title('Transit {}'.format(g+1))

            # ------ ALL Y --------
            ay1 = plt.subplot(2,gs,gs+1+g, sharex = ax1) # define the plotting area
            minf = np.nanmin(np.concatenate((np.array(centroid1[1])[mask],np.array(centroid2[1])[mask]), axis=0) )
            maxf = np.nanmax(np.concatenate((np.array(centroid1[1])[mask],np.array(centroid2[1])[mask]), axis=0) )

            height = maxf - minf

            if i == 0:
                ay1.plot(np.array(alltime12)[mask],np.array(centroid2[1])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='Local Motion')
                ay1.plot(np.array(alltime12)[mask],np.array(centroid1[1])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='Flux-weighted Centroid')
            else:
                ay1.plot(np.array(alltime12)[mask],np.array(centroid2[1])[mask],color ='r',marker = 'o', lw = 0, ms=3,alpha=1,label='_nolegend_')
                ay1.plot(np.array(alltime12)[mask],np.array(centroid1[1])[mask],color ='k',marker = 'o', lw = 0, ms=2,alpha=0.5,label='_nolegend_')


            plt.setp(ay1.get_yticklabels(), visible=False)
            ay1.axes.yaxis.set_ticks([])
            ax1.axes.yaxis.set_ticks([])

            plt.subplots_adjust(wspace=0.02, hspace=0.04)
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = main_pot_color)
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)', fontsize = main_plot_fontsize)
            if g == 0:
                plt.ylabel('y centroid shift', fontsize = main_plot_fontsize)


    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_centroids.png'.format(indir, tic, tic), format='png')

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()


# plot the LCs in two different aperture sizes
def plot_aperturesize(tic, indir, TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, transit_list, args, ql = False):
    '''correspond to: tic,  indir,      alltime,         alltime,         alltime,      allflux_normal, allflux_normal, allflux_small,   transit_list

    LC plot around the time of transit-event extracted in two different aperture sizes. The LC is not corrected for any systematics.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
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
    transit_list   :  list
        list of the marked transit events

    Returns
    -------
        plot of the transit-event extracted in two different aperture sizes. EBs will exhibit different transit shapes with the different aperture sizes.
        The time of the transit event is indicated by the vertical line.

    '''
    if args.FFI == True:
        frame_width = 1.5
    else:
        frame_width = 1.5


    gs = len(transit_list)

    if gs == 1:

        plt.figure(figsize=(7,4))
        plt.tight_layout()
        for g,peak in enumerate(transit_list):

            mask_unb = (np.array(TESS_unbinned_t_l) < peak+frame_width) & (np.array(TESS_unbinned_t_l) > peak-frame_width)
            mask_bin = (np.array(TESS_binned_t_l) < peak+frame_width) & (np.array(TESS_binned_t_l) > peak-frame_width)
            mask_small = (np.array(small_binned_t_l) < peak+frame_width) & (np.array(small_binned_t_l) > peak-frame_width)

            if args.FFI == False:
                minf = np.nanmin(np.array(TESS_unbinned_l)[mask_unb])
                maxf = np.nanmax(np.array(TESS_unbinned_l)[mask_unb])
            else:
                minf = np.nanmin(np.array(TESS_binned_l)[mask_bin])
                maxf = np.nanmax(np.array(TESS_binned_l)[mask_bin])
                diff = maxf - minf
                minf = minf - (diff * 0.01)
                maxf = maxf + (diff * 0.01)

            if args.FFI == False:
                plt.scatter(TESS_unbinned_t_l, TESS_unbinned_l, s = 3, marker = 's',alpha = 0.4, color = 'black', label = 'TESS unbinned')

            plt.scatter(TESS_binned_t_l, TESS_binned_l, s = 11,  marker = 'o', alpha = 1, color = 'blue', label = 'TESS ap')
            plt.scatter(small_binned_t_l, small_binned_l, s = 12, marker = '>', alpha =1, color = 'red', label = 'Small ap')

            plt.title("Detrended LC with various aperture sizes for TIC {}".format(tic), fontsize = main_plot_fontsize)
            plt.tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
            plt.tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
            plt.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')

            #plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
            plt.xlim(peak-frame_width, peak+frame_width)
            plt.axvline(peak, color = main_pot_color, linestyle = '--')
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.ylabel('Normalized Flux')
            plt.title('Aperture Size Test, Transit {}'.format(g+1), fontsize = main_plot_fontsize)


        for transit in transit_list:
            plt.axvline(transit, color = main_pot_color, linestyle = '--', linewidth = 2)

        plt.legend(fontsize = main_plot_fontsize)

        if (args.save == True) and (ql == False):
            plt.subplots_adjust(bottom = 0.13, right = 0.97, top = 0.94)
            plt.savefig('{}/{}/{}_aperture_size.png'.format(indir, tic, tic), format='png')

        if (args.noshow == False) or (ql == True):
            plt.show()
        else:
            plt.close()

    else:
        plt.figure(figsize=(gs*6,4))
        for g,peak in enumerate(transit_list):


            mask_unb = (np.array(TESS_unbinned_t_l) < peak+frame_width) & (np.array(TESS_unbinned_t_l) > peak-frame_width)
            mask_bin = (np.array(TESS_binned_t_l) < peak+frame_width) & (np.array(TESS_binned_t_l) > peak-frame_width)
            mask_small = (np.array(small_binned_t_l) < peak+frame_width) & (np.array(small_binned_t_l) > peak-frame_width)

            if np.sum(mask_unb) != 0:

                if args.FFI == False:
                    minf = np.nanmin(np.array(TESS_unbinned_l)[mask_unb])
                    maxf = np.nanmax(np.array(TESS_unbinned_l)[mask_unb])
                else:
                    minf = np.nanmin(np.array(TESS_binned_l)[mask_bin])
                    maxf = np.nanmax(np.array(TESS_binned_l)[mask_bin])
                    diff = maxf - minf
                    minf = minf - (diff * 0.01)
                    maxf = maxf + (diff * 0.01)


                plt.subplot(1,gs,g+1)

                if args.FFI == False:
                    plt.scatter(TESS_unbinned_t_l, TESS_unbinned_l, s = 3, marker = 's', alpha = 0.4, color = 'black', label = 'TESS Aperture unbinned')
                plt.scatter(TESS_binned_t_l, TESS_binned_l, s = 11,  marker = 'o', alpha = 1, color = 'blue', label = 'TESS Aperture binned')
                plt.scatter(small_binned_t_l, small_binned_l, s = 12, marker = '>', alpha =1, color = 'red', label = 'Small Aperture binned')

                if g > 0:
                    plt.tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
                    plt.yticks([])

                plt.tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
                plt.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')

                plt.tight_layout()
                plt.xlim(peak-frame_width, peak+frame_width)
                plt.axvline(peak, color = main_pot_color, linestyle = '--')
                plt.ylim(minf,maxf)
                plt.xlabel('Time (BJD-2457000)')
                plt.title('Aperture Size Test, Transit {}'.format(g+1), fontsize = main_plot_fontsize)


        if (args.save == True) and (ql == False):
            plt.subplots_adjust(bottom = 0.13, right = 0.97, top = 0.94)
            plt.savefig('{}/{}/{}_aperture_size.png'.format(indir, tic, tic), format='png')

        if (args.noshow == False) or (ql == True):
            plt.show()
        else:
            plt.close()


# plot the background flux around the time of the transit like event.
def plot_background(tic, indir, alltime, allfbkg, transit_list, args, ql = False):

    '''
    LC of the bakcground flux at the time of the transit event.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allfbkg  :  list
        background flux
    transit_list   :  list
        list of the marked transit events

    Returns
    -------
        plot of the background flux at the time of the marked transit-events.
        The backrgound flux should not exhibit any 'spikes' or 'unusual' behaviour at the time of the transit event.
        The time of the transit event is indicated by the vertical line.

    '''


    gs = len(transit_list) # the grid size

    if len(transit_list) == 1:
        plt.figure(figsize=(7,4))
        plt.tight_layout()

        for g,peak in enumerate(transit_list):

            mask_dd = (np.array(alltime) < peak+2) & (np.array(alltime) > peak-2)

            minf = np.nanmin(np.array(allfbkg)[mask_dd])
            maxf = np.nanmax(np.array(allfbkg)[mask_dd])
            height = maxf - minf

            plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
            plt.xlim(peak-2, peak+2)
            plt.axvline(peak, color = main_pot_color)
            plt.ylim(minf,maxf)
            plt.xlabel('Time (BJD-2457000)')
            plt.ylabel('Flux')
            plt.title('Background Flux, Transit {}'.format(g+1))

        if (args.save == True) and (ql == False):
            plt.subplots_adjust(bottom = 0.13, right = 0.97, top = 0.94)
            plt.savefig('{}/{}/{}_background.png'.format(indir, tic, tic), format='png')

        if (args.noshow == False) or (ql == True):
            plt.show()
        else:
            plt.close()

    else:

        gs = len(transit_list) # the grid size
        plt.figure(figsize=(6 * gs, 4))

        for g,peak in enumerate(transit_list):

            mask_dd = (np.array(alltime) < peak+2) & (np.array(alltime) > peak-2)

            if np.sum(mask_dd) != 0:

                minf = np.nanmin(np.array(allfbkg)[mask_dd])
                maxf = np.nanmax(np.array(allfbkg)[mask_dd])
                height = maxf - minf

                plt.subplot(1,gs,g+1)
                plt.plot(np.array(alltime)[mask_dd], np.array(allfbkg)[mask_dd], 'o', markersize = 2, color = 'blue', alpha = 0.7,label='centroid', markerfacecolor='white')
                plt.xlim(peak-2, peak+2)
                plt.axvline(peak, color = main_pot_color)

                try:
                    plt.ylim(minf,maxf)
                except:
                    print('axis limit error (centroid positions)')

                plt.tight_layout()
                plt.xlabel('Time (BJD-2457000)')
                plt.ylabel('Flux')
                plt.title('Background Flux, Transit {}'.format(g+1))

        if (args.save == True) and (ql == False):
            plt.subplots_adjust(bottom = 0.13, right = 0.97, top = 0.94)
            plt.savefig('{}/{}/{}_background.png'.format(indir, tic, tic), format='png', bbox_inches='tight')

        if (args.noshow == False) or (ql == True):
            plt.show()
        else:
            plt.close()

# plot the nearby TESS stars as well as the SDSS cutout - reprojected to have North up.
def plot_TESS_stars(tic,indir, transit_sec, tpf_list, args):

    '''
    Plot of the field of view round the target star showing nearby stars that are brighter than magnitude 17 as well as the SDSS cutout.
    Both images are projected and oriented North for easy comparison.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_sec  :  list or str
        list of the sectors that are being analyse.
    tpf_list   : list
        list of the target pixel files (for each sector)

    Returns
    -------
        Plot of the averaged flux per pixel around the target (left) as well as the SDSS plot (right). The red star on the right plot indicated the location of the target.
        The orange circles show the location of nearby stars with magnitudes brighter than 17 mag where their relative sizes correspond to their relative brightness.
        The location of the target star is hown with the reticle on the right hans side SDSS image.

    Tmag   :   float
        TESS magnitude of the target star
    Teff   :   float
        Effective temperature of the target star (K)
    rad   :   float
        radius of the target star (Solar radii)
    mass   :   float
        mass of the target star (Solar masses)
    '''

    # ----------
    # import the astroplan module that is needed - this is done here and not at the start at this scipt
    # because astroplan cannot be parallelised (issues with python's shelve storage) so import here.
    from astroplan import FixedTarget
    from astroplan.plots import plot_finder_image
    # ----------

    # Query nearby Gaia Stars  --------------------

    sector = str(transit_sec[0])

    starName = "TIC " + str(tic)
    radSearch = 5/60 #radius in degrees

    # this function depends on astroquery working, and sometimes it doesn't.
    # for when it doesn't work (or simply can't connect to it), just skip plotting the other TESS stars.
    try:
        catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
    except:
        print ("Currently cannot connect to Astroquery.")
        # return values that we know aren't real so that we can tell the code that the plotting didn't work
        return -999, -999, -999, 1, -999,-999,-999,-999

    # ra and dec of the target star
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']

    # while we have the astroquery loaded, let's collect some other information about the star
    # these paramaters can help us find out what type of star we have with just a glance

    vmag = catalogData['Vmag'][0] # v magnitude (this migth be more useful than the TESS mag for things such as osbevring)
    logg = catalogData['logg'][0] # logg of the star
    mass = catalogData['mass'][0] # mass of the star
    plx = catalogData['plx'][0]   # parallax

    # sometimes these values aren't accessible through astroquery - so we shoudl just quickly check.
    if not np.isfinite(vmag): vmag = '--' # this is what will appear in the table of the report to indicate that it's unknown
    if not np.isfinite(logg): logg = '--'
    if not np.isfinite(mass): mass = '--'
    if not np.isfinite(plx): plx   = '--'

    # sometimes it's useufl to know if the star has another name
    # check whether it was osberved by one of these four large surveys

    catalogs = ['HIP', 'TYC', 'TWOMASS', 'GAIA']

    for cat in catalogs:
        c_id = str(catalogData[0][cat])
        if c_id != '--':
            c_id = "{} {}".format(cat,c_id)
            break
        else:
            continue

    # ------------------------------------------

    # Create a list of nearby bright stars (tess magnitude less than 17) from the rest of the data for later.
    bright = catalogData['Tmag'] < catalogData['Tmag'][0] + 5

    # ---------------------------------------------
    # Get the data for the SDSS sky viewer --------

    survey = 'DSS2 Red'
    fig, ax = plt.subplots()
    plt.axis("off")

    target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

    try:
        target = FixedTarget(coord=target_coord, name="Survey = {}".format(survey))
        ax, hdu = plot_finder_image(target, survey = survey, reticle='True', fov_radius=5*u.arcmin)

    except: # if DSS2 Red is not available, download the DSS field of view image instead
        try:
            survey = 'DSS'
            target = FixedTarget(coord=target_coord, name="Survey = {}".format(survey))
            ax, hdu = plot_finder_image(target, survey = survey, reticle='True', fov_radius=5*u.arcmin)
        except:
            return -111,-111,-111,-111,-111,-111,-111,-111

    plt.close('all')

    # --------------------------------------------
    # only run this for the first tpf (just needs to be one and might aswell take the first one)

    tpf =  tpf_list[0]

    fig= plt.figure(figsize=(7,5.5))

    #sector =  tpf.header['SECTOR']
    sector = tpf[0].header['SECTOR']

    plt.title('Sector {}'.format(sector))

    if (args.FFI == False) or (args.FFI == 'SPOC'):
        tpf_wcs = extract_wcs(tpf)
        tup = (np.nanmean(tpf[1].data['FLUX'], axis=0),tpf_wcs)
    else:
        tpf_wcs = tpf.wcs
        tup = (np.nanmean(tpf.flux, axis=0),tpf_wcs)
    # create a tupple of the array of the data and the wcs projection of the TESS cutout


    # map the SDSS and TESS image onto each other - the output will be orented NORTH!
    wcs_out, shape_out = find_optimal_celestial_wcs(input_data =[tup, hdu])


    # plot the reprojected TESS image
    ax1 = plt.subplot(1,2,1, projection=wcs_out)
    array, footprint = reproject_interp(tup, wcs_out,shape_out = shape_out,order = 'nearest-neighbor')

    ax1.imshow(array, origin='lower', cmap = plt.cm.YlGnBu_r)

    ax1.coords['ra'].set_axislabel('Right Ascension', fontsize = main_plot_fontsize)
    ax1.coords['dec'].set_axislabel('Declination', fontsize = main_plot_fontsize)
    ax1.grid(color = 'grey', alpha = 0.7)

    # plot the nearby GAIA stars on this image too...
    ra_stars, dec_stars = catalogData[bright]['ra'], catalogData[bright]['dec']
    s = np.maximum((19 - catalogData[bright]['Tmag'])*5, 0)  # the size corresponds to their brightness
    ax1.scatter(ra_stars, dec_stars, s=s, transform=ax1.get_transform('icrs'), color=main_pot_color, zorder=100)

    # plot the target star that we're looking at
    ax1.scatter(ra, dec, s= 200, transform=ax1.get_transform('icrs'), marker = '*', color='red', zorder=100)
    ax1.tick_params(labelsize=main_plot_fontsize)

    # plot the reprojected SDSS image
    ax2 = plt.subplot(1,2,2, projection=wcs_out, sharex=ax1, sharey=ax1)
    array, footprint = reproject_interp(tup, wcs_out,shape_out = shape_out)
    ax2.imshow(hdu.data, origin='lower', cmap = 'Greys')
    ax2.coords['ra'].set_axislabel('Right Ascension', fontsize = main_plot_fontsize)
    #ax2.coords['dec'].set_axislabel('Declination')

    # Draw reticle ontop of the target star
    pixel_width = hdu.data.shape[0]
    inner, outer = 0.03, 0.08

    reticle_style_kwargs = {}
    reticle_style_kwargs.setdefault('linewidth', 1.5)
    reticle_style_kwargs.setdefault('color', 'red')

    ax2.axvline(x=0.5*pixel_width, ymin=0.5+inner, ymax=0.5+outer,
               **reticle_style_kwargs)
    ax2.axvline(x=0.5*pixel_width, ymin=0.5-inner, ymax=0.5-outer,
               **reticle_style_kwargs)
    ax2.axhline(y=0.5*pixel_width, xmin=0.5+inner, xmax=0.5+outer,
               **reticle_style_kwargs)
    ax2.axhline(y=0.5*pixel_width, xmin=0.5-inner, xmax=0.5-outer,
                   **reticle_style_kwargs)
    ax2.grid()
    ax2.tick_params(labelsize=main_plot_fontsize)
    plt.tight_layout(w_pad= 7)

    if args.save == True:
        plt.savefig('{}/{}/{}_star_field.png'.format(indir, tic, tic), format='png', bbox_inches='tight')

    if args.noshow == False:
        plt.show()
    else:
        plt.close()

    return catalogData['Tmag'][0], catalogData['Teff'][0], catalogData['rad'][0], mass, vmag, logg, plx, c_id


# same as plot_TESS_stars but not re-projected
def plot_TESS_stars_not_proj(tic, indir, transit_list, transit_sec, tpf_list, args):

    '''
    Plot of the field of view round the target star showing nearby stars that are brighter than magnitude 15.

    This function does not use astroplan (so no SDSS image) meaning that it needs to be used if the code is parallelised
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_sec  :  list or str
        list of the sectors that are being analyse.
    tpf_list   : list
        list of the target pixel files (for each sector)

    Returns
    -------
        Plot of the averaged flux per pixel around the target. The red star indicated the location of the target.
        The orange circles show the location of nearby stars with magnitudes brighter than 15 mag.

    '''

    # always check whether the file already exists... as to not waste computer power and time
    sector = str(transit_sec[0])

    starName = "TIC " + str(tic)
    radSearch = 4/60 #radius in degrees

    catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")

    # ra and dec of the target star
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']

    # while we have the astroquery loaded, let's collect some other information about the star
    # these paramaters can help us find out what type of star we have with just a glance

    vmag = catalogData['Vmag'][0] # v magnitude (this migth be more useful than the TESS mag for things such as osbevring)
    logg = catalogData['logg'][0] # logg of the star
    mass = catalogData['mass'][0] # mass of the star
    plx = catalogData['plx'][0]   # parallax

    # sometimes these values aren't accessible through astroquery - so we shoudl just quickly check.
    if not np.isfinite(vmag): vmag = '--' # this is what will appear in the table of the report to indicate that it's unknown
    if not np.isfinite(logg): logg = '--'
    if not np.isfinite(mass): mass = '--'
    if not np.isfinite(plx): plx   = '--'
    # sometimes it's useufl to know if the star has another name
    # check whether it was osberved by one of these four large surveys

    catalogs = ['HIP', 'TYC', 'TWOMASS', 'GAIA']

    for cat in catalogs:
        c_id = str(catalogData[0][cat])
        if c_id != '--':
            cat_id = "{} {}".format(cat,c_id)
            break
        else:
            continue

    # ------------------------------------------


    # Create a list of nearby bright stars (tess magnitude less than 14) from the rest of the data for later.
    bright = catalogData['Tmag'] <  catalogData['Tmag'][0] + 5

    start = [np.float64(transit_list[0]) - 0.2]
    end = [np.float64(transit_list[0]) + 0.2]

    # background is just the array of the flux of all the pixels (used for backrgoudn in pixel level LC plot so mildly confusing and should change)

    if (args.FFI == False) or (args.FFI == 'SPOC'): # SC data mode

        for i, tpf in enumerate(tpf_list):

            if (start > np.nanmin(tpf[1].data['TIME']) and start < np.nanmax(tpf[1].data['TIME'])):
                fig, ax = plt.subplots(figsize=(5,5))
                plt.tight_layout()

                sector = tpf[0].header['SECTOR']

                tpf_wcs = extract_wcs(tpf)

                ax = plt.subplot(projection=tpf_wcs)

                plot_cutout(tpf[1].data['FLUX'][0])

                ra_stars, dec_stars = catalogData[bright]['ra'], catalogData[bright]['dec']
                s = np.maximum((19 - catalogData[bright]['Tmag'])*5, 0)
                plt.scatter(ra_stars, dec_stars, s=s, transform=ax.get_transform('icrs'), color=main_pot_color, zorder=100)

                # plot the target
                plt.scatter(ra, dec, s= 200, transform=ax.get_transform('icrs'), marker = '*', color='red', zorder=100)

                plt.title("Sector {}".format(sector), fontsize = main_plot_fontsize)

                plt.xlim(-0.5,10.5)
                plt.ylim(-0.5,10.5)

                #plt.xlabel("RA")
                #plt.ylabel("Dec")

                plt.tight_layout(h_pad= 0.5)

                if args.save == True:
                    plt.subplots_adjust(bottom = 0.1, left = 0.17, top = 0.94)
                    plt.savefig('{}/{}/{}_star_field.png'.format(indir, tic, tic), format='png')

                if args.noshow == False:
                    plt.show()
                else:
                    plt.close()

    else: #FFI mode

        for i, tpf in enumerate(tpf_list):

            if (start > np.nanmin(tpf.time) and start < np.nanmax(tpf.time)):
                fig, ax = plt.subplots(figsize=(5,5))
                plt.tight_layout()

                sector = tpf.header['SECTOR']

                tpf_wcs = tpf.wcs

                ax = plt.subplot(projection=tpf_wcs)

                plot_cutout(np.nanmean(tpf.flux, axis = 0))

                ra_stars, dec_stars = catalogData[bright]['ra'], catalogData[bright]['dec']
                s = np.maximum((19 - catalogData[bright]['Tmag'])*5, 0)
                plt.scatter(ra_stars, dec_stars, s=s, transform=ax.get_transform('icrs'), color=main_pot_color, zorder=100)

                # plot the target
                plt.scatter(ra, dec, s= 200, transform=ax.get_transform('icrs'), marker = '*', color='red', zorder=100)

                plt.title("Sector {}".format(sector), fontsize = main_plot_fontsize)

                plt.xlim(-0.5,10.5)
                plt.ylim(-0.5,10.5)

                #plt.xlabel("RA")
                #plt.ylabel("Dec")

                plt.tight_layout(h_pad= 0.5)

                if args.save == True:
                    plt.subplots_adjust(bottom = 0.1, left = 0.17, top = 0.94)
                    plt.savefig('{}/{}/{}_star_field.png'.format(indir, tic, tic), format='png')

                if args.noshow == False:
                    plt.show()
                else:
                    plt.close()

    return catalogData['Tmag'][0], catalogData['Teff'][0], catalogData['rad'][0], mass, vmag, logg, plx, c_id

# LC per pixel
def plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, tpf_list, apmask_list, arrshape_list, t_list, transit_list, args, ql = False):

    '''
    Plot the LC for each pixel around the time of the transit like event. Each LC is fit with a spline and corrected to flatten.
    each LC is fitted with a 3 order polynomial in order to flatten.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
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
    transit_list  :  int
        list of all the marked transits

    Returns
    -------
        Plot of the normalised LC for each pixel around the time of the transit like event.
        The pixel backrgound colour represents the average flux.
        The time of the transit is highlighted in red/gold for each pixel LC.
    '''

    # loop through the transits and make plot for each ( only the first is currently displayed in the pdf report)
    for idx, X1 in enumerate(X1_list):


        X4 = X4_list[idx]
        oot = oot_list[idx]
        #intr = intr_list[n]
        bkg = np.flip(bkg_list[idx], axis = 0)
        arrshape = arrshape_list[idx]
        t = t_list[idx]
        peak = transit_list[idx]
        tpf = tpf_list[idx]

        if args.FFI != 'QLP':
            apmask = apmask_list[idx]

            mapimg = np.flip(apmask_list[idx], axis = 0)

            ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
            hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

        fig, ax = plt.subplots(arrshape[1], arrshape[2], sharex = True, sharey = False, gridspec_kw={'hspace': 0 ,'wspace': 0}, figsize=(8,8))

        plt.tight_layout()

        # see if the backrgound of this plot can be the average pixel flux (if there are too many nans this will fail and the background will just be black which is also okay)
        try:
            color = plt.cm.viridis(np.linspace(0, 1,int(np.nanmax(bkg))-int(np.nanmin(bkg))+1))
            simplebkg = False
        except:
            simplebkg = True

        for i in range(0,arrshape[1]):
            print ("{}   out of    {} ".format(i+1,arrshape[1] ))
            ii = arrshape[1]-1-i # we want to plot this such that the pixels increase from left to right and bottom to top

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

                if args.FFI == False:
                    binfac = 5

                    N       = len(time)
                    n       = int(np.floor(N/binfac)*binfac)
                    X       = np.zeros((2,n))
                    X[0,:]  = time[:n]
                    X[1,:]  = f1[:n]
                    Xb      = rebin(X, (2,int(n/binfac)))

                    # binned data
                    time_binned    =    np.array(Xb[0])
                    flux_binned    =   np.array(Xb[1])

                else:
                    # binned data -
                    time_binned    =    np.array(time)
                    flux_binned  =   np.array(flux)

                # create a mask that only looks at the times cut around the transit-event
                timemask = (time_binned < peak+1.5) & (time_binned > peak-1.5)

                time_binned = time_binned[timemask]
                flux_binned = flux_binned[timemask]

                # ----------
                # fit a spline to the cut-out of each pixel LC in order to flatten it
                p = np.poly1d(np.polyfit(time_binned, flux_binned, 3))
                flux_binned = flux_binned/p(time_binned)
                # ----------

                intr = abs(peak-time_binned) < 0.1

                if simplebkg == True:
                    ax[ii, j].set_facecolor(color = 'k')
                    linecolor = 'w'
                    transitcolor = 'gold'
                else:
                    ax[ii, j].set_facecolor(color = color[int(bkg[ii,j])-int(np.nanmin(bkg))])

                    if int(bkg[ii,j])-abs(int(np.nanmin(bkg))) > ((np.nanmax(bkg))-abs(int(np.nanmin(bkg))))/2:
                        linecolor = 'k'
                        transitcolor = 'orangered'
                    else:
                        linecolor = 'w'
                        transitcolor = 'gold'


                ax[ii, j].plot(time_binned,flux_binned, color = linecolor, marker = '.', markersize=1, lw = 0)
                ax[ii, j].plot(time_binned[intr],flux_binned[intr], color = transitcolor, marker = '.', markersize=1, lw = 0)

                # get rid of ticks and ticklabels
                ax[ii,j].set_yticklabels([])
                ax[ii,j].set_xticklabels([])
                ax[ii,j].set_xticks([])
                ax[ii,j].set_yticks([])

        # ------------------
        if args.FFI != 'QLP':
            print ("\n Calculating the Aperture Mask...", end =" ")

            for i in range(0,len(ver_seg[1])):
                ax[ver_seg[0][i], ver_seg[1][i]].spines['right'].set_color('red')
                ax[ver_seg[0][i], ver_seg[1][i]].spines['right'].set_linewidth(6)

                top = (ax[ver_seg[0][i], ver_seg[1][i]].get_ylim()[1])
                bottom= (ax[ver_seg[0][i], ver_seg[1][i]].get_ylim()[0])

                change = top - bottom

                ax[ver_seg[0][i], ver_seg[1][i]].spines['right'].set_bounds(bottom,top - (change*0.08) )

                for j in range(0,len(hor_seg[1])):
                    ax[hor_seg[0][j], hor_seg[1][j]].spines['bottom'].set_color('red')
                    ax[hor_seg[0][j], hor_seg[1][j]].spines['bottom'].set_linewidth(6)
                    ax[hor_seg[0][j], hor_seg[1][j]].spines['bottom'].set_bounds(peak-1.5,peak+1.3)

        print ("done.\n")
        # ------------------

        # label the pixels
        if (args.FFI == False) or (args.FFI == 'SPOC'):
            start_column = tpf[1].header['1CRV5P']
            start_row = tpf[1].header['2CRV5P']

        else:
            start_column = tpf.get_keyword('1CRV5P', hdu=1, default=0)
            start_row = tpf.get_keyword('2CRV5P', hdu=1, default=0)

        y_start = '{}'.format(start_row)
        y_end = '{}'.format(start_row + bkg.shape[0])

        x_start = '{}'.format(start_column)
        x_end = '{}'.format(start_column + bkg.shape[1])

        ax[bkg.shape[0]-1, 0].set_xlabel(x_start)
        ax[bkg.shape[0]-1, 0].set_ylabel(y_start)

        ax[0,0].set_ylabel(y_end)
        ax[bkg.shape[0]-1,bkg.shape[1]-1].set_xlabel(x_end)

        fig.text(0.5,0.01, "column (pixel)", ha='center', fontsize = 13)
        fig.text(0.01, 0.5, "row (pixel)", va='center', rotation='vertical', fontsize = 13)

        # - - - - - - - - - -

        plt.subplots_adjust(top=0.95, right = 0.99, bottom = 0.04, left = 0.04)


        plt.suptitle(r"T0 = {} $\pm$ 1.5 d".format(peak ),y=0.98, fontsize = 15)
        plt.xlim(peak-1.5,peak+1.5)

        if (args.save == True) and (ql == False):
            print ("Waiting on plot...")
            plt.savefig('{}/{}/{}_individual_pixel_LCs_{}.png'.format(indir, tic,tic, idx), format='png')

        if (args.noshow == False) or (ql == True):
            plt.show()
        else:
            plt.close()


# full light curve with the momentum dumps
def plot_full_md_old(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, transit_list, args):
    '''
    Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)).

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    allflux_err  :  list
        normalized flux errors (not binned)
    all_md  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    transit_list  :  list
        list of all the marked transit-events

    Returns
    -------
        Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)).
        The plot is always saved as this funtion is only called if the 'save' option was chosen
    '''
    gs = len(transit_list)

    plt.figure(figsize=(15,10))
    plt.tight_layout()

    if (args.FFI == False) or (args.FFI == 'SPOC'):
        # rename things so that the code didn't have to be changed - not a very 'tidy' solution.
        time_dd = alltime
        flux_dd = allflux
        time_dd_binned = alltimebinned
        flux_dd_binned = allfluxbinned

    else:
        # --- do some sigma clipping to make the LC look better ----
        MAD = median_absolute_deviation(allflux,ignore_nan = True)
        madrange = (5 * MAD * 100)
        ymask = (allflux < 1 + madrange) * (allflux > 1 - madrange)

        time_dd = np.array(alltime)[ymask]
        flux_dd = np.array(allflux)[ymask]
        time_dd_binned = np.array(alltime)[ymask]
        flux_dd_binned = np.array(allflux)[ymask]
        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---

    line_dd = all_md

    minf_list = []
    maxf_list = []

    for peak in transit_list:

        mask_dd = (np.array(time_dd) < peak+cutout_window) & (np.array(time_dd) > peak-cutout_window)

        minf0 = np.nanmin(np.array(flux_dd)[mask_dd])
        maxf0 = np.nanmax(np.array(flux_dd)[mask_dd])

        minf_list.append(minf0)
        maxf_list.append(maxf0)

    print (transit_list, minf_list,maxf_list )
    minf = np.nanmin(minf_list)
    maxf = np.nanmax(maxf_list)

    height = maxf - minf

    for g,peak in enumerate(transit_list):

        mask_dd = (np.array(time_dd) < peak+cutout_window) & (np.array(time_dd) > peak-cutout_window)
        mask_dd_binned = (np.array(time_dd_binned) < peak+cutout_window) & (np.array(time_dd_binned) > peak-cutout_window)

        if np.sum(mask_dd) != 0:

            if gs == 1:
                plt.subplot(2,4,(6,7))

            else:
                plt.subplot(2,gs,(gs+g+1))

            height_cut = maxf - minf

            # plot the momentum dumps and the markings.
            # plot the lines first so that they are behind the data - we need to be able to see the data well...
            plt.vlines(line_dd, minf-(height_cut/10), minf + height_cut*0.25 , colors = md_color, label = "Momentum Dump", zorder=1)
            plt.vlines([peak], minf-(height_cut/10),minf + height*0.25 , linewidth = 4, colors = 'darkred', linestyle = '--', zorder=2, alpha = 0.85)


            plt.plot(np.array(time_dd)[mask_dd], np.array(flux_dd)[mask_dd], 'o', markersize = 4, color = main_pot_color, alpha = 0.8, label = "unbinned", markerfacecolor='white', zorder=3)
            plt.plot(np.array(time_dd_binned)[mask_dd_binned], np.array(flux_dd_binned)[mask_dd_binned], marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 5, label = 'binning = 7', markerfacecolor='k', zorder=4)

            plt.xlim(peak-cutout_window, peak+cutout_window)


            try:
                if individual_ylim == False:
                    plt.ylim(minf-(height_cut/10), maxf + (height_cut/10))
            except:
                print ('axis limits error (momentun dumps)')

            plt.title('Transit {}'.format(g+1))

        else:
            print ("size 0 so skip")


    if gs == 1:
        plt.subplot(2,4,(1,4))
    else:
        plt.subplot(2,gs,(1,gs))

    if args.FFI == False:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = main_pot_color, alpha = 0.9, label = "unbinned", markerfacecolor='white', zorder=2)
        plt.plot(np.array(time_dd_binned), np.array(flux_dd_binned), marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 1, label = 'binning = 7', markerfacecolor='k', zorder=3)
    else:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = '#054950', alpha = 0.9, label = "unbinned", markerfacecolor='#054950', zorder=2)

    minf_full = np.nanmin(np.array(flux_dd))
    maxf_full = np.nanmax(np.array(flux_dd))

    height_full = maxf - minf

    plt.vlines(line_dd, minf_full,minf_full + height_full*0.3, colors = md_color, label = "Momentum Dump", zorder=1)
    plt.ylim(minf_full, maxf_full)
    plt.xlim(np.nanmin(np.array(time_dd)), np.nanmax(np.array(time_dd)))
    plt.xlabel('Time (BJD-2457000)')
    plt.ylabel('Normalized Flux')

    plt.vlines(transit_list, minf_full,minf_full + height*0.3 , colors = 'darkred', linestyle = '--', linewidth = 4, zorder=2)

    plt.savefig("{}/{}/{}_fullLC_md.png".format(indir,tic,tic), format='png', bbox_inches = 'tight')

    plt.close()

def plot_full_md(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, transit_list, upper_axis_lim_final, lower_axis_lim_final, args):
    '''
    Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)).

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    allflux_err  :  list
        normalized flux errors (not binned)
    all_md  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    transit_list  :  list
        list of all the marked transit-events

    Returns
    -------
        Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)).
        The plot is always saved as this funtion is only called if the 'save' option was chosen
    '''
    gs = len(transit_list)

    plt.figure(figsize=(15,10))
    plt.tight_layout()

    if args.FFI == False:
        # rename things so that the code didn't have to be changed - not a very 'tidy' solution.
        time_dd = alltime
        flux_dd = allflux
        time_dd_binned = alltimebinned
        flux_dd_binned = allfluxbinned

    else:
        # --- do some sigma clipping to make the LC look better ----
        MAD = median_absolute_deviation(allflux,ignore_nan = True)
        madrange = (5 * MAD * 100)
        ymask = (allflux < 1 + madrange) * (allflux > 1 - madrange)

        time_dd = np.array(alltime)[ymask]
        flux_dd = np.array(allflux)[ymask]
        time_dd_binned = np.array(alltime)[ymask]
        flux_dd_binned = np.array(allflux)[ymask]
        # ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---

    line_dd = all_md

    minf_list = []
    maxf_list = []

    for peak in transit_list:

        mask_dd = (np.array(time_dd) < peak+1) & (np.array(time_dd) > peak-1)

        minf0 = np.nanmin(np.array(flux_dd)[mask_dd])
        maxf0 = np.nanmax(np.array(flux_dd)[mask_dd])

        minf_list.append(minf0)
        maxf_list.append(maxf0)


    minf = lower_axis_lim_final #np.nanmin(minf_list)
    maxf = upper_axis_lim_final #np.nanmin(maxf_list)
    height = maxf - minf

    for g,peak in enumerate(transit_list):

        mask_dd = (np.array(time_dd) < peak+1) & (np.array(time_dd) > peak-1)
        mask_dd_binned = (np.array(time_dd_binned) < peak+1) & (np.array(time_dd_binned) > peak-1)

        if np.sum(mask_dd) != 0:

            if gs == 1:
                plt.subplot(2,4,(6,7))

            else:
                plt.subplot(2,gs,(gs+g+1))

            height_cut = upper_axis_lim_final - lower_axis_lim_final #maxf - minf

            # plot the lines first so that they are behind the data - we need to be able to see the data well...
            plt.vlines(line_dd, minf-(height_cut/10), minf + height_cut*0.25 , colors = md_color, label = "Momentum Dump", zorder=1)
            plt.vlines([peak], minf-(height_cut/10),minf + height*0.25 , linewidth = 4, colors = 'darkred', linestyle = '--', zorder=2, alpha = 0.85)

            # plot the momentum dumps and the markings.
            plt.plot(np.array(time_dd)[mask_dd], np.array(flux_dd)[mask_dd], 'o', markersize = 4, color = main_pot_color, alpha = 0.8, label = "unbinned", markerfacecolor='white', zorder=3)
            plt.plot(np.array(time_dd_binned)[mask_dd_binned], np.array(flux_dd_binned)[mask_dd_binned], marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 5, label = 'binning = 7', markerfacecolor='k', zorder=4)

            plt.xlim(peak-1, peak+1)
            plt.ylim(upper_axis_lim_final, lower_axis_lim_final)

            try:
                plt.ylim(minf-height_cut/10, maxf + height_cut/10)
            except:
                print ('axis limits error (momentun dumps)')

            plt.title('Transit {}'.format(g+1))

        else:
            print ("size 0 so skip")


    if gs == 1:
        plt.subplot(2,4,(1,4))
    else:
        plt.subplot(2,gs,(1,gs))

    if args.FFI == False:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = main_pot_color, alpha = 0.9, label = "unbinned", markerfacecolor='white', zorder=2)
        plt.plot(np.array(time_dd_binned), np.array(flux_dd_binned), marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 1, label = 'binning = 7', markerfacecolor='k', zorder=3)
    else:
        plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = '#054950', alpha = 0.9, label = "unbinned", markerfacecolor='#054950', zorder=2)


    height_full = upper_axis_lim_final - lower_axis_lim_final

    plt.vlines(line_dd, lower_axis_lim_final,lower_axis_lim_final + height_full*0.3, colors = md_color, label = "Momentum Dump", zorder=1)
    plt.ylim(lower_axis_lim_final, upper_axis_lim_final)
    plt.xlim(np.nanmin(np.array(time_dd)), np.nanmax(np.array(time_dd)))
    plt.xlabel('Time (BJD-2457000)')
    plt.ylabel('Normalized Flux')

    plt.vlines(transit_list, lower_axis_lim_final, lower_axis_lim_final + height*0.3 , colors = 'darkred', linestyle = '--', linewidth = 4, zorder=2)

    plt.savefig("{}/{}/{}_fullLC_md.png".format(indir,tic,tic), format='png', bbox_inches = 'tight')

    plt.close()

# full light curve with the momentum dumps
def plot_full_md_notransits(tic, indir, alltime, allflux, all_md, alltimebinned, allfluxbinned, upper_axis_lim_final, lower_axis_lim_final, args):
    '''
    Plot of the full LC with the times of the momentum dumps marked (red lines). Single panel and no transit events - this is for the asteroseismology version of LATTE.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    allflux_err  :  list
        normalized flux errors (not binned)
    all_md  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    transit_list  :  list
        list of all the marked transit-events

    Returns
    -------
        Plot of the full LC with the times of the momentum dumps marked (red lines) and the marked transit events indicated (dashed balck line(s)).
        The plot is always saved as this funtion is only called if the 'save' option was chosen
    '''


    plt.figure(figsize=(15,6))

    line_dd = all_md
    time_dd = alltime
    flux_dd = allflux
    time_dd_binned = alltimebinned
    flux_dd_binned = allfluxbinned

    plt.plot(np.array(time_dd), np.array(flux_dd), 'o', markersize = 2, color = main_pot_color, alpha = 0.9, label = "unbinned", markerfacecolor='white', zorder=2)
    plt.plot(np.array(time_dd_binned), np.array(flux_dd_binned), marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 1, label = 'binning = 7', markerfacecolor='k', zorder=3)

    minf_full = np.nanmin(allflux)
    maxf_full = np.nanmax(allflux)

    height_full = maxf_full - minf_full

    plt.vlines(line_dd, minf_full,minf_full + height_full*0.3, colors = md_color, label = "Momentum Dump", zorder=1)
    plt.ylim(minf_full, maxf_full)
    plt.xlim(np.nanmin(np.array(time_dd)), np.nanmax(np.array(time_dd)))
    plt.xlabel('Time (BJD-2457000)')
    plt.ylabel('Normalized Flux')

    plt.savefig("{}/{}/{}_fullLC_md_notransits.png".format(indir,tic,tic), format='png', bbox_inches = 'tight')

    plt.close()

# plot of the two BLS runs including phase folded on most likely period
def plot_bls(tic, indir, alltime, allflux, alltimebinned, allfluxbinned, model, results, period, duration, t0, mid_transit_t0, args, in_transit = [0], in_transit_notbinned = [0]):
    '''
    Plot the BLS. This functinon is called in data_bls().

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    model :  float
        the transit model at the given period, duration, and phase
    results :  class
        results from the BLS fitting routine
    period :  float
        the period of the 'most-likely' transits
    duration :  float
        the duration of the transit
    t0  :  float
        the mid-transit time of the reference transit
    in_transit = [0] :  float
        if this is [0] (by deafult), the code knows that this is the initial run i.e. no transits have been removes (+ results are plotted in different colors)
    in_transit_notbinned = [0]. :  float
        if this is [0] (by deafult), the code knows that this is the initial run i.e. no transits have been removes (+ results are plotted in different colors)

    Returns
    -------
        Plot the results from the BLS with three pannels: periodgram, best fit model to the transits, phase folded fit.
    '''

    if len(in_transit) == 1:  # conditions for the first 'round' of plotting
        # define the colours of the plot
        color1 = '#DC143C'
        color2 = main_pot_color
        title = 'Initial BLS'
        name = '{}_bls_first.png'.format(tic) # name to be used to save the output file

    else:  # conditions for the second 'round' of plotting once the first event has been removed
        # define the colours of the plot
        color1 = 'deepskyblue'
        color2 = '#4682B4'
        title = 'Initial event removed'
        name = '{}_bls_second.png'.format(tic) # name to be used to save the output file

    fig, axes = plt.subplots(3, 1, figsize=(5, 7))

    # highlight the harmonics of the peak period
    ax = axes[0]
    ax.axvline(period, alpha=0.4, lw=5, color = color1)
    for n in range(2, 15):
        ax.axvline(n*period, alpha=0.4, lw=2, linestyle="dashed", color = color2) # plot the harmonics
        ax.axvline(period / n, alpha=0.4, lw=2, linestyle="dashed", color = color2)

    # ------------
    # plot the periodogram
    ax.plot(results.period, results.power, "k", lw=0.5, label = 'P = %.3f T0 = %.3f' % (period,mid_transit_t0))
    ax.set_title(title)
    ax.set_xlim(results.period.min(), results.period.max())
    ax.set_xlabel("period (days)")
    ax.set_ylabel("log likelihood")
    ax.legend(fontsize = 10, loc = 1)

    # ------------
    # plot the light curve and best-fit model
    ax = axes[1]

    if len(in_transit) == 1:  # for the initial run
        ax.plot(alltime, allflux, marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(alltimebinned, allfluxbinned, marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')
    else:  # for the second run (once the first 'event' has been removed)
        ax.plot(alltime[~in_transit_notbinned], allflux[~in_transit_notbinned], marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(alltimebinned[~in_transit], allfluxbinned[~in_transit], marker ="o", alpha = 0.6, color = 'black',  MarkerFaceColor = 'none', ms=3, lw = 0)

    x = np.linspace(alltimebinned.min(), alltimebinned.max(), 3*len(alltimebinned))
    f = model.model(x, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(alltimebinned.min(), alltimebinned.max())
    ax.set_xlabel("time (days)")
    ax.set_ylabel("de-trended flux (ppt)");

    # ------------
    ax = axes[2]
    if len(in_transit) == 1:  # for the initial run
        x_binned = (alltimebinned - t0 + 0.5*period) % period - 0.5*period
        x = (alltime - t0 + 0.5*period) % period - 0.5*period
    else: # for the second run (once the first 'event' has been removed)
        x_binned = (alltimebinned[~in_transit] - t0 + 0.5*period) % period - 0.5*period
        x = (alltime[~in_transit_notbinned] - t0 + 0.5*period) % period - 0.5*period

    m_binned = np.abs(x_binned) < 0.5
    m = np.abs(x) < 0.5

    # plot the data
    if len(in_transit) == 1:  # for the initial run
        ax.plot(x[m], allflux[m],marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(x_binned[m_binned], allfluxbinned[m_binned], marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')

    else: # for the second run (once the first 'event' has been removed)
        ax.plot(x[m], allflux[~in_transit_notbinned][m],marker =".", alpha = 0.4, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')
        ax.plot(x_binned[m_binned], allfluxbinned[~in_transit][m_binned], marker ="o", alpha = 0.6, color = 'black', ms=3, lw = 0, MarkerFaceColor = 'none')


    x = np.linspace(-0.5, 0.5, 1000)
    f = model.model(x + t0, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("time since transit (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    plt.tight_layout()

    # save the figures
    if args.save == True:
        plt.savefig('{}/{}/{}'.format(indir, tic, name), format='png')

    if args.noshow == False:
        plt.show()
    else:
        plt.close()

# plot of the two BLS runs including phase folded on most likely period for the FFIs
def plot_bls_FFI(tic, indir, alltime, allflux, model, results, period, duration, t0, mid_transit_t0, args, in_transit = [0]):
    '''

    Plot the BLS for the FFIs (these don't have binned data). This functinon is called in data_bls().

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    model :  float
        the transit model at the given period, duration, and phase
    results :  class
        results from the BLS fitting routine
    period :  float
        the period of the 'most-likely' transits
    duration :  float
        the duration of the transit
    t0  :  float
        the mid-transit time of the reference transit
    in_transit = [0] :  float
        if this is [0] (by deafult), the code knows that this is the initial run i.e. no transits have been removes (+ results are plotted in different colors)
    in_transit_notbinned = [0]. :  float
        if this is [0] (by deafult), the code knows that this is the initial run i.e. no transits have been removes (+ results are plotted in different colors)

    Returns
    -------
        Plot the results from the BLS with three pannels: periodgram, best fit model to the transits, phase folded fit.
    '''

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
    ax.plot(results.period, results.power, "k", lw=0.5, label = 'P = %.3f T0 = %.3f' % (period,mid_transit_t0))

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

    # -----------
    # phase fold
    x = (alltime - t0 + 0.5*period) % period - 0.5*period

    m = np.abs(x) < 0.5

    if len(in_transit) == 1:
        ax.plot(x[m], allflux[m],marker = ".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')

    else:
        ax.plot(x[m], allflux[m],marker = ".", alpha = 0.9, color = color2, ms=2, lw = 0, MarkerFaceColor = 'none')


    x = np.linspace(-0.5, 0.5, 1000)
    f = model.model(x + t0, period, duration, t0)
    ax.plot(x, f, lw=2, color = color1) # the line
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("time since transit (days)")
    ax.set_ylabel("de-trended flux (ppt)");
    plt.tight_layout()

    if args.save == True:
        plt.savefig('{}/{}/{}'.format(indir, tic, name), format='png')

    if args.noshow == False:
        plt.show()
    else:
        plt.close()

# plot of the in and out of transit average flux and the difference between them (North NOT up)
def plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, args):

    '''
    Plot the in transit average flux and the out of transit average flux and compare the two (difference image).
    The images are plotted in terms of pixels and are not oriented North.
    To have them oriented North run plot_in_out_TPF_proj (run code with --north in the command line) <-- this one takes longer.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    X4_list  :  list
        PCA corrected flux vs time for each pixel
    oot_list  :  list
        out of transit mask
    t_list  :  list
        time arrays
    intr_list  :  list
        in transit mask
    T0_list  :  list
        list of the peaks
    tpf_filt_list   : list
        list of the filtered tpfs

    Returns
    -------
    plot of the in average in transit flux, average out of transit flux and difference between the tw (out of transit ninu in transit) for each marked transit event.

    '''


    plt.figure(figsize=(16,3.5*len(T0_list)))

    plt.tight_layout()

    count = 0 # keep track of how many images have been plotted to that they appear on a subgrid of plots which has three columns

    # loop through all of the list of PCA corrected flux vs time arrays for each marked transit-event
    for idx, X4 in enumerate(X4_list): # idx is for each maked transit-event

        oot = oot_list[idx] # the out of transit flux mask
        intr = intr_list[idx] # the in-transit flux mask
        T0 = T0_list[idx] # the time of the transit-like event
        t = t_list[idx] # the time array
        tpf_filt  =  tpf_filt_list[idx]  # the filtered target pixel files

        intr = abs(T0-t) < 0.25  # create a mask of the in transit times
        oot = (abs(T0-t) < 0.5) * (abs(T0-t) < 0.3)  # create a mask of the out of transit times
        img_intr = tpf_filt[intr,:,:].sum(axis=0)/float(intr.sum()) # apply the masks and normalize the flux
        img_oot = tpf_filt[oot,:,:].sum(axis=0)/float(oot.sum())
        img_diff = img_oot-img_intr # calculate the diffefence image (out of transit minus in-transit)

        # ---- PLOT -------

        # in transit
        count += 1 # add to the count before each plot
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_intr, cmap=plt.cm.viridis, origin = 'lower')
        plt.colorbar()
        plt.title("t = {} days \n In Transit Flux (e-/candence)".format(T0))

        # out of transit
        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_oot, cmap=plt.cm.viridis, origin = 'lower')
        plt.colorbar()
        plt.title("Out of Transit Flux (e-/candence)")

        # out of transit minus in-transit
        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis('off')
        plt.imshow(img_diff, cmap=plt.cm.viridis, origin = 'lower')
        plt.colorbar()
        plt.title("Difference Flux (e-/candence)")


    plt.subplots_adjust(wspace = 0)

    # save the figure
    if args.save == True:
        plt.savefig('{}/{}/{}_flux_comparison.png'.format(indir, tic, tic), format='png', bbox_inches='tight')

    if args.noshow == False:
        plt.show()
    else:
        plt.close()

# plot of the in and out of transit average flux and the difference between them. Reprojected to have north upward facing.
def plot_in_out_TPF_proj(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, tpf_list, args):
    '''
    Plot the in transit average flux and the out of transit average flux and compare the two (difference image).
    The images are oriented so that North is up. Note that this takes longer to run than if they are not re-projected.

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    X4_list  :  list
        PCA corrected flux vs time for each pixel
    oot_list  :  list
        out of transit mask
    t_list  :  list
        time arrays
    intr_list  :  list
        in transit mask
    T0_list  :  list
        list of the peaks
    tpf_filt_list   : list
        list of the filtered tpfs
    tpf_list   : list
        list of the target pixel files (for each sector)

    Returns
    -------
    plot of the in average in transit flux, average out of transit flux and difference between the tw (out of transit ninu in transit) for each marked transit event.
    All the plots are oriented so that north is upwards.

    '''

    plt.figure(figsize=(15,3.5*len(T0_list)))

    count = 0 # keep track of how many images have been plotted to that they appear on a subgrid of plots which has three columns

    # loop through all of the list of PCA corrected flux vs time arrays for each marked transit-event
    for idx, X4 in enumerate(X4_list):

        oot = oot_list[idx] # the out of transit flux mask
        intr = intr_list[idx] # the in-transit flux mask
        T0 = T0_list[idx] # the time of the transit-like event
        t = t_list[idx] # the time array
        tpf_filt  =  tpf_filt_list[idx] # the filtered target pixel files
        tpf = tpf_list[idx]

        if (args.FFI == False) or (args.FFI == 'SPOC'):
            tpf_wcs = extract_wcs(tpf)
        else:
            tpf_wcs = tpf.wcs
        # create a tupple of the array of the data and the wcs projection of the TESS cutout
        tup = (np.nanmean(tpf[1].data['FLUX'],axis=0),tpf_wcs)

        # map the output so that the image will be oriented NORTH
        # the find_optimal_celestial_wcs function returns new world coordinate system (wcs) orented north that can be used to map the images
        wcs_out, shape_out = find_optimal_celestial_wcs(input_data =[tup], resolution = 0.0002*u.deg)

        # re-project the image with the new wcs
        array, footprint = reproject_interp(tup, wcs_out, shape_out = shape_out, order = 'nearest-neighbor')

        intr = abs(T0-t) < 0.25 # create a mask of the in transit times
        oot = (abs(T0-t) < 0.5) * (abs(T0-t) < 0.3) # the in-transit flux mask
        img_intr = tpf_filt[intr,:,:].sum(axis=0)/float(intr.sum())  # array of the time of the transit-like event
        img_oot = tpf_filt[oot,:,:].sum(axis=0)/float(oot.sum())     # array
        img_diff = img_oot-img_intr                                  # array of the diffefence image (out of transit minus in-transit)

        # ---- RE-PROJECT ----
        # create tupples of the image array and the wcs projection of the TESS cutout
        tup_intr = (img_intr, tpf_wcs)
        tup_oot  = (img_oot, tpf_wcs)
        tup_diff = (img_diff, tpf_wcs)

        # re-project the image with the new (north oriented) wcs
        img_intr, _ = reproject_interp(tup_intr, wcs_out, shape_out = shape_out, order = 'nearest-neighbor')
        img_oot, _ = reproject_interp(tup_oot, wcs_out, shape_out = shape_out, order = 'nearest-neighbor')
        img_diff, _ = reproject_interp(tup_diff, wcs_out, shape_out = shape_out, order = 'nearest-neighbor')

        # --- PLOT ---

        # in transit
        count += 1 # add to the count before each plot
        plt.subplot(len(T0_list), 3, count,projection=wcs_out)
        plt.axis('off')
        plt.imshow(img_intr)
        plt.colorbar()
        plt.xlabel("RA", fontsize = main_plot_fontsize)
        plt.ylabel("Dec", fontsize = main_plot_fontsize)
        plt.title("t = {} days \n In Transit Flux (e-/candence)".format(T0), fontsize = main_plot_fontsize)

        # out of transit
        count += 1
        plt.subplot(len(T0_list), 3, count,projection=wcs_out)
        plt.axis('off')
        plt.imshow(img_oot)
        plt.colorbar()
        plt.xlabel("RA", fontsize = main_plot_fontsize)
        plt.title("Out of Transit Flux (e-/candence)", fontsize = main_plot_fontsize)

        # out of transit minus in transit
        count += 1
        plt.subplot(len(T0_list), 3, count,projection=wcs_out)
        plt.axis('off')
        plt.imshow(img_diff)
        plt.colorbar()
        plt.xlabel("RA", fontsize = main_plot_fontsize)
        plt.title("Difference Flux (e-/candence)", fontsize = main_plot_fontsize)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.45)
    # save the figure
    if args.save == True:
        plt.savefig('{}/{}/{}_flux_comparison.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5)

    if args.noshow == False:
        plt.show()
    else:
        plt.close()


# make a periodogram. This uses lighkurve.
def plot_periodogram(tic, indir, alltime, allflux, teff, input_numax, input_analysis_window, args, ql = False):
    '''
    Plot a periodigram. This function makes use of the lightkurve module

    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    alltime  :  list

    allfux   :  list

    Returns
    -------
    plot of the periodogram.

    '''
    alltime = np.array(alltime)
    allflux = np.array(allflux)
    # make sure that the data doesn't contain any nans
    finite_mask = np.isfinite(alltime) & np.isfinite(allflux)

    #use the lightkurve function
    lc = lk.lightcurve.LightCurve(time = alltime[finite_mask], flux = allflux[finite_mask])

    #input_numax, input_analysis_window

    lc = lc.remove_outliers().remove_nans()
    ls = lc.to_periodogram(method='lombscargle', normalization='psd')

    smooth = ls.smooth(method='boxkernel', filter_width=20.)

    freq = ls.frequency
    power = ls.power

    freq_smooth = smooth.frequency
    power_smooth = smooth.power

    # - - - plot - - -
    plt.figure(figsize=(7,4))

    plt.plot(freq, power, color = 'k', lw = 0.5)
    plt.plot(freq_smooth, power_smooth, color = '#bf4080', lw = 2)

    # set the limit - with TESS data we won't see anything lower than 10 Hz anyway so cute ter=
    plt.xlim(10,np.nanmax(freq.value))

    plt.xlabel(r"$Frequency~(\mu Hz)$", fontsize = main_plot_fontsize)
    plt.ylabel(r"$Power~Spectral~Density~(A^{2} \mu Hz^{-1})$", fontsize = main_plot_fontsize)

    plt.loglog() # make the axis log log
    plt.tight_layout()

    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_periodogram.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5, dpi = 100)

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()


    if ql == False:
        # if the asteroseismic analysis was chosen

        if input_numax != -999:

            # make another periodogram with the filter width already determined! this is needed for the asteroseismic analysis to work

            pg = lc.to_periodogram(method='lombscargle', normalization='psd', minimum_frequency=input_numax-input_analysis_window, maximum_frequency=input_numax+input_analysis_window)
            snr = pg.flatten()
            seis = snr.to_seismology()

            # determine delta nu and numax - these are needed for the diagnostic plots
            # these are LK objects that have lots of information in them - like magic...

            numax = seis.estimate_numax()
            deltanu = seis.estimate_deltanu()

            # call on the two plotting functions that are adaptations of the LK plotting functions (this is LK code, not my own code!)
            plot_diagnostic_numax(tic, indir, numax, deltanu, pg, args)
            plot_scaled_correlation(tic, indir, deltanu, pg, args)
            plot_echelle(tic, indir, seis, args)

            # extract the parameters that we want - need to give it Teff
            mass_ast = seis.estimate_mass(teff)
            radius_ast = seis.estimate_radius(teff)
            logg_ast = seis.estimate_logg(teff)

            return mass_ast.value, radius_ast.value, logg_ast.value, numax.value, deltanu.value

        else:
            # if the asteroseismic mode was not selected
            return -999,-999,-999,-999,-999


def get_fwhm(periodogram, numax):

    # Function from LK to get the fwhm (NOT MY OWN CODE)

    if u.Quantity(periodogram.frequency[-1], u.microhertz) > u.Quantity(500., u.microhertz):
        fwhm = 0.25 * numax
    else:
        fwhm = 0.66 * numax**0.88
    return fwhm


def plot_diagnostic_numax(tic, indir, numax, deltanu, periodogram, args, ql = False):

    # adaption of function from LK (NOT MY OWN CODE)

    if ql == True:
        fig, ax = plt.subplots(3, sharex=True, figsize=(6, 7))
    else:
        fig, ax = plt.subplots(3, sharex=True, figsize=(8, 12))

    periodogram.plot(ax=ax[0], label='', color = 'k')
    ax[0].axvline(deltanu.diagnostics['numax'].value, c=main_pot_color, linewidth=2,
               alpha=1, ls='-')
    ax[0].text(deltanu.diagnostics['numax'].value, periodogram.power.value.max()*0.45,
            '{} ({:.1f} {})'.format(r'$\nu_{\rm max}$', deltanu.diagnostics['numax'].value, deltanu.diagnostics['numax'].unit.to_string('latex')),
            rotation=90, ha='right', color='darkorange', alpha=1, fontsize=10)
    ax[0].text(.025, .9, 'Input Power Spectrum', horizontalalignment='left',
            transform=ax[0].transAxes, fontsize=main_plot_fontsize)


    window_width = 2*int(np.floor(get_fwhm(periodogram, deltanu.diagnostics['numax'].value)))
    frequency_spacing = np.median(np.diff(periodogram.frequency.value))
    spread = int(window_width/2/frequency_spacing)  # spread in indices

    a = np.argmin(np.abs(periodogram.frequency.value - deltanu.diagnostics['numax'].value)) + spread
    b = np.argmin(np.abs(periodogram.frequency.value - deltanu.diagnostics['numax'].value)) - spread

    a = [periodogram.frequency.value[a] if a < len(periodogram.frequency)
         else periodogram.frequency.value[-1]][0]
    b = [periodogram.frequency.value[b] if b > 0
         else periodogram.frequency.value[0]][0]

    #ax[0].axvline(a, c='#bf4080', linewidth=2, alpha=1, ls='--')
    #ax[0].axvline(b, c='#bf4080', linewidth=2, alpha=1, ls='--')

    h = periodogram.power.value.max() * 0.9

    #ax[0].annotate("", xy=(a, h), xytext=(a - (a-b), h),
    #                arrowprops=dict(arrowstyle="<->", color='#bf4080', alpha=1), va='bottom')
    #ax[0].text(a - (a-b)/2, h, r"2 $\times$ FWHM'", color='#bf4080', alpha=1, fontsize=10, va='bottom', ha='center')

    ax[0].set_xlim(b - ((a-b)*0.2), a + ((a-b)*0.2))


    vmin = np.nanpercentile(numax.diagnostics['acf2d'], 5)
    vmax = np.nanpercentile(numax.diagnostics['acf2d'], 95)


    ax[1].pcolormesh(numax.diagnostics['numaxs'],
                     np.linspace(0, numax.diagnostics['window_width'],
                     num=numax.diagnostics['acf2d'].shape[0]),
                     numax.diagnostics['acf2d'],
                     cmap='Greys',
                     vmin=vmin,
                     vmax=vmax)

    ax[1].set_ylabel(r'Frequency lag [{}]'.format(periodogram.frequency.unit.to_string('latex')))
    ax[1].axvline(numax.value, c=main_pot_color, linewidth=2, alpha=1)
    ax[1].text(.025, .9, '2D AutoCorrelation',
                horizontalalignment='left',
                transform=ax[1].transAxes, fontsize=main_plot_fontsize)

    ax[2].plot(numax.diagnostics['numaxs'], numax.diagnostics['metric'], color = 'k')
    ax[2].plot(numax.diagnostics['numaxs'], numax.diagnostics['metric_smooth'],
               lw=2, alpha=1, label='Smoothed Metric', c = '#bf4080')
    ax[2].set_xlabel("Frequency [{}]".format(periodogram.frequency.unit.to_string('latex')))
    ax[2].set_ylabel(r'Correlation Metric')

    ax[2].axvline(numax.value, c=main_pot_color, linewidth=2, alpha=1)
    ax[2].text(.025, .9, 'Correlation Metric',
                        horizontalalignment='left',
                        transform=ax[2].transAxes, fontsize=main_plot_fontsize)
    ax[2].legend(loc='upper right')
    ax[2].set_xlim(numax.diagnostics['numaxs'][0], numax.diagnostics['numaxs'][-1])

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.tight_layout()

    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_diagnostic_numax.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5, dpi = 100)

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()



def plot_scaled_correlation(tic, indir, deltanu, periodogram, args, ql = False):

    # adaption of function from LK (NOT MY OWN CODE)

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(deltanu.diagnostics['lags'][2:], deltanu.diagnostics['acf'][2:], color = 'k',  linewidth=1)
    ax.set_xlabel("Frequency Lag [{}]".format(periodogram.frequency.unit.to_string('latex')))
    ax.set_ylabel(r'Scaled Auto Correlation', fontsize=main_plot_fontsize)

    axin = inset_axes(ax, width="50%",height="50%", loc="upper right")
    axin.set_yticks([])
    axin.plot(deltanu.diagnostics['lags'][deltanu.diagnostics['sel']],
              deltanu.diagnostics['acf'][deltanu.diagnostics['sel']], color = 'k', linewidth=1)
    axin.scatter(deltanu.diagnostics['lags'][deltanu.diagnostics['sel']][deltanu.diagnostics['peaks']],
                 deltanu.diagnostics['acf'][deltanu.diagnostics['sel']][deltanu.diagnostics['peaks']],
                 c='r', s=5)

    mea_label = r'Measured {} {:.1f} {}'.format(
                    r'$\Delta\nu$',
                    deltanu.value,
                    periodogram.frequency.unit.to_string('latex'))
    ax.axvline(deltanu.value,c='r', linewidth=1.5, alpha=1, label=mea_label, zorder = -10)

    emp_label = r'Empirical {} {:.1f} {}'.format(
                    r'$\Delta\nu$',
                    deltanu.diagnostics['deltanu_emp'],
                    periodogram.frequency.unit.to_string('latex'))
    ax.axvline(deltanu.diagnostics['deltanu_emp'], c=md_color, linewidth=1.5,
               alpha=1, ls='--', label=emp_label, zorder = -10)

    axin.axvline(deltanu.value,c='r', linewidth=1.5, alpha=1, zorder = -10)
    axin.axvline(deltanu.diagnostics['deltanu_emp'], c=md_color, linewidth=1.5,
                 alpha=1, ls='--', zorder = -10)

    ax.legend(loc='lower right', fontsize=main_plot_fontsize)
    plt.tight_layout()

    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_scaled_correlation.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5, dpi = 150)

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()


def plot_echelle(tic, indir, seis, args, smooth_filter_width=.3, cmap='Greys', ql = False):

    """Plots an echelle diagram of the periodogram by stacking the
    periodogram in slices of deltanu. THIS IS FUNCTION IS FROM LIGHTKURVE - NOT MY OWN CODE!!!"""

    numax = seis._validate_numax(None)
    deltanu = seis._validate_deltanu(None)

    if smooth_filter_width:
        pgsmooth = seis.periodogram.smooth(filter_width=smooth_filter_width)
        freq = pgsmooth.frequency  # Makes code below more readable below
        power = pgsmooth.power     # Makes code below more readable below
    else:
        freq = seis.periodogram.frequency  # Makes code below more readable
        power = seis.periodogram.power     # Makes code below more readable

    fmin = freq[0]
    fmax = freq[-1]

    # Ensure input numax is in the correct units (if there is one)
    if numax is not None:
        numax = u.Quantity(numax, freq.unit).value

        fwhm = get_fwhm(seis.periodogram, numax)

        fmin = numax - 2*fwhm
        if fmin < freq[0].value:
            fmin = freq[0].value

        fmax = numax + 2*fwhm
        if fmax > freq[-1].value:
            fmax = freq[-1].value


    # Make sure fmin and fmax are Quantities or code below will break
    fmin = u.Quantity(fmin, freq.unit)
    fmax = u.Quantity(fmax, freq.unit)


    # Add on 1x deltanu so we don't miss off any important range due to rounding
    if fmax < freq[-1] - 1.5*deltanu:
        fmax += deltanu

    fs = np.median(np.diff(freq))
    x0 = int(freq[0] / fs)

    ff = freq[int(fmin/fs)-x0:int(fmax/fs)-x0] # Selected frequency range
    pp = power[int(fmin/fs)-x0:int(fmax/fs)-x0] # Power range

    n_rows = int((ff[-1]-ff[0])/deltanu) # Number of stacks to use
    n_columns = int(deltanu/fs)          # Number of elements in each stack

    # Reshape the power into n_rows of n_columns
    ep = np.reshape(pp[:(n_rows*n_columns)], (n_rows, n_columns))

    #Reshape the freq into n_rowss of n_columnss & create arays
    ef = np.reshape(ff[:(n_rows*n_columns)],(n_rows,n_columns))
    x_f = ((ef[0,:]-ef[0,0]) % deltanu)
    y_f = (ef[:,0])

    #Plot the echelle diagram
    fig, ax = plt.subplots(figsize = (7,4))
    extent = (x_f[0].value, x_f[-1].value, y_f[0].value, y_f[-1].value)
    figsize = [7.,4.]
    a = figsize[1] / figsize[0]
    b = (extent[3] - extent[2]) / (extent[1] - extent[0])
    vmin = np.nanpercentile(ep.value, 1)
    vmax = np.nanpercentile(ep.value, 99)
    im = ax.imshow(ep.value, cmap=cmap, aspect=a/b, origin='lower',
              extent=extent, vmin=vmin, vmax=vmax)


    ax.set_xlabel(r'Frequency mod. {:.2f}'.format(deltanu), fontsize = main_plot_fontsize)
    ax.set_ylabel(r'Frequency ({})'.format(freq.unit.to_string('latex')), fontsize = main_plot_fontsize)
    plt.tight_layout()

    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_echelle.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5, dpi = 100)

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()


#plot_echelle(tic, indir, seis, args)


# plot the equivilaten evolutionary phase tracks
def eep_target(tic, indir, syspath, temp, rad, args, ql = False):

    # the eep tracks are downloaded when you install the code. These files are not large ~70 KB each
    phase0 = pd.read_csv('{}/LATTE_eep_data/eep_phase0.csv'.format(syspath)) # these are the main-sequence tracks
    phase2 = pd.read_csv('{}/LATTE_eep_data/eep_phase2.csv'.format(syspath)) # these are the post main-sequence tracks

    # we need the TOI file which has the stellar parameters (unfortunately this is a different file to the other TOI file - both are installed with --new-data)
    infile = "{}/data/TOI_list_star_params.txt".format(indir)
    exoplanets = pd.read_csv(infile, comment = '#')

    # if either the radius or the temperature is unknown, ust ignore it because it can't be plotted.
    # also uknown radii default to 1, so get rid of these as they are wrong.
    exoplanets_r = exoplanets[~(np.isnan(exoplanets['Stellar Radius (R_Sun) err']) & (exoplanets['Stellar Radius (R_Sun)'] == 1))]

    # set up the plot
    fig, ax = plt.subplots(figsize=(7,4))

    for i in range(300,1700, 100): # we have evolutionary tracks from 0.3 M_sun to 1.6 M_sun in steps of 0.1
        colname_phase0_T = (str(i).rjust(6, "0")) + '0M_phase0_T'
        colname_phase0_R = (str(i).rjust(6, "0")) + '0M_phase0_R'

        colname_phase2_T = (str(i).rjust(6, "0")) + '0M_phase2_T'
        colname_phase2_R = (str(i).rjust(6, "0")) + '0M_phase2_R'

        phase0_T = phase0[colname_phase0_T]
        phase0_R = phase0[colname_phase0_R]

        phase2_T = phase2[colname_phase2_T]
        phase2_R = phase2[colname_phase2_R]

        # the -99 values are just there as padding to make the lists the same length -- we dont' need them to mask out.
        mask0 = (phase0_T != -99) & (phase0_R != -99)
        mask2 = (phase2_T != -99) & (phase2_R != -99)

        # plot the 1 Solar Mass track in a different colour to make it stand out - easier for reference.
        if i == 1000:
            plt.plot(phase0_T[mask0], phase0_R[mask0], 'maroon', zorder = -1, lw = 2)
            plt.plot(phase2_T[mask2],phase2_R[mask2], 'maroon', linestyle = ':', zorder = -1, lw = 2)
        else:
            plt.plot(phase0_T[mask0], phase0_R[mask0], 'k', zorder = -1)
            plt.plot(phase2_T[mask2],phase2_R[mask2], 'grey', linestyle = ':', zorder = -1)


    # plot the TOI's and the target in question.

    plt.scatter(temp, rad, c= '#e60063', s = 25, zorder = 1)
    plt.scatter(exoplanets_r['Stellar Eff Temp (K)'], exoplanets_r['Stellar Radius (R_Sun)'], zorder = -2, c = 'navy', alpha = 0.15, s = 4)

    plt.xlim(plt.xlim()[::-1]) # refert the axis.

    try:
        if rad < 3:
            plt.ylim(0.2,3.3)
        else:
            plt.ylim(0.2,5.3)
    except:
        plt.ylim(0.2,5.3)

    plt.xlabel(r"${\rm T_{eff}~(K)}$", fontsize = 15)
    plt.ylabel(r"${\rm Radius~(R_{\odot})}$", fontsize = 15)

    try:
    	if temp < 8000:
    	    plt.xlim(8000,3000)
    	else:
    	    plt.xlim(12000,3000)
    except:
    	plt.xlim(12000,3000)

    plt.subplots_adjust(hspace = 0, top=0.98, right = 0.965, bottom = 0.18, left = 0.1)

    if (args.save == True) and (ql == False):
        plt.savefig('{}/{}/{}_eep.png'.format(indir, tic, tic), format='png', bbox_inches='tight', pad_inches=0.5)

    if (args.noshow == False) or (ql == True):
        plt.show()
    else:
        plt.close()


# --------------------------------------------
#                other functions             #
# --------------------------------------------

def rebin(arr,new_shape):

    ''''
    function used to rebin the data
    '''
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
        new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

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


def find_aperture_FFI(val, ap_count, tpf):
    '''
    function that finds the desired aperture.
    Used with a scipy minimization to find the threshhold value that results in the numbre of pixels wanted for the aperture.
    '''

    flux_array = tpf.flux

    # calculate the mask centerd on the central pixel (on target)
    target_mask = extract_aperture(flux_array, threshold=val)

    # determine how many pixels this threshhold results in and compare to desired mask size
    return abs(np.sum(target_mask) - ap_count)


def find_aperture(val, ap_count, flux_array):
    '''
    function that finds the desired aperture.
    Used with a scipy minimization to find the threshhold value that results in the numbre of pixels wanted for the aperture.
    '''

    # calculate the mask centerd on the central pixel (on target)
    target_mask = extract_aperture(flux_array, threshold=val)

    # determine how many pixels this threshhold results in and compare to desired mask size
    return abs(np.sum(target_mask) - ap_count)

def extract_aperture(flux_array, threshold, reference_pixel='center'):
    """this is a lightkurve function!"""
    if reference_pixel == 'center':
        reference_pixel = (flux_array.shape[2] / 2, flux_array.shape[1] / 2)
    # Calculate the median image

    median_image = np.nanmedian(flux_array, axis=0)
    vals = median_image[np.isfinite(median_image)].flatten()
    # Calculate the theshold value in flux units
    mad_cut = (1.4826 * median_absolute_deviation(vals) * threshold) + np.nanmedian(median_image)
    # Create a mask containing the pixels above the threshold flux
    threshold_mask = np.nan_to_num(median_image) >= mad_cut
    if (reference_pixel is None) or (not threshold_mask.any()):
        # return all regions above threshold
        return threshold_mask
    else:
        # Return only the contiguous region closest to `region`.
        # First, label all the regions:
        labels = label(threshold_mask)[0]
        # For all pixels above threshold, compute distance to reference pixel:
        label_args = np.argwhere(labels > 0)
        distances = [np.hypot(crd[0], crd[1])
                     for crd in label_args - np.array([reference_pixel[1], reference_pixel[0]])]
        # Which label corresponds to the closest pixel?
        closest_arg = label_args[np.argmin(distances)]
        closest_label = labels[closest_arg[0], closest_arg[1]]
        return labels == closest_label


def extract_wcs(tpf):
    '''
    Extract the wcs from the header - function from lighkurve
    '''
    wcs_keywords = {'1CTYP5': 'CTYPE1',
                    '2CTYP5': 'CTYPE2',
                    '1CRPX5': 'CRPIX1',
                    '2CRPX5': 'CRPIX2',
                    '1CRVL5': 'CRVAL1',
                    '2CRVL5': 'CRVAL2',
                    '1CUNI5': 'CUNIT1',
                    '2CUNI5': 'CUNIT2',
                    '1CDLT5': 'CDELT1',
                    '2CDLT5': 'CDELT2',
                    '11PC5': 'PC1_1',
                    '12PC5': 'PC1_2',
                    '21PC5': 'PC2_1',
                    '22PC5': 'PC2_2',
                    'NAXIS1': 'NAXIS1',
                    'NAXIS2': 'NAXIS2'}
    mywcs = {}

    for oldkey, newkey in wcs_keywords.items():
        if (tpf[1].header[oldkey] != pf.Undefined):
            mywcs[newkey] = tpf[1].header[oldkey]
    return WCS(mywcs)


def mastQuery(request, proxy_uri=None):

    host = "mast.stsci.edu"
    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent": "python-requests/" + version,
    }

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    # opening the https connection
    if None == proxy_uri:
        conn = httplib.HTTPSConnection(host)
    else:
        port = 443
        url = urlparse(proxy_uri)
        conn = httplib.HTTPSConnection(url.hostname, url.port)

        if url.username and url.password:
            auth = "%s:%s" % (url.username, url.password)
            headers["Proxy-Authorization"] = "Basic " + str(
                base64.b64encode(auth.encode())
            ).replace("b'", "").replace("'", "")
        conn.set_tunnel(host, port, headers)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request=" + requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode("utf-8")

    # Close the https connection
    conn.close()

    return head, content


def interact_LATTE_test(tic, indir, syspath, sectors_all, sectors, ra, dec, args):

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

    # function needed in the rebinning - changes the shape of the array
    def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
    # ---------------
    global in_sec

    # call function to download the data from MAST
    print ("Start data download.....", end =" ")
    alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = download_data(indir, sectors, tic)
    print ("done.\n")

    plt.close('all') # make sure that all figures are close.

    # -------------------------
    # Plot the interactive plot - uses matplolib
    # -------------------------

    fig, ax = plt.subplots(2, 1, figsize=(10,7))
    plt.tight_layout()

    # Adjust the plots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.24, bottom=0.3)

    fluxmin = np.nanmin(allflux)
    fluxmax = np.nanmax(allflux)

    # function to define the plotting area around the transit event.
    # this needs to be in a function as the area changes with the interactive slider.
    def cutout(transit):
        mask_binned = (np.array(alltimebinned) > transit-1) & (np.array(alltimebinned) < transit+1)
        mask = (np.array(alltime) > transit-1) & (np.array(alltime) < transit+1)

        return [np.array(alltime)[mask], np.array(allflux)[mask], np.array(alltime), np.array(allflux), np.array(alltimebinned)[mask_binned], np.array(allfluxbinned)[mask_binned], np.array(alltimebinned), np.array(allfluxbinned)]

    # ---------------

    # function to define the bin factor
    # this needs to be in a function as binning can be changed with the interactive buttons.
    def binning(binfac):
        N      = len(alltime)
        n      = int(np.floor(N/binfac)*binfac)
        X      = np.zeros((2,n))
        X[0,:]  = alltime[:n]
        X[1,:]  = allflux[:n]
        Xb    = rebin(X, (2,int(n/binfac)))

        time_binned = Xb[0]
        flux_binned = Xb[1]

        return [time_binned, flux_binned]
    # ---------------


    class PageSlider(matplotlib.widgets.Slider):

        def __init__(self, ax, label, numpages = 10, valinit=0, valfmt='%1d',
                     closedmin=True, closedmax=True,
                     dragging=True, **kwargs):

            self.facecolor=kwargs.get('facecolor',"w")
            self.activecolor = kwargs.pop('activecolor',"b")
            self.fontsize = kwargs.pop('fontsize', 10)
            self.numpages = numpages

            super(PageSlider, self).__init__(ax, label, 0, numpages,
                                valinit=valinit, valfmt=valfmt, **kwargs)

            # make it so that numpages is the index of the array of the
            global in_sec

            self.poly.set_visible(False)
            self.vline.set_visible(False)
            self.pageRects = []
            for i in range(numpages):

                facecolor = self.activecolor if i==valinit else self.facecolor
                r  = Rectangle((float(i)/numpages, 0), 1./numpages, 1,
                                    transform=ax.transAxes, facecolor=facecolor)
                ax.add_artist(r)
                self.pageRects.append(r)

                seclist_names = ['All'] + in_sec

                sec_i = seclist_names[i]

                ax.text(float(i)/numpages+0.5/numpages, 0.5, str(sec_i),
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=self.fontsize)
            self.valtext.set_visible(False)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            bax = divider.append_axes("right", size="5%", pad=0.05)
            fax = divider.append_axes("right", size="5%", pad=0.05)
            self.button_back = matplotlib.widgets.Button(bax, label=r'$\blacktriangleleft$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_forward = matplotlib.widgets.Button(fax, label=r'$\blacktriangleright$',
                            color=self.facecolor, hovercolor=self.activecolor)
            self.button_back.label.set_fontsize(self.fontsize)
            self.button_forward.label.set_fontsize(self.fontsize)
            self.button_back.on_clicked(self.backward)
            self.button_forward.on_clicked(self.forward)

        def _update(self, event):
            super(PageSlider, self)._update(event)
            i = int(self.val)
            if i >=self.valmax:
                return
            self._colorize(i)

        def _colorize(self, i):
            for j in range(self.numpages):
                self.pageRects[j].set_facecolor(self.facecolor)
            self.pageRects[i].set_facecolor(self.activecolor)

        def forward(self, event):
            current_i = int(self.val)
            i = current_i+1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)

        def backward(self, event):
            current_i = int(self.val)
            i = current_i-1
            if (i < self.valmin) or (i >= self.valmax):
                return
            self.set_val(i)
            self._colorize(i)



    # define the initial conditions - this is what appears on the plot initially but can change later.
    transit = np.nanmean(alltimebinned)  # start the cut-out in the centre of the plot.
    binfac = 7

    # FIRST PLOT - FULL LC
    # plot the unbinned flux
    [line_full] = ax[0].plot(alltime, allflux , marker='o',lw = 0, markersize = 4, color = main_pot_color, alpha = 0.8,  markerfacecolor='white')
    # plot the binned flux - it is called through the above function so that the binning can be changed.
    [line_full_binned] = ax[0].plot(binning(binfac)[0], binning(binfac)[1],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3,  markerfacecolor='k')
    # ---------------

    # SECOND PLOT - CUT OUT LC
    # plot the cut out around the time of the transit - called with function as it can be changed with the slider.
    [line] =  ax[1].plot(cutout(transit)[0], cutout(transit)[1], marker='o',lw = 0, markersize = 4, color = main_pot_color, alpha = 0.8, markerfacecolor='white')
    [line_binned] =  ax[1].plot(cutout(transit)[4], cutout(transit)[5],marker='o',color = 'k', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='k')
    # ---------------

    global transit_slider_ax
    global transit_slider

    # Define the slider to change the transit-event time (and cut out region)
    transit_slider_ax  = fig.add_axes([0.25, 0.14, 0.65, 0.03])  # location of the slider
    transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=transit, color='teal')

    # Define the slider to change the y axis scale
    scale_slider_ax  = fig.add_axes([0.25, 0.19, 0.65, 0.03])
    scale_slider = Slider(scale_slider_ax, 'Y-Axis Scale', 0.99, 1.01, valinit=1, color='silver')

    # define the slider to change the plotting region
    sector_slider_ax  = fig.add_axes([0.25, 0.10, 0.65, 0.03])
    sector_slider = PageSlider(sector_slider_ax, 'Sector', len(in_sec), activecolor=main_pot_color)

    # ---------------

    # deifne the intial x and y axis limits - y limit can be changed with slider, x limit cannot.
    ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
    ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
    ax[0].set_ylim([fluxmin, fluxmax])

    ax[1].set_xlim([np.nanmean(alltime)-1, np.nanmean(alltime)+1])
    ax[1].set_ylim([fluxmin, fluxmax])

    # ---------------
    # Define an action for acessing the required cut-out data and drawing the line when the slider's value changes
    def sliders_on_changed(val):
        line.set_xdata(cutout(transit_slider.val)[0])
        line.set_ydata(cutout(transit_slider.val)[1])

        line_binned.set_xdata(cutout(transit_slider.val)[4])
        line_binned.set_ydata(cutout(transit_slider.val)[5])

        fig.canvas.draw_idle()  # update figure


    # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
    lver0 = ax[0].axvline(transit, color = 'r', linewidth = 2)
    lver1 = ax[1].axvline(transit, color = 'r', linewidth = 2)

    # ---------------
    # Define an action for modifying the plot region on the second plot
    def update_axis(val):
        # ---------------
        # only ever plot one line so get rid of the old ones before plotting the new ones...
        if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
            ax[0].lines[-1].remove()
            ax[1].lines[-1].remove()

        # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
        lver0 = ax[0].axvline(transit_slider.val, color = 'r', linewidth = 2)
        lver1 = ax[1].axvline(transit_slider.val, color = 'r', linewidth = 2)

        ax[1].set_xlim([transit_slider.val - 1,transit_slider.val + 1])
        lver0.set_xdata(transit_slider.val)
        lver1.set_xdata(transit_slider.val)


    def update_plotting_region_old(val):
        global seclist_names

        index = (seclist_names.index(val))

        if index == 0:
            ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
            ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
        else:
            ax[0].set_title("TIC {}  {}".format(tic,val))
            ax[0].set_xlim([start_sec[index-1][0],end_sec[index-1][0]])


    def update_plotting_region(val):
        global seclist_names

        val = int(val)
        index = val

        if index == 0:
            ax[0].set_title("TIC {}    Sectors {}".format(tic,str(in_sec)[1:-1]))
            ax[0].set_xlim([np.nanmin(alltime), np.nanmax(alltime)])
        else:
            ax[0].set_title("TIC {}  {}".format(tic,val))
            ax[0].set_xlim([start_sec[index-1][0],end_sec[index-1][0]])


    def update_yaxis(val):

        med = 1
        diff = abs(med - (fluxmin * scale_slider.val))

        ax[0].set_ylim([med - diff ,med + diff])
        ax[1].set_ylim([med - diff ,med + diff])


    # ------ ON CLICK -------
    # define the action that lets you click on the image in order to chose the location instead of using the slider - probably more useful
    def onclick(event):

        # just need to know the x (time) position of the click. The y (flux) position doesn't matter.
        val = event.xdata

        # check that the click was within the plotting region. If not ignore the click
        if (event.inaxes == ax[0]) or (event.inaxes == ax[1]):

            # SECOND PLOT - CUT OUT LC
            # plot the cut out around the time of the transit - changed by clicked on the image (either one of them)
            line.set_xdata(cutout(val)[0])
            line.set_ydata(cutout(val)[1])

            line_binned.set_xdata(cutout(val)[4])
            line_binned.set_ydata(cutout(val)[5])

            # ---------------
            # only ever plot one line so get rid of the old ones before plotting the new ones...
            if (len (ax[0].lines) > 1) or (len(ax[1].lines) > 1):
                ax[0].lines[-1].remove()
                ax[1].lines[-1].remove()

            # draw a line at the time of the chosen transit-event (from slider) on the top and bottom plot
            lver0 = ax[0].axvline(val, color = 'r', linewidth = 2)
            lver1 = ax[1].axvline(val, color = 'r', linewidth = 2)

            # if clicked on the image to select the time to zoom in on then update the plot in this function.
            ax[1].set_xlim([val - 1,val + 1])
            lver0.set_xdata(val)
            lver1.set_xdata(val)

            # also update the location of the slider (aesthetic reasons only)
            # Define the slider to change the transit-event time (and cut out region) - this might slow the code down as it re-plots the slider everytime but I can't find an efficient way to do this differently for now.
            fig.canvas.draw_idle()  # update figure

            # these values need to be global as are needed but the clicking and the slider events
            global transit_slider_ax
            global transit_slider

            transit_slider_ax.remove() # remove the old bar and plot a new one which has the colour bar in the right place
            transit_slider_ax  = fig.add_axes([0.25, 0.14, 0.65, 0.03])  # location of the slider
            transit_slider = Slider(transit_slider_ax, 'Transit', np.nanmin(alltimebinned), np.nanmax(alltimebinned), valinit=val, color='teal')

        # if clicked outside of the region activate the slider bar
        else:
            # Link the functions (actions) to the corresponding sliders
            transit_slider.on_changed(update_axis)
            scale_slider.on_changed(update_yaxis)
            transit_slider.on_changed(sliders_on_changed)
            #sector_slider.on_changed(update_plotting_region)

    # alternatively one can click on the image in order to change the zoom in region of the plot.
    fig.canvas.mpl_connect('button_press_event', onclick)

    # ---------------------

    # Define buttons which allow the user to interactively choose options:
    # run simple code, BLS, model the data usign Pyaneti, save the data and figures, creata DV report


    # only give the model option if pyaneti has been sucessfully installed
    if pyaneti_installed == True:
        var_ax = fig.add_axes([0.02, 0.52, 0.119, 0.19]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Hide plots', 'North', 'BLS', 'model', 'Save', 'Report'), (False, args.noshow, args.north, False, False, True, True))
    else:
        var_ax = fig.add_axes([0.02, 0.52, 0.119, 0.18]) # x, y, width, height
        save_var = CheckButtons(var_ax, ('Simple','Hide plots', 'North', 'BLS', 'Save', 'Report'), (False, args.noshow, args.north, False, True, True))


    aster_ax = fig.add_axes([0.02, 0.22, 0.119, 0.18]) # x, y, width, height
    aster_var = RadioButtons(aster_ax, ('transit', 'asteroseismology + transit', 'asteroseimology only'), active=0)


    # -------------------
    # make a check button to select the sectors to display in the top pannel
    # start_sec, end_sec, in_sec,
    global seclist_names
    seclist_names = tuple(['Sec {}'.format(i) for i in in_sec])
    init_secbuttons = tuple([False] * len(seclist_names))

    seclist_names = ('All',)+ seclist_names
    init_secbuttons  = (True,) + init_secbuttons

    box_add = (len(init_secbuttons) - 4) * 0.025

    # if we're looking at more than one sector give the option to change sector
    if len(init_secbuttons) > 2:

        if len(init_secbuttons) < 5:
            secax = plt.axes([0.02, 0.32, 0.1, 0.15])
        else:
            secax = plt.axes([0.02, 0.32 - (box_add), 0.1, 0.15 + box_add])

        save_secax= RadioButtons(secax, seclist_names, init_secbuttons)

        # make sure that the 'dots' check boxes are
        rpos = secax.get_position().get_points()
        fh = fig.get_figheight()
        fw = fig.get_figwidth()
        rscale = (rpos[:,1].ptp() / rpos[:,0].ptp()) * (fh / fw)
        for circ in save_secax.circles:
            circ.height /= rscale

        save_secax.on_clicked(update_plotting_region)
        sector_slider.on_changed(update_plotting_region)

    # Initial values for each option
    simple = False
    hide = args.noshow
    north = args.north
    BLS = False
    model = False
    save = True
    DV = True

    # function to get the status of each button and to save it
    def variables(label):
        status = save_var.get_status()
        simple = status[0]
        hide = status[1]
        north = status[2]
        BLS = status[3]

        if pyaneti_installed == True:
            model = status[4]
            save = status[5]
            DV = status[6]
        else:
            model = False
            save = status[4]
            DV = status[5]

    # ---------------

    # Add a set of radio buttons for changing the binning of the data
    binning_ax = fig.add_axes([0.02, 0.77, 0.1, 0.15]) # x, y, width, height
    binning_radios = RadioButtons(binning_ax, ('2', '5', '7', '10'), active=0)

    # this function accesses the binning functino.
    def binning_button(label):
        line_full_binned.set_xdata(binning(int(label))[0])
        line_full_binned.set_ydata(binning(int(label))[1])
        fig.canvas.draw_idle()

    binning_radios.on_clicked(binning_button)
    # ---------------

    # define the paramaters of the plot
    minf = np.nanmin(np.array(allflux))
    maxf = np.nanmax(np.array(allflux))
    height = maxf - minf

    ax[0].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[0].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[0].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[0].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[0].vlines(all_md, minf-1,minf + height*0.3 , colors = md_color, label = "Momentum Dump")

    ax[1].tick_params(axis="y",direction="inout", labelsize = main_plot_fontsize) #, pad= -20)
    ax[1].tick_params(axis="x",direction="inout", labelsize = main_plot_fontsize) #, pad= -17)
    ax[1].tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')
    ax[1].set_xlabel("Time (BJD-2457000)", fontsize = main_plot_fontsize)
    ax[1].set_ylabel("Normalised Flux", fontsize = main_plot_fontsize)
    ax[1].vlines(all_md, minf-1,minf + height*0.3, lw =1,  colors = md_color, label = "Momentum Dump")


    # ---------------
    # Create a label for the 'Sector Selection' box to clarify what it does
    plt.text(0.02, 1.1, "Sector Selection", fontsize=10, verticalalignment='center')
    # ---------------
    # Create a label for the 'settings' box to clarify what it does
    plt.text(0.02, -0.24, "Settings", fontsize=10, verticalalignment='center')

    # ---------------

    if len(init_secbuttons) > 3:
        # Create a label for the 'binning' box to clarify what it does
        plt.text(0.02, -1.9, "Binning Factors", fontsize=10, verticalalignment='center')
        # ---------------

    # define a box to enter the nickname of the target - this is useful if you analyse a lot of different candidates and you need a way of identifying them easily

    # ---------------
    # define a 'text box' that lets you enter the transit times

    initial_text = ""  # no initial text (the text box is empty)
    nick_name = []  # list of the entered transit times

    # function to store the entererd transit times.
    def submit(text):
        nick_name.append(text)

    axbox = plt.axes([0.25, 0.03, 0.22, 0.04]) # x, y, width, height
    text_box = TextBox(axbox, 'Memorable name (optional): ', initial=initial_text)
    text_box.on_submit(submit)
    # ---------------


    # Create a button to close the figure and more onto the next stage of the code.
    ebx = plt.axes([0.77, 0.03, 0.13, 0.04])
    exit = Button(ebx, 'Done', color=main_pot_color)

    # make button to exit the plot and continue with the code.
    # pop up warning if no transit time is entered


    def close(event):


        if len(transit_times) > 0 and (asteroseismology == False):
            plt.close('all')
        elif (len(transit_times) == 0) and (asteroseismology == True):
            print ("asteroseismic analysis only")
            transit_times = [-999]
            plt.close('all')
        else:
            print ("Must enter at least one transit time!")
            ax[1].text(0.3,0.45, "Please enter at least one transit time!", color = 'red', size=9, ha="center", transform=ax[1].transAxes)
            fig.canvas.draw()


    exit.on_clicked(close)
    # ---------------

    global ttxt
    # text that will be updated to list the 'wanted' transit event times
    ttxt = plt.text(-4.01,1.55, "Transit times: ", weight='bold')

    # Create a button to store the 'current' value of the slider
    stx = plt.axes([0.53, 0.03, 0.11, 0.04])
    store_val = Button(stx, 'Add time', color='gold')

    transit_times = []

    # function to store a transit event time using a button
    def storeval(time):

        storetime = (transit_slider.val)
        transit_times.append(round(storetime,2))
        global ttxt

        ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))

        plt.draw()


    store_val.on_clicked(storeval)
    # ---------------

    # Create a button to delete the last entered value
    delx = plt.axes([0.65, 0.03, 0.11, 0.04])
    del_val = Button(delx, 'Remove time', color='grey')

    transit_times = []

    # function to delete the last entry if using a button
    def deleval(time):
        global ttxt
        if len (transit_times)>0:
            del transit_times[-1]

            ttxt.set_text("Transit times: {} ".format(str(transit_times)[1:-1]))
            plt.draw()


    del_val.on_clicked(deleval)

    # ---------------

    plt.show()
    plt.close('all')

    # ------ END OF INTERACTIVE PART OF CODE ------

    # save the status of all the buttons
    end_status = save_var.get_status()

    if pyaneti_installed == True:
        simple = end_status[0]
        hide = end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = end_status[4]
        save = end_status[5]
        DV = end_status[6]

    else:
        simple = end_status[0]
        hide = end_status[1]
        north = end_status[2]
        BLS = end_status[3]
        model = False
        save = end_status[4]
        DV = end_status[5]



    # ---------------

    # update the arguments if they are different from the ones entered in the command line:
    if args.noshow != hide:
        args.noshow = hide

    if args.north != north:
        args.north = north

    # if the save button was not selected, then change the iput argument to false
    if save == False:
        args.save = save

    # check whether the user identified a nickname, or memorable name, for the candidate
    if len(nick_name) != 0: # if the box was left untouched do nothing and continue without setting a new nickname - will be ignored
        if len(nick_name[-1]) > 0: # [-1] in case the user pressed enter and then deleted the nickname again
            args.nickname = str(nick_name[-1])

    # ---------------
    # get the entered transit times and make sure that all the values are floats
    transit_list = transit_times

    if len(transit_list) > 0:
        print ("Transit times : {}".format(str(transit_list)[1:-1]))

    #  -----  BREW  ------
    brew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, [-111], args)
