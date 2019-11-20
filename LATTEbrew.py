
import os
import ast
import csv
import sys
import warnings
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser

#froms from standards
from tkinter import simpledialog
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter

#froms from non-standards
from os.path import exists

#custom modules to import
import LATTEutils as utils

warnings.filterwarnings('ignore')



def brew_LATTE(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, show = False):
	'''
	This function that runs LATTE - makes the plots, saves them, runs the BLS model and the pyaneti model before 
	making a PHT DV report (if told to do so.)
	
	Parameters
    ----------
	tic  :   str
		target TIC ID
	indir  :  str
		path to directory where all the plots and data will be saved. 
	peak_list   : list
		list of the transit-like events
	simple   :   boolean
		whether or not to run the simple version
	BLS   :   boolean
		whether or not to run the BLS routine
	model   :   boolean
		whether or not to model the transit using pyaneti
	save   :    boolean
		whether or not to save the figures and data
	DV   :   boolean
		whether or not to write and save a DV report
	sectors_all  :   list
		all the sectors in which the target has been/ will be observed
	show  :  list
		the sectors which will be analysed

    Returns
    -------
    	runs the brew_LATTE code...
	'''
	# -------------------
	# SAVE THE DATA FILES 
	# -------------------

	if (save == True) or (DV == True):
		save = True

		newpath = '{}/{}'.format(indir,tic)
		# if this folder doesn't exist then create it...
		if not exists(newpath):
			os.makedirs(newpath)

		with open('{}/{}/{}_data.txt'.format(indir, tic, tic), "w") as f:
			# save the data 
			# get rid of nan values first - this is used for the pyaneti code
			
			good_mask = np.isfinite(np.array(alltime)) * np.isfinite(np.array(allflux)) * np.isfinite(np.array(allflux_err)) 

			alltime_ar = np.array(alltime)[good_mask]
			allflux_ar = np.array(allflux)[good_mask]
			allflux_err_ar = np.array(allflux_err)[good_mask]

			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['time', 'flux', 'flux_err'])
			writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar))
		
		if model == True:
			with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "w") as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['#time', 'flux', 'flux_err'])


			if (len(peak_list) > 1) and ((peak_list[1] - peak_list[0]) < 2): # if there are LOTS of transit events on short period (if so it's probably a TOI but let's keep it here as a condition)
				with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
					writer = csv.writer(f, delimiter='\t')
					writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar)) # save all the data

			# else create a cut out of the data around the time of the transit event
			else:
				for transit in peak_list:
					# save the data 
					# get rid of nan values first - this is used for the pyaneti code
					pyaneti_mask = (alltime_ar > (transit - 1)) * (alltime_ar < (transit + 1))

					with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
						writer = csv.writer(f, delimiter='\t')
						writer.writerows(zip(alltime_ar[pyaneti_mask],allflux_ar[pyaneti_mask],allflux_err_ar[pyaneti_mask]))


		utils.plot_full_md(tic, indir, alltime,allflux,allline,alltimebinned,allfluxbinned, peak_list)


	# --------------
	# START PLOTTING
	# --------------

	## get the sectors that have a transit marked in them
	peak_sec = utils.peak_sec(in_sec,start_sec, end_sec, peak_list)
	
	utils.plot_centroid(tic, indir,alltime12, allx1, ally1, allx2, ally2, peak_list, save = save, show = show)
	utils.plot_background(tic, indir,alltime, allfbkg, peak_list, save = save, show = show)
	print ("Centroid and background plots... done")
	if simple == True:
		print ("Simple option was selected, therefore end analysis here.")
		exit([])


	TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l = utils.tpf_data(indir, sectors, tic)
	utils.plot_aperturesize(tic,indir,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, peak_list, save = save, show = show)
	print ("Aperture size plots... done")
	# nearby stars
	utils.plot_TESS_stars(tic,indir,peak_list, peak_sec, save = save, show = show)
	print ("Star Aperture plots... done")
	# plot the phase folded
	if len (peak_list) > 1:

		period = peak_list[1] - peak_list[0]
		t0 = peak_list[0]
		
		fig, ax = plt.subplots(figsize=(10,6))
		phased = np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in alltimebinned])
		
		ax.plot(phased, allfluxbinned,marker='o',color = 'navy', alpha = 0.7, lw = 0, markersize = 2, label = 'binning = 7', markerfacecolor='white')
		plt.title("Phase folded LC")
		ax.set_xlabel("Phase (days)")
		ax.set_ylabel("Normalized Flux")
		plt.plot()

		if save == True:
			plt.savefig('{}/{}/{}_phase_folded.png'.format(indir, tic, tic), format='png')
		
		if show == True:
			plt.show()

		print ("Phase folded plot... done")

	else:
		print ("\n Only one transit marked - therefore can't be phase folded. \n")
		
	
	X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = utils.download_data_tpf(indir, peak_sec, peak_list, tic)
	

	# plot the in and out of transit flux comparison
	utils.plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, save = save, show = show)
	print ("In and out of aperture flux comparison... done")

	# plot one lightcurve for each pixel extracted around the time of the transit
	utils.plot_pixel_level_LC(tic, indir,X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list,t_list, T0_list, save = save, show = show)
	print ("Pixel level LCs plot... done")

	# nearest neighbour light curve comparison plot
	ticids, distance, target_ra, target_dec = utils.nn_ticids(indir, peak_sec, tic)
	alltime_nn, allflux_nn, allline_nn, alltimebinned_nn, allfluxbinned_nn,outtics,tessmag_list, distance = utils.download_data_neighbours(indir, peak_sec[0], ticids, distance)
	utils.plot_nn(tic, indir,alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, peak_list, outtics, tessmag_list, distance, save = save, show = show)
	print ("Nearest neighbour plot... done")

	if BLS == True:
		print ("Running BLS")
		bls_stats1, bls_stats2 = utils.data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, save = save, show = show)
		

	# run pyaneti - this takes a while to run so carefully cnosider whether you want to run this
	
	# first check if Pyaneti is even installed...
	if os.path.exists("./pyaneti_LATTE.py"):

		if model == True:
			print ("Running Pyaneti modelling - this could take a while so be patient...")
	
			peak_list_model =  ("{}".format(str(np.asarray(peak_list)))[1:-1]) # change the list into a string and get rid of the brackets
			os.system("python3 pyaneti_LATTE.py {} {} {}".format(tic, indir, peak_list_model))
	
	else:
		print ("Pyaneti has not been installed so you can't model anything yet. Ask Nora or Oscar for the LATTE version of the Pyaneti code.")
		model = False


	if DV == True: 
		import LATTE_DV as ldv

		if BLS == True:
			ldv.LATTE_DV(tic, indir, peak_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, bls = True, model = model)
		else:
			ldv.LATTE_DV(tic, indir, peak_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, [0], [0], bls = False, model = model)
