
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
from LATTE import LATTEutils as utils
warnings.filterwarnings('ignore')

'''
Overview of LATTE scipts:

__main__.py      : Intitialises the parameters, what TIC ID, sector, checks for downloads of the data, FFI or not? 

LATTEutils.py    : All the functions needed to download data and text files, runs the interactive gui, all of the plotting and data handling. 

LATTE_DV.py      : Scipt to combine all of the results from LATTEbrew in order to generate a pdf data validation report.

LATTEbrew.py     : Calls the LATTEutils functions in turn in order to store the results and keeps track of what has been generated. Calls the LATTE_DV.py function to collate all the results.


This scipt has two functions, one for the short cadence data and one for the FFI data. 
These functions are either called from within __main__.py (if run with input target list), or from within the interact functions in LATTEutils.py (normal mode).

Brew keeps track if any part of the code crashes (e.g. the downloading of the TPF - so that these are omitted in the generation of the DV report).
'''

def brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, args):
	'''
	This function combines all the results from LATTE and calls all the different functions - 
	it makes the plots, saves them, runs the BLS model and the pyaneti model before making a PHT DV report (if this option is selected.) 
	
	Parameters
	----------
	tic  :   str
		target TIC ID
	indir  :  str
		path to directory where all the plots and data will be saved. 
	transit_list   : list
		list of the transit-like events
	simple   :   boolean
		whether or not to run the simple version
	BLS   :   boolean
		whether or not to run the BLS routine
	model   :   boolean
		whether or not to model the transit using pyaneti
	save   :	boolean
		whether or not to save the figures and data
	DV   :   boolean
		whether or not to write and save a DV report
	sectors_all  :   list
		all the sectors in which the target has been/ will be observed
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
	teff  :  float
		effective temperature of the tagret star (K)
	srad  :  float
		radius of the target star (solar radii)
	 ra	:	float 
		the right ascension of the target stars
	 dec	:   float 
		the declination of the target star

	'''

	# -------------------
	# SAVE THE DATA FILES 
	# -------------------

	if (save == True) or (DV == True):
		save = True
		# if this folder doesn't exist then create it. These are the folder where the images, data and reports for each TIC ID will be stored.
		newpath = '{}/{}'.format(indir,tic)
		
		if not exists(newpath):
			os.makedirs(newpath)

		# save the data used as a text file - these often come in use later for a quick re-analysis.
		with open('{}/{}/{}_data.txt'.format(indir, tic, tic), "w") as f:

			# get rid of nan values first using a mask
			good_mask = np.isfinite(np.array(alltime)) * np.isfinite(np.array(allflux)) * np.isfinite(np.array(allflux_err)) 
			alltime_ar = np.array(alltime)[good_mask]
			allflux_ar = np.array(allflux)[good_mask]
			allflux_err_ar = np.array(allflux_err)[good_mask]

			# save
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['time', 'flux', 'flux_err'])
			writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar))
		
		'''
		if the modelling option was also chose, save another data file with slightly different formatting to be called by Pyaneti.
		Pyaneti requires a very specific data format.
		
		Furhermore, in order for Pyaneti to run more efficiently (it has a Bayesian backend which scales with number of data points)
		we create a cutout of the times around the times of the marked transit events. 
		'''

		if model == True:
			with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "w") as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['#time', 'flux', 'flux_err'])

			# If the dip separations are too small, then don't create cut outs and save the whole dataset
			if (len(transit_list) > 1) and ((transit_list[1] - transit_list[0]) < 2): # if there are LOTS of transit events on short period (if so it's probably a TOI but let's keep it here as a condition)
				with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
					writer = csv.writer(f, delimiter='\t')
					writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar)) # save all the data

			# else create a cut out of the data around the time of the transit events
			else:
				for transit in transit_list:
					# save the data 
					# get rid of nan values first - this is used for the pyaneti code
					pyaneti_mask = (alltime_ar > (transit - 1)) * (alltime_ar < (transit + 1))

					with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
						writer = csv.writer(f, delimiter='\t')
						writer.writerows(zip(alltime_ar[pyaneti_mask],allflux_ar[pyaneti_mask],allflux_err_ar[pyaneti_mask]))


	# -----------------------------------
	#			START PLOTTING		  	 - calls functions from LATTEutils.py
	# -----------------------------------
	
	# create a plot of the fulllighcurves with the momentum dumps (MDs) marked and a zoom-in of the marked transits
	# this plit is saved but not shown (as already shown in the interact part fo the code)
	utils.plot_full_md(tic, indir, alltime,allflux,all_md,alltimebinned,allfluxbinned, transit_list, args)

	# Get a list of the sectors that have transit marked in them
	# this is so that we no longer have to loop through all of the sectors, and can focus on the ones which are important.
	transit_sec = utils.transit_sec(in_sec, start_sec, end_sec, transit_list)
	
	# -----------

	# plot how the centroids moved during the transit event
	utils.plot_centroid(tic, indir,alltime12, allx1, ally1, allx2, ally2, transit_list, args)
	# plot the background flux at the time of the transit event.
	utils.plot_background(tic, indir,alltime, allfbkg, transit_list, args)
	
	print ("Centroid and background plots... done.")
	# -----------

	# if the 'simple' option is chosen in the GUI, then the code will end here - this is designed to provide a quick analysis requiring no TPFs.
	if simple == True:
		print ("Simple option was selected, therefore end analysis here.")
		sys.exit('')

	# -----------

	# call function to extract the Target Pixel File information 
	# this is needed in order to extract the LCs in different aperture sizes.
	# the data is extracted using the open source Lightkurve package as they a built in function to extract LCs using different aperture sizes
	TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, tpf_list = utils.download_tpf_lightkurve(indir, transit_list, sectors, tic)
	
	# if the TPF wasn't corrupt then make the TPF files (only very ocassionally corrupt but don't want code to crash if it is corrrupt)
	if (TESS_unbinned_t_l != -111):

		tpf_corrupt = False
		# plot the LCs using two different aperture sizes.
		utils.plot_aperturesize(tic,indir,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, transit_list, args)
		
		print ("Aperture size plots... done.")
		# ------------
	
		'''
		Plot the average pixel brightness of the cut-out around the target star and the corresponding SDSS field of view.
		Both are oriented so that North is pointing upwards.
		The former also shows the nearby stars with TESS magnitude brighter than 17. Queried from GAIA using astroquery.
		The function returns the mass of the star (also output from astroquery)- this is a useful input for the Pyaneti modelling		
		'''
		if args.mpi == False:
			test_astroquery, _, _, mstar = utils.plot_TESS_stars(tic,indir, transit_list, transit_sec, tpf_list, args)
		else:
			test_astroquery, _, _, mstar = utils.plot_TESS_stars_not_proj(tic,indir, transit_list, transit_sec, tpf_list, args)

		# keep track of whether astroquery is working (sometimes the site is down and we don't want this to stop us from running the code)
		astroquery_corrupt = False

		if test_astroquery == -999:
			astroquery_corrupt = True
			print ("Star Aperture plots... failed.")
		
		else:
			print ("Star Aperture plots... done.")
		
		# ------------
		
		# Download the Target Pixel File using the raw MAST data - this comes in a different format as the TPFs extracted using Lightkurve
		# This data is then corrected using Principal Component Analysis is orderto get rid of systematics.
		X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = utils.download_tpf_mast(indir, transit_sec, transit_list, tic)
		
		# ------------
		
		'''
		plot the in and out of transit flux comparison.
		By default the images are NOT oriented north - this is because the reprojection takes longer to run and for a simple
		analysis to check whether the brightest pixel moves during the transit this is not required.
		The orientation towards north can be defined in the command line with '--north'.
		'''
		if args.north == True:
			utils.plot_in_out_TPF_proj(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, tpf_list, args)
			print ("In and out of aperture flux comparison with reprojection... done. ")
		
		else:
			utils.plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, args)
			print ("In and out of aperture flux comparison... done.")
		# ------------
		
		# For each pixel in the TPF, extract and plot a lightcurve around the time of the marked transit event.
		utils.plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list,t_list, T0_list, args)
		print ("Pixel level LCs plot... done.")
		# ------------
		
	else:
		tpf_corrupt = True
		mstar = 1 # need to define mstar otherwise pyaneti will complain - just make it one as an approximation.
		astroquery_corrupt = True
	
	# ------------
	# end of plots that require target pixel files
	# ------------

	# If more than one transit has been marked by the user, the LC is phase folded based on the period of the separation of the first two maarked peaks.
	# These plots are saved but do not feature in the DV report.
	if len (transit_list) > 1: # needs to know a period so can only do this if more than one transit has been marked.

		period = transit_list[1] - transit_list[0]
		t0 = transit_list[0] # time of the first marking
		
		# calculate the phase
		phased = np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in alltimebinned])

		fig, ax = plt.subplots(figsize=(5.55,5))
		ax.plot(phased, allfluxbinned, marker='.',color = 'k', alpha = 1, lw = 0, markersize = 4, label = 'None', markerfacecolor='k')
	
		#ax.plot(phased, allflux,marker='o',color = 'navy', alpha = 0.7, lw = 0, markersize = 2, label = 'binning = 7', markerfacecolor='white')
		plt.title("Phase folded LC")
		ax.set_xlabel("Phase (days)")
		ax.set_ylabel("Normalized Flux")
		plt.plot()


		if save == True:
			plt.savefig('{}/{}/{}_phase_folded.png'.format(indir, tic, tic), format='png')
		
		if args.noshow == False:
			plt.show()

		print ("Phase folded plot... done.")

	else:
		print ("\n Only one transit marked - therefore can't be phase folded. \n")
	# ------------

	'''
	Plot LCs of the six closest TESS target stars. This allows us to check whether the transit-like events 
	also appear in other nearby LCs which would be a sign that this is caused by a background event.
	'''

	# get the tic IDs of the six nearest stars
	ticids, distance, target_ra, target_dec = utils.nn_ticids(indir, transit_sec, tic)

	# download the data for these stars
	alltime_nn, allflux_nn, all_md_nn, alltimebinned_nn, allfluxbinned_nn,outtics,tessmag_list, distance = utils.download_data_neighbours(indir, transit_sec[0], ticids, distance)
	
	# plot the LCs
	utils.plot_nn(tic, indir,alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, transit_list, outtics, tessmag_list, distance, args)
	print ("Nearest neighbour plot... done.")
	# ------------

	# if the BLS option is chose, a BLS search is run. The LCs are first detrended and smoothed using a moving average. 
	# The corrected and uncorrected LCs are saves as a single plot for comparison and to verify that the correction worked well - saved but do not feature in the DV report. 
	if BLS == True:
		print ("Running BLS algorithm...", end =" ")
		bls_stats1, bls_stats2 = utils.data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, args)
		print ("done.")
	# ------------

	# SKIP FROM HERE....

	'''
	NOTE: CURRENTLY ONLY WORKS ON NORA'S COMPUTER - WILL BE AVAILABLE IN NEXT RELEASE SO PLEASE SCIP THIS PART OF THE CODE
	If the modelling option is selected (in the GUI), model the transit event using Pyaneti (Barragan et al 2018)
	which uses an Bayesian approach with an MCMC sampling to best fit and model the transit.
	The code runs slightly differently depending on whether one or multiple transits have been marked. 
	This is because with multiple transits the code has information about the possible orbital period.
	Need to ensure that the code has compiled correctly on the users computer. 
	'''

	# First check if Pyaneti is installed...
	if os.path.exists("{}/pyaneti_LATTE.py".format(syspath)):

		if model == True:
			print ("Running Pyaneti modelling - this could take a while so be patient...")
	
			transit_list_model =  ("{}".format(str(np.asarray(transit_list)))[1:-1]) # change the list into a string and get rid of the brackets
			# the code is usually run through the command line so call it using the os.system function.
			
			os.system("python3 {}/pyaneti_LATTE.py {} {} {} {} {} {} {}".format(syspath, tic, indir, syspath, mstar, teff, srad, transit_list_model))
	
	else:
		print ("Pyaneti has not been installed so you can't model anything yet. Contact Nora or Oscar for the LATTE version of the Pyaneti code.")
		model = False

	# ... UNTIL HERE
	# ------------


	# Finally, create a DV report which summarises all of the plots and tables.
	if DV == True: 
		from LATTE import LATTE_DV as ldv

		if BLS == True:
			ldv.LATTE_DV(tic, indir, syspath, transit_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, tpf_corrupt, astroquery_corrupt, FFI = False,  bls = True, model = model, mpi = args.mpi)
		else:
			ldv.LATTE_DV(tic, indir, syspath, transit_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, [0], [0], tpf_corrupt, astroquery_corrupt, FFI = False,  bls = False, model = model, mpi = args.mpi)


def brew_LATTE_FFI(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux_normal, allflux_small, allflux, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, args):
	
	'''
	This function that runs LATTE - makes the plots, saves them, runs the BLS model and the pyaneti model before 
	making a PHT DV report (if told to do so.)
	This function is very similar to brew_LATTE - 
	except that it is designed to analyse the FFIs and therefore the data download is different.
	
	Parameters
	----------
	tic  :   str
		target TIC ID
	indir  :  str
		path to directory where all the plots and data will be saved. 
	transit_list   : list
		list of the transit-like events
	simple   :   boolean
		whether or not to run the simple version
	BLS   :   boolean
		whether or not to run the BLS routine
	model   :   boolean
		whether or not to model the transit using pyaneti
	save   :	boolean
		whether or not to save the figures and data
	DV   :   boolean
		whether or not to write and save a DV report
	sectors_all  :   list
		all the sectors in which the target has been/ will be observed

	alltime  :  list
		times
	allflux_normal  :  list
		normalized flux extracted with the larger aperture (PCA corrected)
	allflux_small  : list
		normalized flux extracted with the smaller aperture (PCA corrected)
	allflux   : list 
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
	ra   :   float
		right ascension of the target star
	dec   :   float
		declination of the target star

	'''

	# ----------------------------------------
	#		   SAVE THE DATA FILES 
	# ----------------------------------------


	if (save == True) or (DV == True):
		save = True

		# if this folder doesn't exist then create it. These are the folder where the images, data and reports for each TIC ID will be stored.
		newpath = '{}/{}'.format(indir,tic)
		
		if not exists(newpath):
			os.makedirs(newpath)
		
		# save the data used as a text file - these often come in use later for a quick re-analysis.
		with open('{}/{}/{}_data.txt'.format(indir, tic, tic), "w") as f:
			
			# get rid of nan values first using a mask
			good_mask = np.isfinite(np.array(alltime)) * np.isfinite(np.array(allflux))
			alltime_ar = np.array(alltime)[good_mask]
			allflux_ar = np.array(allflux)[good_mask]
			allflux_err_ar = allflux_ar * 0.001

			# save data 
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['time', 'flux', 'flux_err'])
			writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar))

		'''
		if the modelling option was also chose, save another data file with slightly different formatting to be called by Pyaneti.
		Pyaneti requires a very specific data format.
		
		Furhermore, in order for Pyaneti to run more efficiently (it has a Bayesian backend which scales with number of data points)
		we create a cutout of the times around the times of the marked transit events. 
		'''

		if model == True:
			with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "w") as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['#time', 'flux', 'flux_err'])

			# if the dip separations are too small, then don't create cut outs and save the whole dataset
			if (len(transit_list) > 1) and ((transit_list[1] - transit_list[0]) < 2): # if there are LOTS of transit events on short period (if so it's probably a TOI but let's keep it here as a condition)
				with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
					writer = csv.writer(f, delimiter='\t')
					writer.writerows(zip(alltime_ar,allflux_ar,allflux_err_ar)) # save all the data

			# else create a cut out of the data around the time of the transit event
			else:
				for transit in transit_list:
					# save the data 
					# get rid of nan values first - this is used for the pyaneti code
					pyaneti_mask = (alltime_ar > (transit - 1)) * (alltime_ar < (transit + 1))

					with open('{}/{}/{}_data_pyaneti.dat'.format(indir, tic, tic), "a") as f:
						writer = csv.writer(f, delimiter='\t')
						writer.writerows(zip(alltime_ar[pyaneti_mask],allflux_ar[pyaneti_mask],allflux_err_ar[pyaneti_mask]))


	# ------------------------------------------
	#			   START PLOTTING
	# ------------------------------------------
	
	# create a plot of the fulllighcurves with the momentum dumps (MDs) marked and a zoom-in of the marked transits
	utils.plot_full_md(tic, indir, alltime, allflux, all_md, alltime, allflux, transit_list, args)

	# get the sectors that have a transit marked in them
	transit_sec = utils.transit_sec(in_sec,start_sec, end_sec, transit_list)
	
	# -----------

	# plot the background flux at the time of the transit event.
	utils.plot_background(tic, indir, allfbkg_t, allfbkg, transit_list, args)
	print ("Background plots... done.")
	# -----------

	if simple == True:
		print ("Simple option was selected, therefore end analysis here.")
		sys.exit('')
	# -----------

	# plot the LC using different aperture sizes
	utils.plot_aperturesize(tic,indir,alltime, alltime, alltime, allflux_normal, allflux_normal, allflux_small, transit_list, args)
	
	print ("Aperture size plots... done.")
	# -----------

	'''
	Plot the average pixel brightness of the cut-out around the target star and the corresponding SDSS field of view.
	Both are oriented so that North is pointing upwards.
	The former also shows the nearby stars with TESS magnitude brighter than 17. Queried from GAIA using astroquery.
	The function returns the TESS magnitude, effective temperature (K), the radius of the star (solar radii)and the 
	mass of the star (solar mass) (also output from astroquery)- this is a useful input for the Pyaneti modelling		
	'''

	if args.mpi == False:
		tessmag, teff, srad, mstar = utils.plot_TESS_stars(tic,indir,transit_list, transit_sec, tpf_list, args)
	else:
		tessmag, teff, srad, mstar = utils.plot_TESS_stars_not_proj(tic,indir,transit_list, transit_sec, tpf_list, args)

	# keep track of whether astroquery is working (sometimes the site is down and we don't want this to stop us from running the code)
	astroquery_corrupt = False

	if tessmag == -999:
		astroquery_corrupt = True
		print ("Star Aperture plots... failed.")
	
	else:
		print ("Star Aperture plots... done.")
		

	# -----------

	# If more than one transit has been marked by the user, the LC is phase folded based on the period of the separation of the first two maarked peaks.
	# These plots are saved but do not feature in the DV report.

  # needs to know a period so can only do this if more than one transit has been marked.
	if len(transit_list) > 1:
		# find the period
		period = transit_list[1] - transit_list[0]
		t0 = transit_list[0] # time of the first marking
		
		# calculate the phase
		#phased = np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in alltime])

		# phased plot where odd and even are different colours
		phased = (np.array(alltime) - t0) % period
		phased2 = (np.array(alltime) - (t0)) % (period*2)
		
		index = phased>0.5*period
		index2 = phased2 >0.5*(period *2)
		
		phased[index] -= period
		phased2[index2] -= (period *2)
	
		fig, ax = plt.subplots(figsize=(5.55,5))

		ax.plot(phased,np.array(allflux), 'k.', markersize=5)
		ax.plot(phased2,np.array(allflux), 'r.', markersize=5)
	
		#ax.plot(phased, allflux,marker='o',color = 'navy', alpha = 0.7, lw = 0, markersize = 2, label = 'binning = 7', markerfacecolor='white')
		plt.title("Phase folded LC")
		ax.set_xlabel("Phase (days)")
		ax.set_ylabel("Normalized Flux")
		plt.plot()

		if save == True:
			plt.savefig('{}/{}/{}_phase_folded.png'.format(indir, tic, tic), format='png')
		
		if args.noshow == False:
			plt.show()

		print ("Phase folded plot... done.")

	else:
		print ("\n Only one transit marked - therefore can't be phase folded. \n")	
	# ------------

	'''
	In order to create the plots that compare the in and out of transit average flux 
	we need to create TPF masks for the in and out of transit. 
	This couldn't be done at the time of the extraction of the TPF arrays as the transit times 
	had not yet been defined. 
	'''

	oot_list = [] # out of transit
	intr_list = []  # in transit

	X1_list_n = []
	X4_list_n = []
	bkg_list_n = []
	apmask_list_n = []
	arrshape_list_n = []
	t_list_n = []

	tpf_filt_list_n = []
	tpf_list_n = []

	T0_list = []  # list of the transit times - make a new list to ensure that thet are in the same order.
	

	for T0 in transit_list:
		for idx,t in enumerate(t_list):
			# these times were found to give good approximations for out of transit and in transit
			# Note: this will NOT be 'perfect' for all events, as it depends on the transit duration which the code doesnt' know at this time but it's a good approximation.
			if (T0 > np.nanmin(t)) and (T0 < np.nanmax(t)): # if T0 is within the time... (should be)
				
				oot_list.append((abs(T0-np.array(t)) < 0.55) * (abs(T0-np.array(t)) < 0.3)) 
				intr_list.append(abs(T0-np.array(t)) < 0.1)
				
				X1_list_n.append(X1_list[idx])
				X4_list_n.append(X4_list[idx])
				bkg_list_n.append(bkg_list[idx])
				apmask_list_n.append(apmask_list[idx])
				arrshape_list_n.append(arrshape_list[idx])
				t_list_n.append(t)

				tpf_filt_list_n.append(tpf_filt_list[idx])
				tpf_list_n.append(tpf_list[idx])

				T0_list.append(T0)

	X1_list	 	= 	X1_list_n
	X4_list	 	= 	X4_list_n
	bkg_list		= 	bkg_list_n
	apmask_list 	= 	apmask_list_n
	arrshape_list   =   arrshape_list_n
	t_list 			=   t_list_n
	tpf_filt_list 	= 	tpf_filt_list_n
	tpf_list 		=   tpf_list_n
	
	# ------------


	# Plot the in and out of transit flux comparison 
	# By default the images are NOT orented north - this is because the reprojectio takes longer to run and for a simple
	# analysis to chekc whether the brightest pixel moves during the transit this is not required.
	# The orientation towards north can be defined in the command line with '--north'.

	if args.north == True:
		utils.plot_in_out_TPF_proj(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, tpf_list, args)
		print ("In and out of aperture flux comparison with reprojection... done. ")
	else:
		utils.plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, args)
		print ("In and out of aperture flux comparison... done.")
	# ------------

	# For each pixel in the TPF, extract and plot a lightcurve around the time of the marked transit event.
	utils.plot_pixel_level_LC(tic, indir, X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, args)
	print ("Pixel level LCs plot... done.")
	# ------------

	# If the BLS option is chose, a BLS search is run. The LCs are first detrended and smoothed using a moving average. 
	# The corrected and uncorrected LCs are saves as a single plot for comparison and to verify that the correction worked well - saved but do not feature in the DV report. 
	if BLS == True:
		print ("Running BLS")
		bls_stats1, bls_stats2 = utils.data_bls_FFI(tic, indir, alltime, allflux, args)
	# ------------

	'''
	If the modelling option is selected (in the GUI), model the transit event using Pyaneti (Barragan et al 2018)
	which uses an Bayesian approach with an MCMC sampling to best fit and model the transit.
	The code runs slightly differently depending on whether one or multiple transits have been marked. 
	This is because with multiple transits the code has information about the possible orbital period.
	Need to ensure that the code has compiled correctly on the users computer. 

	SKIP FOR NOW - this will become available in the next version of LATTE
	'''

	# first check if Pyaneti is even installed...

	if os.path.exists("{}/pyaneti_LATTE.py".format(syspath)):

		if model == True:
			print ("Running Pyaneti modelling - this could take a while so be patient...")
			
			# the code is usually run through the command line so call it using the os.system function.
			transit_list_model =  ("{}".format(str(np.asarray(transit_list)))[1:-1]) # change the list into a string and get rid of the brackets
			os.system("python3 {}/pyaneti_LATTE.py {} {} {} {} {} {} {}".format(syspath, tic, indir, syspath, mstar, teff, srad, transit_list_model))
	
	else:
		print ("Pyaneti has not been installed so you can't model anything yet. Ask Nora or Oscar for the LATTE version of the Pyaneti code.")
		model = False
	
	# SKIP until HERE
	
	# ------------
	# Finally, create a DV report which summarises all of the plots and tables.
	if DV == True: 
		from LATTE import LATTE_DV as ldv

		if BLS == True:
			ldv.LATTE_DV(tic, indir, syspath, transit_list, sectors_all, ra, dec, tessmag, teff, srad, bls_stats1, bls_stats2, False, astroquery_corrupt, FFI = True, bls = True, model = model, mpi = args.mpi)
		else:
			ldv.LATTE_DV(tic, indir, syspath, transit_list, sectors_all, ra, dec, tessmag, teff, srad, [0], [0], False, astroquery_corrupt, FFI = True, bls = False, model = model, mpi = args.mpi)

	else:
		print ("\n  Complete! \n ")
