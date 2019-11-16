# last updated 12 Nov 2019

#imports
import os
import ast
import csv
import sys
import warnings
import matplotlib
import numpy as np
import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
from argparse import ArgumentParser

#froms from standards
from tkinter import simpledialog
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter

#froms from non-standards
from os.path import exists
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons

#custom modules to import
import LATTEutils as lu
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


		lu.plot_full_md(tic, indir, alltime,allflux,allline,alltimebinned,allfluxbinned, peak_list)


	# --------------
	# START PLOTTING
	# --------------

	## get the sectors that have a transit marked in them
	peak_sec = lu.peak_sec(in_sec,start_sec, end_sec, peak_list)
	
	lu.plot_centroid(tic, indir,alltime12, allx1, ally1, allx2, ally2, peak_list, save = save, show = show)
	lu.plot_background(tic, indir,alltime, allfbkg, peak_list, save = save, show = show)
	
	if simple == True:
		print ("Simple option was selected, therefore end analysis here.")
		exit([])


	TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l = lu.tpf_data(indir, sectors, tic)
	lu.plot_aperturesize(tic,indir,TESS_unbinned_t_l, TESS_binned_t_l, small_binned_t_l, TESS_unbinned_l, TESS_binned_l, small_binned_l, peak_list, save = save, show = show)
	
	# nearby stars
	lu.plot_TESS_stars(tic,indir,peak_list, peak_sec, save = save, show = show)

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

	else:
		print ("\n Only one transit marked - therefore can't be phase folded. \n")
		
	
	X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list, t_list, T0_list, tpf_filt_list = lu.download_data_tpf(indir, peak_sec, peak_list, tic)
	

	# plot the in and out of transit flux comparison
	lu.plot_in_out_TPF(tic, indir, X4_list, oot_list, t_list, intr_list, T0_list, tpf_filt_list, save = save, show = show)


	# plot one lightcurve for each pixel extracted around the time of the transit
	lu.plot_pixel_level_LC(tic, indir,X1_list, X4_list, oot_list, intr_list, bkg_list, apmask_list, arrshape_list,t_list, T0_list, save = save, show = show)
	

	# nearest neighbour light curve comparison plot
	ticids, distance, target_ra, target_dec = lu.nn_ticids(indir, peak_sec, tic)
	alltime_nn, allflux_nn, allline_nn, alltimebinned_nn, allfluxbinned_nn,outtics,tessmag_list, distance = lu.download_data_neighbours(indir, peak_sec[0], ticids, distance)
	lu.plot_nn(tic, indir,alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, peak_list, outtics, tessmag_list, distance, save = save, show = show)

	if BLS == True:
		print ("running BLS")
		bls_stats1, bls_stats2 = lu.data_bls(tic, indir, alltime, allflux, allfluxbinned, alltimebinned, save = save, show = show)
		

	# run pyaneti - this takes a while to run so carefully cnosider whether you want to run this
	
	# first check if Pyaneti is even installed...
	if os.path.exists("./pyaneti_LATTE.py"):

		if model == True:
			print ("Running Pyaneti modelling - this could take a while so be patient...")
	
			peak_list_model =  ("{}".format(str(np.asarray(peak_list)))[1:-1]) # change the list into a string and get rid of the brackets
			os.system("python pyaneti_LATTE.py {} {} {}".format(tic, indir, peak_list_model))
	
	else:
		print ("Pyaneti has not been installed so you can't model anything yet. Ask Nora or Oscar for the LATTE version of the Pyaneti code.")
		model = False


	if DV == True: 
		import LATTE_DV as ldv

		if BLS == True:
			ldv.LATTE_DV(tic, indir, peak_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, bls = True, model = model)
		else:
			ldv.LATTE_DV(tic, indir, peak_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, [0], [0], bls = False, model = model)


def interact_LATTE(tic, indir, sectors_all, sectors):
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
	alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = lu.download_data(indir, sectors, tic)
	print ("Done.\n")
	
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
		N	   = len(alltime)
		n	   = int(np.floor(N/binfac)*binfac)
		X	   = np.zeros((2,n))
		X[0,:]  = alltime[:n]
		X[1,:]  = allflux[:n]
		Xb	  = rebin(X, (2,int(n/binfac)))
		
		time_binned	= Xb[0]
		flux_binned	= Xb[1]
	
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
	
	plt.text(4.5, -3.15, "(Press Enter, then close window)", fontsize=10, verticalalignment='center')
	plt.text(0.05, 1.1, "Binning Factor", fontsize=10, verticalalignment='center')
	
	initial_text = ""
	
	transit_times = []

	def submit(text):
		ydata = eval(text)
		transit_times.append(ydata)
	
	axbox = plt.axes([0.25, 0.05, 0.65, 0.03])
	text_box = TextBox(axbox, 'Enter transit-event times', initial=initial_text)
	text_box.on_submit(submit)
	
	plt.show()
	
	end_status = save_var.get_status()

	simple = end_status[0]
	BLS = end_status[1]
	model = end_status[2]
	save = end_status[3]
	DV = end_status[4]


	if len(transit_times) == 0:
		print ("\n WARNING: You can't continue without entering a transit-time.\n")
		print ("Exit.\n")
		sys.exit()
	
	if type(transit_times[-1]) == tuple:
		peak_list = list(transit_times[-1])
		peak_list = [float(i) for i in peak_list]
	
	else:
		peak_list = [float(i) for i in transit_times]
		peak_list = [peak_list[-1]] #if you entered it twice

	print ("Transits you have entered:	  {}   \n".format(str(peak_list))[1:-1])
	print ("Check that these are the transits that you want")
	
	
	# END OF INTERACTIVE PART OF CODE

	#  -----  BREW  ------
	brew_LATTE(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, show = args.noshow)


if __name__ == '__main__':
	ap = ArgumentParser(description='Lightcurve Analysis Tool for Transiting Exoplanets')
	ap.add_argument('--new-data', action='store_true')
	ap.add_argument('--targetlist',type=str, help='the link to the target file list', default='no')
	ap.add_argument('--noshow', action='store_false', help='if you want to NOT show the plots write --noshow in the command line')
	ap.add_argument('--o', action='store_true', help='if you call this old files will be overwriten in the non-interactive version')
	
	args = ap.parse_args()

	# This is the only thing that needs to be changes in this program to make it link to your own Computer
	# ----------

	indir = "/Users/Nora/Documents/research/TESS/planethunters/LATTE"  # CHANGE THIS
	
	# ----------------
	# ---- START -----
	# ----------------

	if args.new_data != False: 
		# This checks to download the data reference files
		# This will only run if you tell the program to run this - needs to be run the first time that one wants to download data
		# Also run if new data is released - the program will check what is already available and only download new things.
	
		# ---- REFERENCE FILES DOWNLOAD -----
		lu.data_files(indir)
		lu.nn_files(indir)
		lu.TOI_TCE_files(indir)
		# ----
	
	# ------  INTERACTIVE VERSION -------
	# -----------------------------------
	'''
	NOTE: this requires you to have Tkinter installed. Tkinter does not work with certain new Mac operating systems.
	In order to run LATTE with an input list and not interatcively, state the path to the csv file when running the program.
	csv file must have format: "TICID, sectors, transits, BLS, model" - see example.
	'''
	
	if args.targetlist == 'no': 


		# make a GUI interface with TKinter
		ROOT = tk.Tk()
		ROOT.withdraw()
		
	
		# first prompt window which asks for TIC ID	
		tic = simpledialog.askstring(title="TIC",
										  prompt="Enter TIC ID:")
		
		# returns all of the sectors in which TESS observed the given TIC id - this uses TESS-point
		sectors_all = lu.tess_point(indir, tic) 
		
	
		#  Ask to define the sectors
		sectors = simpledialog.askstring(title="Sectors",
										  prompt="TIC {} was observed in sector(s):\n {}. \n Enter the ones you wish to look at (e.g. 1,4) or 'all' for all of them.".format(tic, sectors_all))
		
		ROOT.quit()
		ROOT.destroy()
		
		if len(sectors) == 0:
			sectors = 'all'
		
		if sectors != 'all':
			sectors = sectors.split(",")
			sectors = [int(i) for i in sectors]
		
		print ("\n")
	
		if sectors == 'all':
			print ("Will look at  sector(s) ({}) \n".format(sectors_all))
		else:
			print ("Will look at available sectors: {} \n".format(sectors))
	
	


		interact_LATTE(tic, indir, sectors_all, sectors)
		
		# Done


	# --------------------------------
	#		RUN WITH INPUT FILE
	# --------------------------------
	
	else:

		try:
			targetlist = pd.read_csv("{}".format(args.targetlist)) # If a path is defined, open the input file
		except:
			print ("This target list can't be found. Check the path you have given and the name and format of the file.")

		# process each target individually. 
		for index, row in targetlist.iterrows():
			
			# ---- INPUT PARAMETERS ------

			tic = str(row['TICID'])

			# Check whether this file already exists
			if (os.path.exists("{}/{}".format(indir, tic)))  and (args.o != True): 
				print ("This file already exists therefore SKIP. To overwrite files run this code with --o in the command line.")
				continue
			# -------------

			# --- WHAT SECTORS IS IT OBSERVED IN? ---
			sectors_all = lu.tess_point(indir, tic) 
			sectors_in = row['sectors']
			
			try:
				sectors_in = ast.literal_eval(sectors_in)
				if (type(sectors_in) == int) or (type(sectors_in) == float):
					sectors_in = [sectors_in]
				else:
					sectors_in = list(sectors_in) 
				
				# Sucessfully enteredt sectors
				# check that the target was actually observed in the stated sector
				sectors = list(set(sectors_in) & set(sectors_all))
				if len(sectors) == 0:
					print ("The target was not observed in the sector(s) you stated ({}). \
							Therefore take all sectors that it was observed in: {}".format(sectors, sectors_all))
					sectors =sectors_all
			except:
				sectors = sectors_all
	


			# ---- IF NO TRANSIT MARKED RUN WITH INTERFACE ----
			if type(row['transits']) == float:
				interact_LATTE(tic, indir, sectors_all, sectors)



			else:
				peak_list_in = (row['transits'])
				peak_list = ast.literal_eval(peak_list_in)
				
				# convert the input transit times and sectors into peak_list in the form of a list
				if type(peak_list) == float:
					peak_list = [peak_list]
				else:
					peak_list = list(peak_list)
				
	
				BLS_in = row['BLS']
				model_in = row['model']
	
				# ---- do you want to run a BLS? ----
				if (BLS_in == 'yes') or (BLS_in == 'True') or (BLS_in == 'true') or (BLS_in == '1'):
					BLS = True
				else:
					BLS = False
				
				# ---- do you want to model this transit? ----
				if (model_in == 'yes') or (model_in == 'True') or (model_in == 'true') or (model_in == '1'):
					model = True
				else:
					model = False
				
				simple = False  # we don't want to run the simple version - that's option is simply to do quick test
				save = True  # always save the files - no point running this if you're not going to save them
				DV = True   # we'll make a DV report for all of them
				

				# -------------------
				# DOWNLOAD DATA 
				# -------------------
	
				alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = lu.download_data(indir, sectors, tic)
				
				# --------------------------
				#	  START BREWING ....
				# --------------------------
	
				brew_LATTE(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, show = False)


# End.

