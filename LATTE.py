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
from glob import glob
#import matplotlib.pyplot as plt
from argparse import ArgumentParser

#froms from standards
from tkinter import simpledialog
#from matplotlib.ticker import AutoMinorLocator
#from matplotlib.ticker import FormatStrFormatter

#froms from non-standards
from os.path import exists

#custom modules to import
import LATTEutils as utils
import LATTEbrew as brew  
warnings.filterwarnings('ignore')


if __name__ == '__main__':
	ap = ArgumentParser(description='Lightcurve Analysis Tool for Transiting Exoplanets')
	ap.add_argument('--new-data', action='store_true')
	ap.add_argument('--tic', type=str, help='the target TIC id, if blank will be prompted to enter', default = 'no')
	ap.add_argument('--sector', type=str, help='the sector(s) to look at', default = 'no')
	ap.add_argument('--targetlist', type=str, help='the link to the target file list', default='no')
	ap.add_argument('--noshow', action='store_false', help='if you want to NOT show the plots write --noshow in the command line')
	ap.add_argument('--o', action='store_true', help='if you call this old files will be overwriten in the non-interactive version')
	#ap.add_argument('--mstar', type=float, help='the mass of the star in case it is known', default=1)
	ap.add_argument('--nickname', type=str, help='give the target a memorable name', default='no')
	ap.add_argument('--FFI', action='store_true', help='is this an FIIs?')

	args = ap.parse_args()

	# ------- CHANGE THIS --------
	#indir = "/Users/Nora/Documents/research/TESS/planethunters/LATTE"  # CHANGE THIS
	
	indir = "./LATTE_output"

	# ----------------------------

	if not os.path.exists("{}".format(indir)):
		os.makedirs(indir)

	if not os.path.exists("{}/data".format(indir)):
		os.makedirs("{}/data".format(indir))

	# ----------------
	# ---- START -----
	# ----------------

	if args.new_data != False: 
		# This checks to download the data reference files
		# This will only run if you tell the program to run this - needs to be run the first time that one wants to download data
		# Also run if new data is released - the program will check what is already available and only download new things.
	
		# ---- REFERENCE FILES DOWNLOAD -----
		utils.data_files(indir)
		utils.nn_files(indir)
		utils.TOI_TCE_files(indir)
		utils.momentum_dumps_info(indir)
		# ----
	
	# ------  INTERACTIVE VERSION -------
	# -----------------------------------
	'''
	NOTE: this requires you to have Tkinter installed. Tkinter does not work with certain new Mac operating systems.
	In order to run LATTE with an input list and not interatcively, state the path to the csv file when running the program.
	csv file must have format: "TICID, sectors, transits, BLS, model" - see example.
	'''
	
	if args.targetlist == 'no': 
		#if both the sectors and the tic ID are already entered then TKinter does not need to be loaded
		if args.tic != 'no' and args.sector != 'no':
			tic = str(args.tic)
			sectors = str(args.sector)

			if args.FFI == True:
				newpath = '{}/{}'.format(indir,tic)
				# if this folder doesn't exist then create it...
				if not exists(newpath):
					os.makedirs(newpath)

			sectors_all, ra, dec = utils.tess_point(indir, tic) 
		
		else:
			# make a GUI interface with TKinter
			ROOT = tk.Tk()
			ROOT.withdraw()
			
			# has the tic already been defined?
			if args.tic == 'no':

				# first prompt window which asks for TIC ID	
				tic = simpledialog.askstring(title="TIC",
												  prompt="Enter TIC ID:")
			else:
				tic = str(args.tic)

			if args.FFI == True:
				newpath = '{}/{}'.format(indir,tic)
				# if this folder doesn't exist then create it...
				if not exists(newpath):
					os.makedirs(newpath)

			# has the sector already been defined? 
			if args.sector == 'no':
				# returns all of the sectors in which TESS observed the given TIC id - this uses TESS-point
				sectors_all, ra, dec = utils.tess_point(indir, tic) 
		
				#  Ask to define the sectors
				sectors = simpledialog.askstring(title="Sectors",
												  prompt="TIC {} was observed in sector(s):\n {}. \n Enter the ones you wish to look at (e.g. 1,4) or 'all' for all of them.".format(tic, sectors_all))
			else:
				sectors_all, ra, dec = utils.tess_point(indir, tic) 
				sectors = str(args.sector)

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
	
		if args.FFI == False:
			utils.interact_LATTE(tic, indir, sectors_all, sectors, ra, dec, args.noshow)  # the argument of whether to shos the images or not 
		else:
			utils.interact_LATTE_FFI(tic, indir, sectors_all, sectors, ra, dec, args.noshow)
		
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

			tic = str(int(row['TICID']))

			existing_files = glob("{}/*{}*".format(indir, tic))
			# Check whether this file already exists
			if (len(existing_files) > 0)  and (args.o != True): 
				print ("This file already exists therefore SKIP. To overwrite files run this code with --o in the command line.")
				continue

			# -------------

			# --- WHAT SECTORS IS IT OBSERVED IN? ---

			sectors_all = utils.tess_point(indir, tic) 
			sectors_in = row['sectors']
			
			try:
				sectors_in = ast.literal_eval(sectors_in)
				if (type(sectors_in) == int) or (type(sectors_in) == float):
					sectors_in = [sectors_in]
				else:
					sectors_in = list(sectors_in) 
				
				# Sucessfully entered sectors
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
				utils.interact_LATTE(tic, indir, sectors_all, sectors, args.noshow)


			else:
				peak_list_in = (row['transits'])
				peak_list = ast.literal_eval(peak_list_in)
				
				# convert the input transit times and sectors into peak_list in the form of a list
				print (peak_list)
				print (type(peak_list))

				if (type(peak_list) == float) or (type(peak_list) == int):
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
	
				alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = utils.download_data(indir, sectors, tic)
				
				# --------------------------
				#	  START BREWING ....
				# --------------------------
	
				brew.brew_LATTE(tic, indir, peak_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, allline, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, show = False)

	# end by changing the name of the folder to include the nicknane if so defined in the input functions
	# this allows users to keep track of their targets more easily. We name our candidates after pastries. 
	
	if (not args.nickname == 'no') and (args.FFI == True):
		os.system("mv {}/{} {}/FFI_{}_{}".format(indir, tic, indir, tic, args.nickname))

	elif not args.nickname == 'no':
		os.system("mv {}/{} {}/{}_{}".format(indir, tic, indir, tic, args.nickname))
	
	elif args.FFI == True:
		os.system("mv {}/{} {}/FFI_{}".format(indir, tic, indir, tic))
	
# End.

