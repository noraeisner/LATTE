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

from os.path import exists
from tkinter import simpledialog
from argparse import ArgumentParser

#custom modules to import
from LATTE import LATTEbrew as brew 
from LATTE import LATTEutils as utils
#sys.tracebacklimit = 0
warnings.filterwarnings('ignore')

'''
NOTE: a warning message currently appears when the code is exited (i.e. at the end of the code). 
This is due to an error in the current verion of matplolib with the interactive widgets but has been adressed
in future releases of matplolib (see https://github.com/matplotlib/matplotlib/issues/13660)
The above link shows how to change the matplolib cbook __init__ file in order to make it disapear.

# REQUIRES MATPLOLIB 3.2 (there is a bug in 3.1 the current stable version - as of January 2020)). 
 --- pip install matplotlib==3.2.0rc1 (still in testing?)


Overview of LATTE scipts:

__main__.py	  : Intitialises the parameters, what TIC ID, sector, checks for downloads of the data, FFI or not? 

LATTEutils.py	: All the functions needed to download data and text files, runs the interactive gui, all of the plotting and data handling. 

LATTE_DV.py	  : Scipt to combine all of the results from LATTEbrew in order to generate a pdf data validation report.

LATTEbrew.py	 : Calls the LATTEutils functions in turn in order to store the results and keeps track of what has been generated. Calls the LATTE_DV.py function to collate all the results.

'''


# -------------------- START ---------------------
# ------------------------------------------------

if __name__ == '__main__':
	# these are all the input parameters that the user can (but doesn't have to) chose to identify in the command line.
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
	ap.add_argument('--mpi', action = 'store_true', help ='enter if the program is run with mpi')

	args = ap.parse_args()


	# ------------------------------------------------
	# Check what the current path is - when the program is first downloaded the path is set to 'no/path/set/yet' and the user is automatically prompted to change'no/path/set/yet'
	# Get the path to where the package is installed. This is where the configuration file and the images (for the DV report are stored)

	syspath = str(os.path.abspath(utils.__file__))[0:-14]
	# ----------

	def yes_or_no():
		'''
		Yes/No command line option to verify that the user really wants to change the output/input path
		'''
		print ('\n \n WARNING: if you have already downloded the input files (with --new-data) then these will remain in the location set by your previous path, so you will have to redowload the data (not recommended) or move the data to the new location set by this path. \n \n ')
	
		reply = str(input('Are you sure that you want to change the path?' + '(yes/no): '))
	
		if (reply == 'y') or (reply == 'yes') or (reply == 'yep') or (reply == 'yeah'):
			return True
	
		else: # if anything else is entered assume that this is a 'no' and continue with the old path
			return False	 
	
	# ------------------------------
	# ---- INITIALISE THE CODE -----
	# ------------------------------

	'''
	This first part of the code 1) ensures that an input and output path has been defined (and allows the user to change it if necessary)
	The path links to where the curl scipt files are stored and to where the output files will be stores. 2) Creates the output file if it doens't
	already exist, and 3) downloads the text files if run for the first time or if stated to download. 
	'''

	# check whether a path already exists
	if not os.path.exists("{}/_config.txt".format(syspath)):

		# if it doesn't exist ask the user to put it in the command line
		indir = input("\n \n No output path has been set yet. \n \n Please enter a path to save the files (e.g. ./LATTE_output or /Users/yourname/Desktop/LATTE_output) : " )

		#try three times to find the output path
		intry = 0
		worked = True

		while (intry < 3 ) and worked:
			if not os.path.exists("{}".format(indir)):
				indir = input("\n \n This path does not exist on your computer, please enter a valid  path: " )
	
				intry += 1

			else:
				worked = False
	
		if intry == 3:
			print ("You have entered an invalid path 3 times. \n Make sure that the path exists on your computer and that you have access to it before trying again.")
			sys.exit('')

		# SAVE the new output path
		with open("{}/_config.txt".format(syspath),'w') as f:
			f.write(str(indir))
		
		print("\n New path: " + indir)
	
		# this is also the first time that the program is being run, so download all the data that is required.
		print ("\n Download the text files required ... " )
		print ("\n Only the manifest text files (~325 M) will be downloaded and no TESS data." )
		print ("\n This step may take a while but luckily it only has to run once... \n" )

		# ------------------------------------------------
		#check whether the chosen (output) directory already exists, and if it doesn't create the directory.
		
		if not os.path.exists("{}/data".format(indir)):
			os.makedirs("{}/data".format(indir))
	
		# ------------------------------------------------

		# ----- REFERENCE FILES DOWNLOAD -----
		utils.data_files(indir)
		utils.tp_files(indir)
		utils.TOI_TCE_files(indir)
		utils.momentum_dumps_info(indir)
		
		# -----

	# if the user chooses to redefine the path
	elif args.new_path == True: 
	
		reply = yes_or_no() # double check that they really want to do that.
	
		if reply == True:
			indir = input("\n \n Please enter a path to save the files (e.g. ./LATTE_output or /Users/yourname/Desktop/LATTE_output) : " )

			#try three times to find the output path
			intry = 0
			worked = True
	
			while (intry < 3 ) and worked:
				if not os.path.exists("{}".format(indir)):
					indir = input("\n \n This path does not exist on your computer, please enter a valid  path: " )
		
					intry += 1
	
				else:
					worked = False
		
			if intry == 3:
				print ("You have entered an invalid path 3 times. \n Make sure that the path exists on your computer and that you have access to it before trying again.")
				sys.exit('')
	
			# SAVE the new output path
			with open("{}/_config.txt".format(syspath),'w') as f:
				f.write(str(indir))	
			
			print("\n New path: " + indir)
	
			if not os.path.exists("{}/data".format(indir)):
				os.makedirs("{}/data".format(indir))


		else:
			with open("{}/_config.txt".format(syspath), 'r') as f:
				indir = str(f.readlines()[-1])
				
			print ("LATTE will continue to run with the old path: {}".format(indir))

			#try three times to find the output path
			intry = 0
			worked = True
	
			while (intry < 3 ) and worked:
				if not os.path.exists("{}".format(indir)):
					indir = input("\n \n This old path does not exist on your computer, please enter a valid  path: " )
		
					intry += 1
	
				else:
					worked = False
		
			if intry == 3:
				print ("You have entered an invalid path 3 times. \n Make sure that the path exists on your computer and that you have access to it before trying again.")
				sys.exit('')
	
			if intry > 0:
				# SAVE the new output path
				with open("{}/_config.txt".format(syspath),'w') as f:
					f.write(str(indir))	

			if not os.path.exists("{}/data".format(indir)):
				os.makedirs("{}/data".format(indir))
				

	else:
		with open("{}/_config.txt".format(syspath), 'r') as f:
			indir = str(f.readlines()[-1])
	

		#try three times to find the output path
		intry = 0
		worked = True
	
		while (intry < 3 ) and worked:
			if not os.path.exists("{}".format(indir)):
				indir = input("\n \n Please enter a path to save the files (e.g. ./LATTE_output or /Users/yourname/Desktop/LATTE_output) : " )
		
				intry += 1
	
			else:
				worked = False
		
		if intry == 3:
			print ("You have entered an invalid path 3 times. \n Make sure that the path exists on your computer and that you have access to it before trying again.")
			sys.exit('')
	
		if intry > 0:
			# SAVE the new output path
			with open("{}/_config.txt".format(syspath),'w') as f:
				f.write(str(indir))	


		if not os.path.exists("{}/data".format(indir)):
			os.makedirs("{}/data".format(indir))
	
		# ------------------------------------------------

	'''
	Check whether to download the data reference files
	This will only run if you tell the program to run this (with the args)- needs to be run the first time that one wants to download data
	This should also be called if new TESS data is released 
	The program will check what data has already been downloaded and only download new data files.
	'''

	if (args.new_data != False) and (os.path.exists("{}/_config.txt".format(syspath))): 

		# ----- REFERENCE FILES DOWNLOAD -----
		utils.data_files(indir)
		utils.tp_files(indir)
		utils.TOI_TCE_files(indir)
		utils.momentum_dumps_info(indir)
		# -----
	


	# -----------  INTERACTIVE VERSION  ------------
	# ---------------------------------------------
	
	'''
	This section starts the interactive versio of the code:

	1) asks to enter the TIC ID
	2) identifies in which sectors the target was observed and asks the user to identify which sectors should be analysed

	NOTE: this part of the code requires you to have Tkinter installed. Tkinter currently does not work with certain new Mac operating systems.
	In order to run LATTE with an input list and not interatcively, state the path to the csv file when running the program.
	csv file must have format: "TICID, sectors, transits, BLS, model" - see example.
	The TIC ID and Sectors to look at can also be stated in the command line as an argument to save time
	The interactive tool is for ease of use and to avoid having to understand how to use the command line.
	'''
	
	# Check whether the a target list has been defined in the command line. If so, the code will not ask for a TIC ID or a sector 
	# as these would be listed in the input file and the code will run automatically. 
	if args.targetlist == 'no': 

		# Check whether the tic ID and the sector have already been entered
		# If both the sectors and the TIC ID are already entered then TKinter does not need to be loaded
		if args.tic != 'no' and args.sector != 'no':
			tic = str(args.tic)
			sectors = str(args.sector)

			# Check whether we are looking at an FFI
			# The FFI information needs to be stored straight away so a folder needs to be created to store them. 
			# The folder for the non-FFIs is created later after the user choses whether to 'save' data or not.
			if args.FFI == True:
				newpath = '{}/{}'.format(indir,tic)

				# if this folder doesn't exist then create it...
				if not exists(newpath):
					os.makedirs(newpath)

			# --------
			# Run a function called TESS-point. This returns the sectors in which the target has 
			# been observed as well as the RA and DEC of the target. 
			sectors_all, ra, dec = utils.tess_point(indir, tic) 
			# --------

		# --------------------------------
		# If either the TIC or the sectors or both have not already been identified, run Tkinter (interactive tools)
		
		else:

			# make a GUI interface with TKinter
			ROOT = tk.Tk()
			ROOT.withdraw()

			# if this hasn't already been identified as an FFI through the command line, give the option to chose this when entering the TIC ID
			if args.FFI == False:
				# -----------
				class TICprompt(simpledialog.Dialog):
				
					def body(self, master):
					
						# make a text box for the TIC ID
						tk.Label(master, text="Enter TIC ID:").grid(row=0)
						self.e1 = tk.Entry(master)
						self.e1.grid(row=0, column=1)
						
						# make a check button with the option to run this in FFI mode
						self.FFIbox = tk.IntVar()
						self.answer = tk.Checkbutton(master,text="Check for FFI mode", variable=self.FFIbox)
						self.answer.grid(row=1, column=1,  columnspan=2)
				
					def apply(self):
						# make these global variables so that they can be used outside of this class and applied to the rest of the program
						global tkTIC
						global tkFFI
						
						ROOT.form=(self.FFIbox.get())
						tkTIC = (self.e1.get())
						tkFFI =  (self.FFIbox.get())
				# -----------

				TICprompt(ROOT)
				# make the TIC a string
				tic = str(tkTIC)

				# If the FFI button was checked, change the FFI argument 

				if tkFFI == 1:
					args.FFI = True

			else: # if this is a FFI, don't give the FFI option again
				# has the TIC already been defined? If not, ask for the TIC ID with a text box. 
				if args.tic == 'no':
	
					# load first prompt window which asks for TIC ID	
					tic = simpledialog.askstring(title="TIC",
													  prompt="Enter TIC ID:")
				else:
					tic = str(args.tic)

			# -----------

			# Is this an FFI? - if it is, need to create a folder to store the FFI data in (one for each TIC ID)
			if args.FFI == True:
				newpath = '{}/{}'.format(indir,tic)
				# if this folder doesn't exist then create it...
				if not exists(newpath):
					os.makedirs(newpath)

			# has the sector already been defined? 
			if args.sector == 'no':
				# returns all of the sectors in which TESS observed the given TIC id - this uses TESS-point

				sectors_all, ra, dec = utils.tess_point(indir, tic) 
		
				# The user is shown a list of the sectors in which the target is observed 
				# and asked to state which ones they want to assess

				sectors = simpledialog.askstring(title="Sectors",
												  prompt="TIC {} was observed in sector(s):\n {}. \n Enter the ones you wish to look at (e.g. 1,4) or 'all' for all of them.".format(tic, sectors_all))
			else:
				# still need to run tess point even if the targets are already defined as we need to check whether target appears in the stated sector - sanity check
				sectors_all, ra, dec = utils.tess_point(indir, tic) 
				sectors = str(args.sector)

			# close the TKinter windows.
			ROOT.quit()
			ROOT.destroy()
		
		# if no sector is defined or the word 'all' is written in the box, analyse all of the given sectors.
		if len(sectors) == 0:
			sectors = 'all'

		
		# if not all of them are chosen, convert the input list (given as a string) into a python readable list
		if sectors != 'all':
			sectors = sectors.split(",")
			sectors = [int(i) for i in sectors]
		
		print ("\n")
		
		# print out the information that has been chosen to the command line.
		if sectors == 'all':
			print ("Will look at sector(s): {} \n".format(str(sectors_all)[1:-1]))
			sectors = sectors_all
		else:
			print ("Will look at sector(s):  {} \n".format(str(sectors)[1:-1]))
	

		# -------------------------------------
		# ---- OPEN INTERACTIVE MATPLOTLIB ----
		# -------------------------------------

		''' Start up LATTE interactive where the transit times can be chosen manually 
		 	this works differently for FFI data and target files as the data has a different format. 
		 	At this point the code continues in the LATTE_utils scipt. 
		 	This scipt (__main__.py) is only to set up the parameters for the target and to intitialise the code. 
		'''

		if args.FFI == False:
			utils.interact_LATTE(tic, indir, syspath, sectors_all, sectors, ra, dec, args)  # the argument of whether to shos the images or not 
		else:
			utils.interact_LATTE_FFI(tic, indir, syspath, sectors_all, sectors, ra, dec, args)
		
		# Done.

	# ---------------------------------------
	# ---------------------------------------
	#			RUN WITH INPUT FILE
	# ---------------------------------------
	#The below code is executed if the 'input target list' option has been chose - this is defined in the command line. 
	#The code is run with input targetlist - either with phase fold information or with transit time information.

	else:
		try:
			targetlist = pd.read_csv("{}".format(args.targetlist)) # If a path is defined, open the input file
		except:
			print ("This target list can't be found. Check the path you have given and the name and format of the file.")

		# check what kind of target list it is - with phase fold information or with period and t0 information
		if 'period' in targetlist.columns:
			pp = True # pp = phase fold
		else:
			pp = False

		# process each target individually
		for index, row in targetlist.iterrows():
			
			# ---- INPUT PARAMETERS ------
			try:
				tic = str(int(row['TICID']))
			except:
				continue

			# check whether this file already exist
			# if it already exists it will only be overwritten if --o function has been enabled to avoid file loss.
			existing_files = glob("{}/*{}*".format(indir, tic))
			
			if (len(existing_files) > 0)  and (args.o != True): 
				print ("This file already exists therefore SKIP. To overwrite files run this code with --o in the command line.")
				failed_tics = [] #keep track of the TIC IDs that failed to complete
				continue


			# --- WHAT SECTORS WAS IT OBSERVED IN? ---

			# load tess-point to identify the sectors tht it was observed in.
			# need to know whether the target actually was observed in the stated sector. 
			sectors_all, ra, dec = utils.tess_point(indir, tic)
			
			# convert the input list of sectors (string) to a list of numbers.
			try:
				sectors_in = str(int(row['sectors']))
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
					sectors = sectors_all
			except:
				sectors = sectors_all

			failed_tics = [] #keep track of the TIC IDs that failed to complete


			# ---- IF NO TRANSIT MARKED RUN WITH INTERFACE ----
			if (pp == False) and (type(row['transits']) == float):
				utils.interact_LATTE(tic, indir, syspath, sectors_all, sectors, ra, dec, args)

			else:

				if pp == False: # if the transits were identified and no period information
					
					# extract the information of the tiems of the transit like events. 
					transit_list_in = (row['transits'])
					transit_list = ast.literal_eval(transit_list_in)
					
					# convert the input transit times and sectors into transit_list in the form of a list
					
					if (type(transit_list) == float) or (type(transit_list) == int):
						transit_list = [transit_list]
					else:
						transit_list = list(transit_list)
				
				# if the user entered a T0 adn a period, the code will calcualte the times if the transit like events and use that. 

				else:
					# get up to 3 markings - more than that will just be cluttered and will take too long - this can be changed later if more or less are desired. 
					transit_list = []
					period = row['period']
					t0 = row['t0']
					for i in range(0,4):
						transit = t0 + (i*period)
						if transit < (t0 + 20):
							transit_list.append(float(transit))

				# extract the other parameters from the input file.
				BLS_in = row['BLS']
				model_in = row['model']
				FFI_in = row['FFI']
	
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
				
				# ---- run in FFI mode? ----
				if (FFI_in == 'yes') or (FFI_in == 'True') or (FFI_in == 'true') or (FFI_in == '1'):
					args.FFI = True


				# --- other SETTINGS ---
				simple = False  # we don't want to run the simple version - that's option is simply to do quick test
				save = True  # always save the files - no point running this if you're not going to save them
				DV = True   # we'll make a DV report for all of them
				args.noshow = True  # don't show these, just save them
				args.auto = True # the code will automatically choose it's own apertures in this mode

				if args.FFI == False:
					try:
						# ----------------------------------------
						# 			 DOWNLOAD DATA 
						# ----------------------------------------
						alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = utils.download_data(indir, sectors, tic)
					
						# ----------------------------------------
						#			   START BREWING ....
						# ----------------------------------------
						
						brew.brew_LATTE(tic, indir, syspath,transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, args)
					except:

						failed_tics.append(tic)
						print ("TIC {} failed to complete. Continue anyway. You can find a list of the failed IDs stored in the output folder.".format(tic))
						continue
				else:
					try:
						# ----------------------------------------
						# 			 DOWNLOAD DATA 
						# ----------------------------------------

						alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t,start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list = utils.download_data_FFI(indir, sectors, syspath, sectors_all, tic, args)
	
						# ----------------------------------------
						#			   START BREWING ....
						# ----------------------------------------
						brew.brew_LATTE_FFI(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime0, allflux_list, allflux_small, allflux0, all_md, allfbkg, allfbkg_t, start_sec, end_sec, in_sec, X1_list, X4_list, apmask_list, arrshape_list, tpf_filt_list, t_list, bkg_list, tpf_list, ra, dec, args)
					except:
						failed_tics.append(tic)
						print ("TIC {} failed to complete. Continue anyway. You can find a list of the failed IDs stored in the output folder.".format(tic))
						continue	


		# ----- ALL FILES IN INFILE PROCESSED -----
		# save a list of the failed TIC IDs

		if len(failed_tics) > 0: # if any of the targets failed to complete
			# check if there is already a file with failed tic-ids - we don't want to overwrite it but create a new one
			failed_outpath = "{}/failed_files".format(indir)
		
			# --------
			if not os.path.exists(failed_outpath): # if this folder doesn't already exist, make it
				os.makedirs(failed_outpath)
			# --------
		
			failedfiles = np.sort(glob("{}/failed_tics*".format(failed_outpath)))
			
			if len(failedfiles) > 0:
				last_number = int((failedfiles[-1].split('_')[-1][0:-4]))
			else:
				last_number = 1
		
			with open('{}/failed_tics_{}.csv'.format(failed_outpath,str(last_number+1)), "w") as f:
				# save
				writer = csv.writer(f)
				for val in failed_tics:
						writer.writerow([val])


	# ---------------------------------------------------------------
	# ------------------------- LAST STEP ---------------------------
	# ---------------------------------------------------------------

	# end by changing the name of the folder to include the nicknane if defined
	# this allows users to keep track of their targets more easily e.g. Planet Hunters TESS candidates are named after pastries.
		
	if not args.nickname == 'noname':
		os.system("mv {}/{} {}/{}_{}".format(indir, tic, indir, tic, args.nickname))

	print ("\n  Complete! \n ")

# End.

