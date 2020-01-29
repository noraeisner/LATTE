#from __future__ import print_function, absolute_import, division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from glob import glob
import pandas as pd
from matplotlib import cm

#from bls import BLS
import seaborn as sb
import os
from astroquery.mast import Catalogs

from astropy.table import Table
from os.path import basename, exists

import time
import ast
import csv
from argparse import ArgumentParser

import LATTE
from LATTE import LATTEutils, LATTEbrew
import sys

# get the system path
syspath = str(os.path.abspath(LATTEutils.__file__))[0:-14]


# set up for supercomputer
try:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	mpi_rank = comm.Get_rank()
	mpi_size = comm.Get_size()
	with_mpi = True

except ImportError:
	mpi_rank = 0
	mpi_size = 1
	with_mpi = False

mpi_root = 0
# ---------

# --- IMPORTANT TO SET THIS ----
out = 'pipeline'  # or TALK or 'pipeline'
fitting = 'no'  # to fit this needs to say 'yes'

detrending = 'no' # if ='yes' then apply exotrending before fitting wit pyaneti
ttran = 0.1

# -------------------

def process(indir, sectors, tic, transit_list, args):

	# download the data 
	alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad = LATTEutils.download_data(indir, sectors, tic)

	# in my input file the the times start at 0 for each sector so I need the line below
	#transit_list = list(np.array(transit_list) + np.nanmin(alltime))
	# ------------

	simple = False
	BLS = False
	model = False
	save = True
	DV = True

	sectors_all, ra, dec = LATTEutils.tess_point(indir, tic)

	# generate the plots
	try:
		LATTEbrew.brew_LATTE(tic, indir, syspath, transit_list, simple, BLS, model, save, DV, sectors, sectors_all, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, args)
		tp_downloaded = True
	except:
		# see if it made any plots - often it just fails on the TPs as they are very large
		if exists("{}/{}/{}_fullLC_md.png".format(indir,tic,tic)):
			print ("couldn't download TP but continue anyway")
			tp_downloaded = False
		else:
			mnd = {}
			mnd['TICID'] = -99
			return mnd

	# check again whether the TPs downloaded - depending on where the code failed it might still have worked. 
	if exists("{}/{}/{}_aperture_size.png".format(indir,tic,tic)):
		tp_downloaded = True
	else:
		tp_downloaded = False
		print ("code ran but no TP -- continue anyway")
		
	# -------------
	# check whether it's a TCE or a TOI

	# TCE -----
	lc_dv = np.genfromtxt('{}/data/tesscurl_sector_all_dv.sh'.format(indir), dtype = str)

	TCE_links = []

	for i in lc_dv:
		if str(tic) in str(i[6]):
			TCE_links.append(i[6])

	if len(TCE_links) == 0:
		TCE = " - "
		TCE = False
	else:
		TCE_links = np.sort(TCE_links)
		TCE_link = TCE_links[0]  # this link should allow you to acess the MAST DV report
		TCE = True
	# TOI -----
	TOI_planets = pd.read_csv('{}/data/TOI_list.txt'.format(indir), comment = "#")
	
	TOIpl = TOI_planets.loc[TOI_planets['TIC'] == float(tic)]
		
	if len(TOIpl) == 0:
		TOI = False

	else:
		TOI = True
		TOI_name = (float(TOIpl["Full TOI ID"]))
	
	# -------------

	# return the tic so that it can be stored in the manifest to keep track of which files have already been produced
	# and to be able to skip the ones that have already been processed if the code has to be restarted.

	mnd = {}
	mnd['TICID'] = tic
	mnd['MarkedTransits'] = transit_list
	mnd['Sectors'] = sectors_all
	mnd['RA'] = ra
	mnd['DEC'] = dec
	mnd['SolarRad'] = srad
	mnd['TMag'] = tessmag
	mnd['Teff'] = teff
	mnd['thissector'] = sectors

	# make empty fields for the test to be checked
	if TOI == True:
		mnd['TOI'] = TOI_name
	else:
		mnd['TOI'] = " "

	if TCE == True:
		mnd['TCE'] = "Yes"
		mnd['TCE_link'] = TCE_link
	else:
		mnd['TCE'] = " "
		mnd['TCE_link'] = " "


	mnd['EB'] = " "
	mnd['Systematics'] = " "
	mnd['TransitShape'] = " "
	mnd['BackgroundFlux'] = " "
	mnd['CentroidsPositions'] = " "
	mnd['MomentumDumps'] = " "
	mnd['ApertureSize'] = " "
	mnd['InoutFlux'] = " "
	mnd['Keep'] = " "
	mnd['Comment'] = " "
	mnd['starttime'] = np.nanmin(alltime)

	return mnd

# ------------------------------------------------
# ------------------------------------------------


if __name__ == '__main__':
	ap = ArgumentParser(description='Script to pLot TESS LCs for Zooniverse project')
	ap.add_argument('--sector', type=int, help='The sector to analyse', default=None)
	ap.add_argument('--number', type=int, help='The number of randomly selected subjects (default all in indir folder)', default=None)

	# these just have to be here
	#ap.add_argument('--sector', type=str, help='the sector(s) to look at', default = 'no')
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

	sector = args.sector
	sectors = [args.sector]

	# make sure that they don't show...
	args.noshow = True
	args.mpi = True

	indir = 'path/to/input/file'


	df = pd.read_csv('link to input file goes here')

	# get rid of duplicates (why would there be suplicate...?)
	TIC_wanted0 = list(set(df['TIC_ID']))  # get rid of duplicates
	TIC_wanted = TIC_wanted[0:args.number] # take the top XX as defined in input


	## -----------
	if mpi_rank == mpi_root:
		
		print ("with MPI?:  {}".format(with_mpi))
		print ("MPI size = {}".format(mpi_size))

		nlc = len(TIC_wanted)

		print ("nlc length: {}".format(nlc))
		print ('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector)))

		if exists('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector))):
			print("Existing manifest file found, will skip previously processed LCs and append to end of manifest file")
		
		else:
			print("Creating new manifest file")
			metadata_header = ['TICID', 'Marked Transits', 'Sectors', 'RA', 'DEC', 'Solar Rad', 'TMag', 'Teff', 'thissector', 'TOI', 'TCE', 'TCE link', 'EB', 'Systematics', 'Background Flux', 'Centroids Positions','Momentum Dumps','Aperture Size','In/out Flux','Keep','Comment', 'starttime']


			with open('{}/manifest_sec{}.csv'.format(str(indir), str(sector)), 'w') as f: # save in the photometry folder
				writer = csv.writer(f, delimiter=',')
				writer.writerow(metadata_header)


		## Without MPI or running with a single node
		## =========================================
		if (not with_mpi) or (mpi_size==1) or (nlc==1):
			print ("not with MPI")
			for f in TIC_wanted:
				# check the existing manifest to see if I've processed this file!
				manifest_table = pd.read_csv('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector)))

				# get a list of the current URLs that exist in the manifest
				urls_exist = manifest_table['TICID']

				if not np.isin(f,urls_exist):

					# get the transit time list 
					transit_list = ast.literal_eval(((df.loc[df['TIC_ID'] == f]['db_peak']).values)[0])
					
					res = process(indir, sectors, f, transit_list, args)

					if res['TICID'] == -99:
						print ('something went wrong')
						continue
					# make sure the file is opened as append only
					with open('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector)), 'a') as f: # save in the photometry folder

						writer = csv.writer(f, delimiter=',')

						metadata_data = [res['TICID']]
						metadata_data.append(res['MarkedTransits'])
						metadata_data.append(res['Sectors'])
						metadata_data.append(res['RA'])
						metadata_data.append(res['DEC'])
						metadata_data.append(res['SolarRad'])
						metadata_data.append(res['TMag'])
						metadata_data.append(res['Teff'])
						metadata_data.append(res['thissector'])
						metadata_data.append(res['TOI'])
						metadata_data.append(res['TCE'])
						metadata_data.append(res['TCE_link'])
						metadata_data.append(res['EB'])
						metadata_data.append(res['Systematics'])
						metadata_data.append(res['BackgroundFlux'])
						metadata_data.append(res['CentroidsPositions'])
						metadata_data.append(res['MomentumDumps'])
						metadata_data.append(res['ApertureSize'])
						metadata_data.append(res['InoutFlux'])
						metadata_data.append(res['Keep'])
						metadata_data.append(res['Comment'])
						metadata_data.append(res['starttime'])
						
						writer.writerow(metadata_data)
				else:
					print('TIC {} already done - skipping'.format(f))

		else:
			## Master node
			## -----------

			if mpi_rank == 0:
				free_workers = list(range(1,mpi_size))
				active_workers = []
				n_finished_items = 0

				while TIC_wanted or active_workers:
					## Send a file
					while TIC_wanted and free_workers:
						f = TIC_wanted.pop()

						# check the existing manifest to see if I've processed this file!
						manifest_table = pd.read_csv('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector)))
						 # get a list of the current URLs that exist in the manifest
						urls_exist = manifest_table['TICID']
						if not (np.isin(f, urls_exist)):
							w = free_workers.pop()
							comm.send(f, dest=w, tag=0)
							active_workers.append(w)
							print('Sending TIC {} to worker {}'.format(str(f), w))
						else:
							print('TIC {} already done - skipping'.format(str(f)))
							n_finished_items += 1

					## Receive the results
					for w in active_workers:
						if comm.Iprobe(w, 2):
							res = comm.recv(source=w, tag=2)
							if res['TICID'] == -99:
								print ('something went wrong')
							else:
								print('Worker {} finished processing TIC {}'.format(w, res['TICID']))
								with open('{}/manifest_CTC_sec{}.csv'.format(str(indir), str(sector)), 'a') as f: # save in the photometry folder
									writer = csv.writer(f, delimiter=',')

									metadata_data = [res['TICID']]
									metadata_data.append(res['MarkedTransits'])
									metadata_data.append(res['Sectors'])
									metadata_data.append(res['RA'])
									metadata_data.append(res['DEC'])
									metadata_data.append(res['SolarRad'])
									metadata_data.append(res['TMag'])
									metadata_data.append(res['Teff'])
									metadata_data.append(res['thissector'])
									metadata_data.append(res['TOI'])
									metadata_data.append(res['TCE'])
									metadata_data.append(res['TCE_link'])
									metadata_data.append(res['EB'])
									metadata_data.append(res['Systematics'])
									metadata_data.append(res['BackgroundFlux'])
									metadata_data.append(res['CentroidsPositions'])
									metadata_data.append(res['MomentumDumps'])
									metadata_data.append(res['ApertureSize'])
									metadata_data.append(res['InoutFlux'])
									metadata_data.append(res['Keep'])
									metadata_data.append(res['Comment'])
									metadata_data.append(res['starttime'])
									writer.writerow(metadata_data)

							free_workers.append(w)
							active_workers.remove(w)
							n_finished_items += 1

				print ('Completed {} of {} files'.format(n_finished_items,nlc))
				if n_finished_items < nlc:
					print ('Failed on {} files'.format(nlc-n_finished_items))
				for w in free_workers:
					comm.send(-1, dest=w, tag=0)


	## Worker node
	## -----------
	else:
		while True:
			filename = comm.recv(source=mpi_root, tag=0)
			if filename == -1:
				break

			
			transit_list = ast.literal_eval(((df.loc[df['TIC_ID'] == filename]['db_peak']).values)[0])

			res = process(indir, sectors, filename, transit_list, args)
			comm.send(res, dest=mpi_root, tag=2)









