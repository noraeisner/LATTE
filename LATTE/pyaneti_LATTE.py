#-----------------------------------------------------------
#                         pyaneti.py
#                     Main pyaneti file
#                   Barragan O, March 2016
#-----------------------------------------------------------

#Load libraries
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sys



#-------------------------------------------------------------
#                   INITIALIZATION
#-------------------------------------------------------------

#Load the input file
#You have to run the program as ./pyaneti star_name
#tic = str(sys.argv[1])
#peak_list = str(sys.argv[2])

tic = str(sys.argv[1])
star = tic  # confusing but some places require the star to be called the tic and vice versa...

indir = str(sys.argv[2])

syspath = str(sys.argv[3])

mstar = str(sys.argv[4])

teff = str(sys.argv[5])

srad = str(sys.argv[6])


## --------------------
## now that we know where these file are - fownload the pyaneti files that we need
sys.path.append('{}/pyaneti_LATTE/'.format(syspath))
import pyaneti as pti #FORTRAN module
## --------------------

# find the transit list from the input params
peak_list = []

for i in range(7,len(sys.argv)): # the first agrument is the name of the program, the second is the tic ID, the third the directory, the fourth mstar, fifth teff and sixth rad
	peak_list.append(float(sys.argv[i]))

print ("-----------------------")
print (peak_list)
print ("-----------------------")

#Read the file with all the python functions
#execfile('pyaneti_LATTE/src/todo-py.py')
exec(open('{}/pyaneti_LATTE/src/todo-py.py'.format(syspath)).read())

#Read the file with the default values
#execfile('pyaneti_LATTE/src/default.py')
exec(open('{}/pyaneti_LATTE/src/default.py'.format(syspath)).read())

# -------------------
# PRIORS THAT CHANGE  - LATTE adaption
# -------------------

if len(peak_list) == 1:
	is_single_transit = True   # TOGGLE
else:
	print ("Pyaneti: this is a multi transit so we can constrain the orbital period...")
	is_single_transit = False 


T0b = peak_list[0]        # input variable

el = 100

min_t0  = [T0b - el*0.005]
max_t0  = [T0b + el*0.005]

if is_single_transit == False:
	Pb  = abs(peak_list[1] - peak_list[0])          # input varibale for period if there is one.
	print ("Period estimate: {} days".format(Pb))
	
	# only if not single transit
	min_P   = [Pb - el*0.001]
	max_P   = [Pb + el*0.001]
	
	fit_P = ['u']

	print (min_P, max_P)


#Define the star parameters to calculate the planet parameters
#Final values
mstar_mean  = float(mstar) #0.870
mstar_sigma = 0.1 #0.040

rstar_mean  = float(srad)
rstar_sigma = 0.1 

vsini_mean = 8.

tstar_mean  = float(teff)
tstar_sigma = 150.

print ("VALUE")
print (mstar_mean, rstar_mean, tstar_mean)

# ---------------
#Read input file
#execfile(inf_name)
#exec(open(inf_name).read())

#Prepare data
#execfile('pyaneti_LATTE/src/prepare_data.py')
exec(open('{}/pyaneti_LATTE/src/prepare_data.py'.format(syspath)).read())

#Create ouput directory
outdir = indir + '/' + tic + '/model_out'
if not os.path.exists(outdir):
  os.makedirs(outdir)

#Obtain smart priors based on iput data
if is_smart_priors:
  smart_priors()

print_init()

#-------------------------------------------------------------
#                   FITTING ROUTINES
#-------------------------------------------------------------

joint_fit()

#-------------------------------------------------------------
#             	PRINT AND PLOT ROUTINES
#-------------------------------------------------------------

#execfile('pyaneti_LATTE/src/output.py')

#exec(open('{}/pyaneti_LATTE/src/output.py'.format(syspath)).read())


#Print the values
#execfile('src/print_values.py')
exec(open('{}/pyaneti_LATTE/src/print_values.py'.format(syspath)).read())

#Create plots
#execfile('src/plot_data.py')
#execfile('src/plot_tr.py')
#execfile('src/plot_rv.py')
exec(open('{}/pyaneti_LATTE/src/plot_data.py'.format(syspath)).read())
exec(open('{}/pyaneti_LATTE/src/plot_tr.py'.format(syspath)).read())
exec(open('{}/pyaneti_LATTE/src/plot_rv.py'.format(syspath)).read())

if ( is_plot_chains ):
  plot_chains()
#  plot_postiter()


if ( is_corner_plot ):
  create_corner_plot()
else:
  if ( is_plot_posterior ):
    plot_posterior()

if ( is_plot_correlations ):
    plot_correlations()

 #PLOT TRANSIT
if ( total_tr_fit ):
#  plot_lightcurve_timeseries()
  create_folded_tr_plots()
#  if (nbands == 1):
#    plot_all_transits()
    #clean_transits(sigma_clean)
    #create_tango_input()

#PLOT RV CURVE
if ( total_rv_fit ):
  plot_rv_timeseries()
  plot_rv_phasefolded()

#-------------------------------------------------------------
#             	 END pyaneti.py FILE
#-------------------------------------------------------------
