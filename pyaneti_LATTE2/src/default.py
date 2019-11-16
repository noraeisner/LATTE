#-----------------------------------------------------------
#                       default.py
#  This file contains the defalut values to run pyaneti.
#  You can change all the values here. You can also control
#	 the values inside the input_file.py
#                 O. Barragan, March 2016
#-----------------------------------------------------------
#----------------------------------------------------------------
#                         CONSTANTS
#----------------------------------------------------------------
#Solar and planetary constants according to IAU 2015
# http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1605.09788
#Sun
S_radius_SI    = 6.957e8       	  #m
S_radius_cgs   = S_radius_SI*1e2  #cm
S_GM_SI        = 1.3271244e20     #m^3 s^{-1} 
S_GM_cgs       = S_GM_SI*1e6      #cm^3 s^{-1} 
G_SI           = 6.67408e-11      #m^3 kg^{-1} s^{-2}
G_cgs          = G_SI*1e3         #cm^3 g^{-1} s^{-2}
S_vol_SI       = 4./3.*np.pi*S_radius_SI**3  #m^3 
S_vol_cgs      = 4./3.*np.pi*S_radius_cgs**3 #cm^3
S_den_SI       = S_GM_SI/G_SI/S_vol_SI    #kg/m^3
S_den_cgs      = S_GM_cgs/G_cgs/S_vol_cgs #g/cm^3
S_Teff         = 5772.0 # K
#Earth
E_GM_SI        = 3.986004e14      #m^3 s^{-2}
E_radius_e_SI  = 6.3781e6         # ecuatorial radius [m]
E_radius_p_SI  = 6.3568e6         # polar radius [m]
#Jupiter
J_GM_SI        = 1.2668653e17     #m^3 s^{-1} 
J_radius_e_SI  = 7.1492e7         # ecuatorial radius [m]
J_radius_p_SI  = 6.6854e7         # polar radius [m]
#Other constants
AU_SI = 1.4960e11 # m
#from https://www.cfa.harvard.edu/~dfabricant/huchra/ay145/constants.html
m_hyd = 1.673534e-24 #g (hydrogen mass, cgs)
k_btz = 1.380658e-16 #erg/K (Boltzmann constant, cgs, 1 erg = g cm^2 / s^2)

#Default stellar parameters
mstar_mean = 1.0
rstar_mean = 1.0
tstar_mean = S_Teff
vsini_mean = 2.0
mstar_sigma = 1.e-1
rstar_sigma = 1.e-1
tstar_sigma = 100.
vsini_sigma = 0.1
#mag in J for the star
mag_j = 5.

#----------------------------------------------------------------
#                     MCMC CONTROLS
#----------------------------------------------------------------

#method can be 'mcmc' or 'plot'
#'mcmc' fits the data while 'plot' creates the plots using data of a previous run
method = 'mcmc'
#Number maximum of iterations
maxi = int(1e8)
#Thin factor to be used in the run
thin_factor = 10
#Number of iterations to be used in the run
#The code checks for convergence each niter \times thin_factor iterations
niter   = 100
#Number of chains or walkers to be used to exoplore parameter space
nchains = 100

#Indicate the number of planets to be fitted
nplanets = 1

#----------------------------------------------------------------
#                     FIT CONTROLS
#----------------------------------------------------------------

#Flags to control if we want to fit RV or transit models
#For multiplanet case, fit_rv and fit_tr are list in which each element
#refers to each fitted planet
fit_rv = [False]
fit_tr = [False]

#Parametrizations

#If True, the code sample for log10(P) instead of P
is_log_P    = False
#If True, the code sample for sqrt(e) sin w and sqrt(e) cos w instead of e and w
is_ew       = False
#If True, the code sample for impact parameter insted of inclination
is_b_factor = True
#If True, the code sample for the stellar density to the 1/3 power (rho^{1/3})
#instead of scaled semi major axis
#for multiplanet fits, it fits, the scaled semi-major axis of all fitted planets
#are calculated from the fitted rho^{1/3}
is_den_a    = False
#If True, the code sample for log10(k) instead of K
is_log_k    = False
#NOT WORK NOW! If True, the code sample for log10(v0) instead of v0
is_log_rv0  = False

#flat to control paramter priors
#For multiplanet fits, this variable has to have N elemens are N planets we are fitting
fit_t0  = ['f']
fit_P   = ['f']
fit_e   = ['f']
fit_w   = ['f']
fit_i   = ['f']
fit_a   = ['f']
fit_rp  = ['f']
fit_k   = ['f']
fit_v0  = ['f']
fit_q1  = ['f']
fit_q2  = ['f']

#if is_ew = True
fit_ew1 = ['f']
fit_ew2 = ['f']

#if is_b_factor = True
fit_b   = ['f']

is_linear_trend = 'f'
is_quadratic_trend = 'f'

#Default priors rages
#For multiplanet fits, this variable has to have N elemens are N planets we are fitting
min_t0  = [0.0]          #days
max_t0  = [1e6]         #days
min_P   = [0.1]       #days
max_P   = [1e4]      #days
min_e   = [1.e-10]      #zero
max_e   = [0.999]       #one
min_w   = [0.0]      #rad
max_w   = [2.0*np.pi]      #rad
min_ew1 = [-1.0]
max_ew1 = [1.0]
min_ew2 = [-1.0]
max_ew2 = [1.0]
#transit fit
min_i   = [0.0 ]    # 70 degrees
max_i   = [np.pi/2.0]  # 90 degrees
min_b   = [0.0]
max_b   = [1.0]
min_a   = [1.5]       # The planet is outside the star
max_a   = [1.e3]       # The planet is really far
min_q1  = [0.0]   #
max_q1  = [1.0]     #
min_q2  = [0.0]     #
max_q2  = [1.0]      #
min_rp  = [0.0]       # Earth size planet / sun
max_rp  = [0.99]      # a really big planet
#rv fit
min_k   = [1.e-6]      # m/s amplitudes
max_k   = [30.]     # a really big planet

jtr_prior_flag = []
jtr_prior_vals = []
jrv_prior_flag = []
jrv_prior_vals = []

#Flag to control single transit fits
is_single_transit = False

#TRANSIT FIT

#Number of points to be used to integrate the data
n_cad = 1
#time of integrations in days
t_cad = 30. / 60. / 24.
#Default lc_data is free, this means no re-binning of light curve model
lc_data = 'free'
#other options are:
#'kepler_lc', with default n_cad = 10 and t_cad = 30 minutes
#'kepler_sc', with default n_cad = 1 and t_cad = 1.5 minutes

#bin light curve each bin_lc time (units of days)
#this is useful for light curve fits with GPs
bin_lc = 0

#New parameters
nbands = 1
nldc   = 2
bands = ['']

#If True, it will set a Gaussian prior for the scaled semi-major axis based
#on the input stellar paramters
#BE AWARE THAT YOU HAVE TO BE CONFIDENT ABOUT YOUR STELLAR PARAMETERS TO USE IT!
#For multi-planet fits, a_from_kepler has to be N elements are planets are in the fit
a_from_kepler = [False]

#Default input columns for light curve file
#first column -> time
#second -> flux
#third ->  flux_error
columns_tr = [0,1,2]

#If True, it fits for a jitter term for the light curve data
is_jitter_tr = False

#If the input file only have two columns, time and flux, we can specify the
#error bar for the whole data set with this variable
my_tr_err = 0


#this flag allows to fit for a different radius for each bands
is_multi_radius = False
nradius = 1

#RADIAL VELOCITY FIT

#Spectrometer controls
#Label in the RV input file for each telescope
telescopes = []
#Label to appear in the plots for each telescope
telescopes_labels = ['']

#If True, it fits for a jitter term for each instrument
is_jitter_rv = False
#if we want to assign the jitter term fitting in a different way than the nominal
#one, we can specify it with the is_special_jitter flag
#if true, it will take a fifth column in the RV data file, in which the labels
#in that file will be taken as the jitter for a different instrument.
#jrvvec = [] has to contain the labels in the fifth column of the RV file
is_special_jitter = False
jrvvec = []


#----------------------------------------------------------------
#          Gaussian Process controls
#----------------------------------------------------------------
kernel_rv = 'None'
kernel_tr = 'None'


#----------------------------------------------------------------
#                     PLOT CONTROLS
#----------------------------------------------------------------

#If True it creates a posterior distribution plot
is_plot_posterior = True
#If True, it plots the prior disribution over the posterior distribution if is_plot_posterio == True
is_plot_prior = True
#how many colums do we want in the posterior plot?
n_columns_posterior = 3
#If True it creates a correlation plot
is_plot_correlations = False
#If True it creates a iterations vs chains plot
is_plot_chains = True
#If True it creates a corner plot and the posterior and correlation plots are not created
is_corner_plot = False
#We can select which parameters we want to plot in the correlation and posterior plots
#OSCAR: more description on how to use it is needed
plot_parameters = []
#If True, it will create the plots using seaborn library
is_seaborn_plot = False
#Default color palette for seaborn plots
#More details about palettes at https://seaborn.pydata.org/tutorial/color_palettes.html
seaborn_palette = 'deep'


#Figure size default, python takes inches by default
figure_size_x = 23./2.56 #23cm
figure_size_y = 23./2.56/1.618
#Font size label for the plots
font_size_label = 18

#TRANSIT PLOTS

#If True, it plots each individual transit with the fit model
#For multiplanet fits, this variable has to have N elemens are N planets we are fitting
is_plot_all_tr = [False]
#If True it plots percentiles in the phase pholded Transit plot
is_special_plot_tr = False
#If True, it overplots an unbinned model in the light curve model
plot_unbinned_model = False
#If True, it overplots an unbinned model in the light curve model
plot_binned_data = False
#If True, the transit plot will have each data point with an error bar
#If False, all points will be filled circles with the error bar size showed in the down right corner
plot_tr_errorbars = False
#If True, the stardard deviation of the residuals is shown in the transit plot
is_plot_std_tr = False

#If is_jitter_tr = True, the error bar size will be changed in the plot accordingly
resize_tr = True
#We can select the window size for the transit plot, it has to be given in units of days
#For multiplanet fits, this variable has to have N elemens are N planets we are fitting
span_tr  = [0.0]
#if True, we can select the window size for the tr plot given y_lim_min and y_lim_max
select_y_tr = False
#if select_y_tr = True, is the maximum limit for the tr plot
y_lim_max = 1.05
#if select_y_tr = True, is the minimum limit for the tr plot
y_lim_min = 0.95
#When creating timeseries plot of light curve data, the code performs a sigma clipping
#algorithm to remove outliers, default sigma is 10, but it can be changed
sigma_clean = 10
#Label to appear in the timeseries plot of the light curve
tr_xlabel = "BJD - 2450000 (days)"

tr_colors = ['r','b','g','k','y','#ffbf00','#ff1493']

#RV PLOTS
#If True it plots percentiles in the phase pholded RV plot
is_special_plot_rv = False

#RVlabels and residuals
rv_labels = ['RV (m/s)']
rv_res    = ['Residuals (m/s)']

#If is_jitter_rv = True, the error bar size will be changed in the plot accordingly
#The jitter term is displayed as a gray extension to the nominal error bars
resize_rv = True

#if True, we can select the window size for the tr plot given rv_lim_min and rv_lim_max
select_y_rv = False
#if select_y_rv True, is the maximum limit for the rv plot
rv_lim_max = -100. #m/s
#if select_y_rv True, is the minimum limit for the rv plot
rv_lim_min =  100. #m/s
#Label to appear in the RV timeseries plot
rv_xlabel = "BJD - 2450000 (days)"
#if True, the RV plots contain a legend indicating the color and symbol of each instrument
is_rv_legend = True
#Default markers for the different instruments of the RV plots
mark = ['o', 'D', 's', 'p', 'h', '8', '^', '<', '*', 'v','>','.', 'H', 'd','+']
mark_tr = ['o', 'D', 's', 'p', 'h', '8', '^', '<', '*', 'v','>','.', 'H', 'd','+']
#Default colors for the different instruments of the RV plots
rv_colors = ['#007fff','r','g','y','m','c']

#marker size for the data points of the RV plot
rv_markersize  = 8
#Type of marker in the RV plot
rv_fillstyle = 'full'


#----------------------------------------------------------------
#                     OTHER CONTROLS
#----------------------------------------------------------------

#unit_mass controls the units in which the planet parameters are displayed
#'solar' prints it in solar units
#'earth' prints it in Earth units
#'jupiter' prints it in Jupiter units
unit_mass = 'solar'

#if True, it will print the mode and 99% limits of the credible interval distribution
is_print_mode = False

#If True it scales the error bars for the case when chi^2/dof < 1
scale_error_bars = False

#If True, it will check for clustering chains in order to remove lost chain in
#local minima
is_clustering = True
#The clustering_delta parameter checks for a threshold from the median value of the chains
#all chains with chi2/dof above median + 3*clustering_delta are rejected
#BE WARE THAT a small value may reject converged chains in the local minima
clustering_delta = 0.5

#If True, the code creates authomatic priors for the RV offests
#And also checks for maximum values for Rp/R*
#BE AWARE THAT If False, the code may not work!
is_smart_priors = True

#If True, it creates a LaTeX file with all fitted and derived parameters
latex_values = True

#Default output directory is outpy, this can be changed with outdir variable
outdir = 'outpy/'

#For a joint fit, in the case in which RV and Light curve data have not the same
#units of time, the variable textra adds a constant time to the light curve time
#series in order to match both time sets.
#Example:
#if RV data is given in BJD - 2,450,000 and
#light curve data is given in BJD - 2,454,833
#textra = 4833 will add 4833. days to the light curve time series in such a way that
#both data sets will be now in BJD - 2,450,000
textra = 0.0

#The default planet labels
#If you need to add more labels, you have made a great discovery! ;-)
plabels = ['b','c','d','e','f','g','h','i','k','l']

#Select how the values are taken to create plots
#'median' -> it uses the median
#'mode' -> it uses the mode
get_value = 'median'


# ---------
# THESE WERE PREVIOUSLY IN THE INPUT_FIT.PY FILE
# ---------

#Input file for test problem for pyaneti
#Created by Barragan O.
#This file contains the basic parameters to change for simple pyaneti run
#There are more flags which can control the input/output of the code
#They can be found inside src/default.py, all the flags inside this file
#can be changed inside the input file.

fname_tr = ['{}_data_pyaneti.dat'.format(tic)]

#bands = ['LC','SC']

#t_cad = [30./60./24,2./60./24.]
#n_cad = [10,1]


#fname_tr = ['sec2_detrended.dat']
#fname_tr = ['sec5_detrended.dat']
#fname_tr = ['sec8_detrended.dat']


#MCMC controls
#the thin factor for the chains
thin_factor = 1
#The number of iterations to be taken into account
#The TOTAL number of iterations for the burn-in phase is thin_factor*niter
niter       = 500
#Number of independent Markov chains for the ensemble sampler
nchains     = 100

span_tr = [2]

#Choose the method that we want to use
# mcmc -> runs the mcmc fit program
# plot -> this option create the plots only if a previus run was done
method = 'mcmc'
#method = 'plot'


plot_binned_data = False

is_seaborn_plot = True

is_print_mode = True

#Define the star parameters to calculate the planet parameters
#Final values
mstar_mean  = 1.320 #0.870
mstar_sigma = 0.063 #0.040
rstar_mean  = 1.942
rstar_sigma = 0.1 
vsini_mean = 8.
#
tstar_mean  = 5907.
tstar_sigma = 150.
#is_clustering = False

#What units do you prefer for your planet parameters?
# earth, jupiter or solar
unit_mass = 'earth'

#If we want posterior, correlation and/or chain plots these options have to be set True
is_plot_posterior    = True
is_plot_correlations = False
is_plot_chains       = False

#Are we setting gaussian priors on the semi-major axis based on the stellar parameters?
a_from_kepler = [False]*nplanets

#For a pure RV fit, fit_tr has to be False
#For a pure TR fit, fit_rv has to be False
#For multi-planet fits fit_rv and fit_tr have the form [True,True,False,...]
#one element for each planet.
fit_rv = [False]*nplanets
fit_tr = [True]*nplanets

#is_ew controls the parametrization sqrt(e)sin w and sqrt(e) cos w
#if True we fit for the parametrization parameters
#if False we fit for e and w
#Default is True
is_ew = False

#Prior section
# f -> fixed value
# u -> Uniform priors
# g -> Gaussian priors
fit_t0 = ['u','u','u']   #We fit for t0 with uniform priors
fit_P  = ['f','u','u']   #We fit for P with uniform priors
fit_e  = ['f','f','f']   #We fix e, it works only if is_ew = False
fit_w  = ['f','f','f']   #We fix w, it works only if is_ew = False
fit_ew1 = ['u']*nplanets
fit_ew2 = ['u']*nplanets
fit_b  = ['u','u','u']   #We fix the impact factor
fit_a  = ['u','u','u']   #We fit a with gaussian priors (given by the stellar parameters)
fit_rp = ['u','u','u']   #We fit rp with uniform priors
fit_k  = ['f','f','f']   #We fit k with uniform priors

fit_v0 = 'f'     #We fit systemc velicities with uniform priors
fit_q1 = ['u']*2     #We fit q1 with gaussian priors
fit_q2 = ['u']*2#We fit q2 with gaussian priors

min_q1 = [0.]*2
max_q1 = [1.]*2
min_q2 = list(min_q1)
max_q2 = list(max_q1)

#Prior ranges for a parameter A
#if 'f' is selected for the parameter A, A is fixed to the one given by min_A
#if 'u' is selected for the parameter A, sets uniform priors between min_A and max_A
#if 'g' is selected for the parameter A, sets gaussian priors with mean min_A and standard deviation max_A


min_a   = [1.2,31.68,5]
max_a   = [1000,0.70,34]

val = 0.99
min_ew1 = [-val]*nplanets
max_ew1 = [ val]*nplanets
min_ew2 = [-val]*nplanets
max_ew2 = [ val]*nplanets

min_b   = [0.0]*nplanets
max_b   = [1.0]*nplanets

min_k   = [0.000]*nplanets
max_k   = [0.010]*nplanets

min_rp  = [0.0]*nplanets
max_rp  = [0.3]*nplanets


min_e = [0.0]*nplanets
max_e = [1.0]*nplanets
min_w = [0.0]*nplanets
max_w = [2*np.pi]*nplanets

if fit_e[0] == 'b':
  min_e   = [0.867]
  max_e   = [3.03]

is_multi_radius = False

is_jitter_tr=True

is_plot_all_tr = [False]*nplanets

is_den_a = False

if is_den_a:
  min_a   = [3.]*nplanets        
  max_a   = [20.]*nplanets    






