import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.colors import LogNorm
sys.path.append('../../') #Path to lard pyaneti module
import pyaneti as pti
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit
#sns.set_color_codes()
#sns.set(style='ticks')

#Call the default values file
execfile('../../src/default.py')


#Here is the input
#Read the time-stamps and errors, it can be modified to read any file with this data
x,y,errs = np.loadtxt('RV.C10_5255.dat',unpack=True, usecols=(0,1,2))

#Parameters definitions
T0 = 7588.2850480 
P  = 6.5692094
e  = 0.0
w  = np.pi/2
K  = 98.9875886*1.e-3
g  = 1.2455757 #Gamma
l  = 0.0 #Linear trend
q  = 0.0 #Quadratic trend

fsize = 20

ew_par = True

#Number of times that the data is modeled for the Monte Carlo experiment
niterations = 10000

#-----------------------------------------------------------------
#                          The magic begins here
#-----------------------------------------------------------------

#The idea is:
#Create a for cycle which creates synthetic data for each set of timestamps
#add gaussian noise
#fit the data using LM method with the pyaneti orbits
#Create a sample of parameters to see if the eccentricity is significant

#Create the model to perturb
rv_good = pti.rv_curve_mp(x,g,T0,K,P,e,w,l,q)

#It assumes no gamma, no linear or quadratic trends
def frv(x,K,ee,ww,g):
  global T0, P, ew_par #Assumes constant T0 and P
  if ( ew_par ):
    e  = ee*ee + ww*ww
    w  = np.arctan2(ee,ww)
  else:
    e  = ee
    w  = ww
  #Create function which allows parametrization of e and w
  #Here add flags for the parametrization
  rv = pti.rv_curve_mp(x,g,T0,K,P,e,w,0.0,0.0)
  
  return rv

if ( ew_par ):
  myp0 = [K,e*np.sin(w),e*np.cos(w),g]
  mybounds = ([0.,-0.5,-0.5,g-0.1],[1.0,0.5,0.5,g+0.1])
else:
  myp0 = [K,e,w]
  mybounds = ([0.,0.,0.],[1.0,1.0,2*np.pi])

#Create the for cycle
k_vec = [None]*niterations
ee_vec = [None]*niterations
ww_vec = [None]*niterations
for i in range(0,niterations):

  #This creates gaussian noise to the original model 
  rv_perturbed = np.random.normal(rv_good,errs)

  #Let us make the fit using curve fit
  fpars, dummy = curve_fit(frv,x,rv_perturbed,p0=myp0,bounds=mybounds,sigma=errs)
  k_vec[i] = fpars[0]
  ee_vec[i] = fpars[1]
  ww_vec[i] = fpars[2]

#Until This point we have all our parameters 

mydata = [None]*niterations 
for o in range(0,niterations):
  mydata[o] = [ee_vec[o],ww_vec[o]]

plt.figure(figsize=(fsize,fsize))
df = pd.DataFrame(mydata, columns=["x", "y"])
g = sns.jointplot(x="x",y="y",data=df,kind="kde",space=0,stat_func=None)
g.set_axis_labels("$\sqrt{e} \, \sin \, \omega_\star$", "$\sqrt{e} \, \cos \, \omega_\star$",fontsize=fsize)
g.savefig('ewprob.pdf')
plt.show()


mye = [0.0]*niterations
for o in range(0,niterations):
  mye[o] = ee_vec[o]**2 + ww_vec[o]**2

plt.hist(mye,bins=50,normed=True) #This should be the eccentricity
plt.xlabel('$e$',fontsize=fsize)
plt.ylabel('Probability',fontsize=fsize)
plt.tick_params(labelsize=fsize,direction='in')
plt.savefig('eprob.pdf')
plt.show()

#plt.figure(figsize=(fsize,fsize))
#plt.xlabel('$\sqrt{e} \, \sin \, \omega_\star$',fontsize=fsize)
#plt.ylabel('$\sqrt{e} \, \cos \, \omega_\star$',fontsize=fsize)
#plt.tick_params(labelsize=fsize,direction='in')
#plt.hist2d(ee_vec,ww_vec,bins=50,normed=True,norm=LogNorm())
#cbar = plt.colorbar(orientation='horizontal')
#cbar.ax.set_xlabel('$\log P$',fontsize=fsize)
#plt.savefig('ewprobnorm.pdf')
#plt.show()
