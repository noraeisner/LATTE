import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../') #Path to lard pyaneti module
import pyaneti as pti
import seaborn as sns
sns.set_color_codes()
#sns.set(style='ticks')


#Call the default values file
execfile('../../src/default.py')

#period days
#M solar masses
#R solar radii
def get_scaled_a(P,M,R):
  a = S_GM_SI*M*(P*3600.*24.)**2
  a = a / (4. * np.pi**2)
  a = a**(1./3.)
  #a in meters
  a = a / (R * S_radius_SI)
  print a
  return a

def get_K(P,i,a,e,mp,M,R):
    n = 2*np.pi/(P*3600.*24.) 
    K = ( mp*E_GM_SI + M*S_GM_SI ) * ( np.sqrt( 1.0 - e*e ) ) 
    K = mp*E_GM_SI * n * a * R * S_radius_SI * np.sin(i) / K
    print K, 'm/s'
    return K
  

#Script to genereate syntetic data for a multiplanet  model

#stellar parameters
mstar = 1.0
rstar = 1.0
#Limb darkening coefficients
u1 = 0.54321
u2 = 0.12345

#Define here the properties of your planet
#Period
P1 =  5.6
P =[P1]
#Transit epoch
T01 = 10.0 
T0 = [T01]
#Eccentricity
e1 = 0
e = [e1]
#Periastron
w1 = np.pi/2
w = [w1]
#Scaled planet radius
pz1 = (3.5*E_radius_e_SI/S_radius_SI)
rp = [pz1]
#masses
mp = [12.3] #Earth masses
#Scaled Semi-major axis
#calculated from period and stellar parameters
#inclination
ii1 = 87.*np.pi/180.
ii = [ii1]

ndays = 30
tr_ppday = 24*20.
rv_points = 100
err_lc = 5.e-5
rv_err = 1. #m/s
rv_offset = 0 #km/s
rv_inst = 'A'
colors = ['g','y','r']
lsv = ['dashdot','dashed','solid']

#------------------------------------------#
#            The magic begins              #
#------------------------------------------#

#Calculate scaled semi-major axis
npl = len(P)
a = [None]*npl
for o in range(0,npl):
  a[o] = get_scaled_a(P[o],mstar,rstar)


#CREATE TR DATA
#Create the vector for the time-stamps for the TR data
xde = [None]*(int(tr_ppday)*ndays)

#Create the time-stamps
for i in range(0,len(xde)):
    xde[i] = i%tr_ppday/tr_ppday + int(i/tr_ppday) 
    #xde[i] = xde[i]# + T0[0] - P1/4

flux_tot = 0.0
#Create the flux vector
for o in range(0,npl):
  pars = [T0[o],P[o],e[o],w[o],ii[o],a[o]]
  z    = pti.find_z(xde,pars,[0,0,0,0])
  flux0, sdf = pti.occultquad(z,u1,u2,rp[o])
  flux_tot = flux_tot + flux0

#Calculate the whole light curve
flux_tot = flux_tot  - npl + 1

#Add the gaussian noise
flux_total_sim = np.random.normal(flux_tot,scale=err_lc)
fos = 3.5

plt.figure(1,figsize=(fos,fos/1.618))
plt.plot(xde,flux_tot,'-',linewidth=5)
plt.plot(xde,flux_total_sim,'ro',markersize=4)
plt.xlabel('Time [arbitrary units]',fontsize=fos*3)
plt.ylabel('Relative flux',fontsize=fos*3)
#
plt.tick_params(labelsize=fos)
plt.ticklabel_format(useOffset=False, axis='y')
plt.tick_params(labelsize=fos*3)
plt.tick_params(labelsize=fos*3)
#
plt.xlim(T0[0]-3./24,T0[0]+3/24.)
plt.ylim(0.9988,1.0002)
plt.xticks([T0[0]],' T$_0$')
#plt.legend(loc=4, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize=fos*0.7)
plt.savefig('transit.pdf',bbox_inches='tight')
plt.savefig('transit.png',bbox_inches='tight',dpi=500)
plt.show()

trd = open('tr_data.dat','w')
for i in range(0,len(xde)):
  trd.write('%8.8f %8.8f %8.8f \n'%(xde[i], flux_total_sim[i], err_lc))
trd.close()

#CREATE RV DATA

#Calcualte the semi-amplitudes
k = [None]*npl
for o in range(0,npl):
  k[o] = get_K(P[o],ii[o],a[o],e[o],mp[o],mstar,rstar)
  k[o] = k[o]*1e-3 #Transform to km/s
 
#Create RV timestamps
#min_rv = np.min(xde) 
#max_rv = np.max(xde) 
min_rv = T0[0]
max_rv = T0[0] + P[0]
#rv_time = np.random.uniform(min_rv,max_rv,rv_points)
rv_time = np.arange(min_rv,max_rv,(max_rv - min_rv)/(rv_points-1))
#Create the data
rv = pti.rv_curve_mp(rv_time,rv_offset,T0,k,P,e,w,0.0,0.0)
#Add the Gaussian noise
rv_noise = np.random.normal(rv,rv_err*1e-3)


#Plot the data
plt.figure(1,figsize=(fos,fos/1.618))
plt.xlim(T0[0],T0[0]+P[0])
plt.plot(rv_time,rv*1e3,'-',linewidth=5)
plt.plot(rv_time,rv_noise*1e3,'ro',markersize=4)
#
plt.tick_params(labelsize=fos)
plt.ticklabel_format(useOffset=False, axis='y')
plt.xticks([T0[0],T0[0]+P[0]/2,T0[0]+P[0]],('0','0.5','1'))
plt.tick_params(labelsize=fos*3)
#
plt.xlabel('Phase',fontsize=fos*3)
plt.ylabel('RV [arbitrary units]',fontsize=fos*3)
plt.savefig('rv.pdf',bbox_inches='tight')
plt.savefig('rv.png',bbox_inches='tight',dpi=500)
plt.show()

