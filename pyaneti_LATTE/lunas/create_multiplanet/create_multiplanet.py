import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../') #Path to load pyaneti module
import pyaneti as pti
import seaborn as sns
sns.set_color_codes()
sns.set(style='ticks')


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
P0 =  0.66
P1 =  5.6
P2 =  12.1
P3 =  20.3
P =[P0,P1,P2]
#Transit epoch
T00 = 0.4 
T01 = T00 + P0*3 
T02 = T00 + P0*3 - 0.05
T03 = T00 + P0
T0 = [T00,T01,T02,T03]
#Eccentricity
e0 = 0
e1 = 0
e2 = 0.3
e3 = 0.2
e = [e0,e1,e2,e3]
#Periastron
w0 = np.pi/2
w1 = np.pi/2
w2 = np.pi/2 + 0.67
w3 = np.pi/2 + 1.7
w = [w0,w1,w2,w3]
#Scaled planet radius
pz0 = (1.5*E_radius_e_SI/S_radius_SI)
pz1 = (3.5*E_radius_e_SI/S_radius_SI)
pz2 = (5.5*E_radius_e_SI/S_radius_SI)
pz3 = (6.0*E_radius_e_SI/S_radius_SI)
rp = [pz0,pz1,pz2,pz3]
#masses
mp = [3.1,12.3,26.0,30.0] #Earth masses
#Scaled Semi-major axis
#calculated from period and stellar parameters
#inclination
ii0 = 89.*np.pi/180.
ii1 = 87.*np.pi/180.
ii2 = 88.*np.pi/180.
ii3 = 84.*np.pi/180.
ii = [ii0,ii1,ii2,ii3]

ndays = 30
tr_ppday = 24*12.
rv_points = 100
err_lc = 5.e-5
rv_err = 5. #m/s
rv_offset = 10 #km/s
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
  z    = pti.find_z(xde,pars)
  flux0, sdf = pti.occultquad(z,u1,u2,rp[o])
  flux_tot = flux_tot + flux0

#Calculate the whole light curve
flux_tot = flux_tot  - npl + 1

#Add the gaussian noise
flux_total_sim = np.random.normal(flux_tot,scale=err_lc)
fos = 25
plt.figure(1,figsize=(fos,fos/3.))
plt.plot(xde,flux_total_sim,'.')
plt.xlabel('Days',fontsize=fos)
plt.ylabel('Relative flux',fontsize=fos)
plt.tick_params(labelsize=fos)
for o in range(0,npl-1):
  myrange = np.arange(T0[o],T0[o]+ndays,P[o])
  plt.vlines(myrange,1.0003,1.0008,colors[o],lsv[o],label='Planet with period '+str(P[o])+' days')
plt.xlim(np.min(xde),np.max(xde))
plt.legend(loc=4, ncol=1,scatterpoints=1,numpoints=1,frameon=False,fontsize=fos*0.7)
plt.savefig('light_curve.pdf',bbox_inches='tight')
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
min_rv = np.min(xde) 
max_rv = np.max(xde) 
rv_time = np.random.uniform(min_rv,max_rv,rv_points)
#Create the data
rv = pti.rv_curve_mp(rv_time,rv_offset,T0,k,P,e,w,0.0,0.0)
#Add the Gaussian noise
rv_noise = np.random.normal(rv,rv_err*1e-3)
#Plot the data
plt.plot(rv_time,rv,'bo',rv_time,rv_noise,'rs')
plt.show()

rv_err = 3.

rvd = open('rv_data.dat','w')
for i in range(0,len(rv_time)):
  rvd.write('%8.8f %8.8f %8.8f %s \n'%(rv_time[i], rv_noise[i], rv_err*1e-3,rv_inst))
rvd.close()
