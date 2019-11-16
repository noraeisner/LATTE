import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../') #Path to lard pyaneti module
import pyaneti as pti

#Call the default values file
execfile('../src/default.py')

#Script to genereate fake RV and Transit data
#for a given planet

#Define here the properties of your planet
#Period
P =  8.872810
#Transit epoch
T0 = 7143.3560292 + 100*P
#Eccentricity
e = 0
#Periastron
w = np.pi/2
#Mass (Earh masses)
Mp = 1.0
#Limb darkening coefficients
#      0.45816881      0.24219662
u1 = 0.458
u2 = 0.242
#Scaled planet radius
pz = (E_radius_e_SI/S_radius_SI)
#Scaled Semi-major axis
a = (149.60e9/S_radius_SI)
#inclination
ii = np.pi/2
#Systemic velocity 
rv0 = -76.6149275

#------------------------------------------#
#            The magic begins              #
#------------------------------------------#

#Calculate the semi-amplitude parameter
k = (2*np.pi*G_SI/(P*24*3600))**(1./3) \
    * (Mp*E_GM_SI/G_SI)/(S_GM_SI/G_SI)**(2./3)
k = 4.1499560
k_km = k * 1e-3 #transform to Km


xderv = np.arange(T0,T0 + P,1.0)
err_rv = 3.5e-3
rv = pti.rv_curve_mp(xderv,rv0,T0,k_km,P,e,w,0.0,0.0)
rv2 = np.random.normal(rv,scale=err_rv)
rvdata = open('c59627_data.dat','w')
for i in range(0,len(xderv)):
    rvdata.write('%4.10f%  4.7f%  4.7f %s \n'%(xderv[i], rv2[i], err_rv, 'HN'))
rvdata.close()

plt.plot(xderv,rv2,'o',xderv,rv)
plt.show()
