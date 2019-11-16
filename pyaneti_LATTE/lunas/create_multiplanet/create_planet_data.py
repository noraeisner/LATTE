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
#Transit epoch
T0 = 2448285.090278 
#Period
P =  365.256
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
rv0 = 22.071987

#------------------------------------------#
#            The magic begins              #
#------------------------------------------#

#Calculate the semi-amplitude parameter
k = (2*np.pi*G_SI/(P*24*3600))**(1./3) \
    * (Mp*E_GM_SI/G_SI)/(S_GM_SI/G_SI)**(2./3)
k_km = k * 1e-3 #transform to Km

xde = [None]*24*2*365*3

for i in range(0,len(xde)):
    xde[i] = i%48./48. + int(i/48.) 
    xde[i] = xde[i] + T0 - 28.4

pars = [T0,P,e,w,ii,a]

z = pti.find_z(xde,pars,[0,0,0,0])

flux, sdf = pti.occultquad(z,u1,u2,pz)

xdetr = xde[28*48:29*48]
xdetr.append(xde[28*48+48*365:29*48+48*365])
xdetr.append(xde[28*48+48*365*2:29*48+48*365*2])

err_lc = 1.0e-5

flux2 = np.random.normal(flux,scale=err_lc)
trans = open('lc_data.dat','w')
for i in range(0,len(xde)):
    trans.write('%4.7f  %4.7f  %4.7f \n'%(xde[i], flux2[i], err_lc))
trans.close()

xderv = xde[0::48*31]
err_rv = 1e-5
rv = pti.rv_curve_mp(xderv,rv0,T0,k_km,P,e,w,0.0,0.0)
rv2 = np.random.normal(rv,scale=5e-6)
rvdata = open('rv_data.dat','w')
for i in range(0,len(xderv)):
    rvdata.write('%4.7f%  4.7f%  4.7f %s \n'%(xderv[i], rv2[i], err_rv, 'S'))
rvdata.close()
