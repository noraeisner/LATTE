import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../') #Path to lard pyaneti module
import pyaneti as pti

#Call the default values file
execfile('../src/default.py')

#period days
#M solar masses
#R solar radii
def get_scaled_a(P,M,R):
  a = S_GM_SI*M*(P*3600.*24.)**2
  a = a / (4. * np.pi**2)
  a = a**(1./3.)
  #a in meters
  a = a / (R * S_radius_SI)
  return a
  

#Script to genereate fake Transit data for a two planet model

#stellar parameters
mstar = 1.0
rstar = 1.0

#Define here the properties of your planet
#Period
P0 =  0.66
P1 =  5.6
P2 =  12.1
#Transit epoch
T00 = 0.4 
T01 = T00 + P0*3 
T02 = T00 + P0*3 - 0.05
#Eccentricity
e0 = 0
e1 = 0
e2 = 0.4
#Periastron
w0 = np.pi/2
w1 = np.pi/2
w2 = np.pi/2 + 0.67
#Limb darkening coefficients
u1 = 0.54321
u2 = 0.12345
#Scaled planet radius
pz0 = (1.5*E_radius_e_SI/S_radius_SI)
pz1 = (3.5*E_radius_e_SI/S_radius_SI)
pz2 = (5.5*E_radius_e_SI/S_radius_SI)
#Scaled Semi-major axis
a0 = get_scaled_a(P0,mstar,rstar)
a1 = get_scaled_a(P1,mstar,rstar)
a2 = get_scaled_a(P2,mstar,rstar)
#inclination
ii0 = 89.*np.pi/180.
ii1 = 87.*np.pi/180.
ii2 = 88.*np.pi/180.

#------------------------------------------#
#            The magic begins              #
#------------------------------------------#

#each 30 min, for 90 days
vd =24*6.0

xde = [None]*int(vd)*28

for i in range(0,len(xde)):
    xde[i] = i%vd/vd + int(i/vd) 
    xde[i] = xde[i] + T01 - P1/4

pars0 = [T00,P0,e0,w0,ii0,a0]
pars1 = [T01,P1,e1,w1,ii1,a1]
pars2 = [T02,P2,e2,w2,ii2,a2]

z0 = pti.find_z(xde,pars0,[0,0,0,0])
z1 = pti.find_z(xde,pars1,[0,0,0,0])
z2 = pti.find_z(xde,pars2,[0,0,0,0])

flux0, sdf = pti.occultquad(z0,u1,u2,pz0)
flux1, sdf = pti.occultquad(z1,u1,u2,pz1)
flux2, sdf = pti.occultquad(z2,u1,u2,pz2)

flux_total = flux1 + flux2 + flux0  - 2.0

err_lc = 5.e-5

flux_total_sim = np.random.normal(flux_total,scale=err_lc)

plt.plot(xde,flux_total_sim)
plt.show()

for i in range(0,len(xde)):
  print xde[i], flux_total_sim[i], err_lc
