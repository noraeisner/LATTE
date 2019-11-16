import numpy as np
import pyaneti as pti
import matplotlib.pyplot as plt

execfile('src/default.py')

#Script to genereate fake 2 planets RV data

P = [176.251,264.377]
e = [0.0,0.0]
w = [np.pi/10.,3.345*np.pi/4.]
#Let us add a random julian date
#+245
T0 = [7303.023611,7567.400694]
#k
k = [None]*2
for i in [0,1]:
 k[i] = (2*np.pi*G_SI/(P[i]*24*3600))**(1./3) * (E_GM_SI/G_SI)/(S_GM_SI/G_SI)**(2./3) * 1e-3
#k
k_km = [k[0] * 317.8, k[1] * 95.16 ]

rv0 = 35.67234

#--------------------------------------------

#xde = [None]*24*2*365
xde = [None]*24*2*365*3

for i in range(0,len(xde)):
    xde[i] = i%48./48. + int(i/48.)
    xde[i] = xde[i] + T0[1] - 28.4

xderv = xde[0::48*30]
rv = pti.rv_curve_mp(xderv,rv0,T0,k_km,P,e,w)
error_median = 5e-3
error_sigma = 1e-3
err_rv = np.random.normal(error_median,scale=error_sigma,size=len(rv))
rv2 = [None]*len(rv)
for i in range(0,len(xderv)):
    rv2[i] = np.random.normal(rv[i],scale=err_rv[i])
    print xderv[i], rv2[i], err_rv[i], 'S'

#plt.plot(xderv,rv,'ro')
#plt.errorbar(xderv,rv2,err_rv)
