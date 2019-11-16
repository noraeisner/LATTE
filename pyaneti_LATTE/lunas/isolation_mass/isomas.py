#!/usr/bin/python
import numpy as np

mstar = 1.25 #solar masses
saxis = 0.424 #AU
mplanet = 42. #Earth masses
#This analysis is based on
#Schlichting, 2014, ApJ, 795, L15

#Calculate isolation mass by knowing
#stellar mass and planet distance
miso = (3. * mstar * 1.99e33 )**(0.5) # gr**(1/2)
miso = ( 10. * np.pi * 7.0 * (saxis)**(-1.5) * \
       (saxis * 1.496e13)**2 )**(1.5) / miso # gr
miso = miso / 5.97e27 #Earth masses

#Let us calculated the required Sigma to form the planet
#in situ
snew = (mplanet / miso )**(2.0/3.0)

#Let us get the Toomre factor with this new surface density
#eq. (6)
#Assuming an isothermal disk, and a gas-to-dust ratio 200
qnew = ( snew * 7.0 * (saxis)**(-1.5) / 1e4 ) #sigma new disk
qnew = 4 * (saxis/0.1)**(-1.5) / qnew

print ''
print 'Isolation mass = ', miso, ' Earth masses'
print ''
print 'Your planet would need a CSD ', snew, 'more dense'
print 'to for the planet in situ'
print ''
print 'This leds to a Toomre factor of Q= ', qnew
print ''

