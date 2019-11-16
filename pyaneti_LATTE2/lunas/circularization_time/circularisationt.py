import numpy as np

#Define parameters
ms = 1.15 #M_Sun
rs = 1.24 #R_Sun
qs = 10**(6.5)
qs = 10**(5)
mp = 21.2/317.8
rp = 3.85/11.2
qp = 10**(5.5)
qp = 35000
a  = 0.0305 #AU

#--------------------------------------------------------------

#Define constants
#Solar and planetary constants according to IAU 2015
# http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1605.09788
#Sun
S_radius_SI    = 6.957e8          #m
S_GM_SI        = 1.3271244e20     #m^3 s^{-1}
G_SI           = 6.67408e-11      #m^3 kg^{-1} s^{-2}
#Jupiter
J_GM_SI        = 1.2668653e17     #m^3 s^{-1}
J_radius_e_SI  = 7.1492e7         # ecuatorial radius [m]
J_radius_p_SI  = 6.6854e7         # polar radius [m]
#Other constants
AU_SI = 1.4960e11 # m
yr    = 60*60*24*365.25 # s
Gyr   = yr*1e9 # s
Myr   = yr*1e6 # s

#Convert all the parameters to SI
ms = ms*S_GM_SI/G_SI
rs = rs*S_radius_SI
mp = mp*J_GM_SI/G_SI
rp = rp*J_radius_e_SI
a  = a*AU_SI

#Eq. (1) Jackson et al., 2008, ApJ, 725, 1995
ededt = 63./4. * np.sqrt(G_SI * ms**3) * rp**5/(qp*mp)
ededt = ededt + 171./16. * np.sqrt(G_SI/ms) *rs**5*mp/qs
ededt = - (ededt) * a**(-13./2.)

tcirc = - 1./ededt

print 'Circularisation time = ', tcirc/Gyr, ' Gyr'
print 'Circularisation time = ', tcirc/Myr, ' Myr'


