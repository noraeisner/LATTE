#Torres relations

import numpy as np

#Gaussian error propagation
logg = [4.15,0.05]
teff = [5635.,50.]
feh  = [-0.5,0.1]

logg = [4.10,0.05]
teff = [5585.,50.]
feh  = [-0.6,0.05]

#logg = [4.20,0.05]
#teff = [5685.,50.]
#feh  = [-0.40,0.05]

ntot = int(1e6)

#------------------------------------------------------------

g_logg = np.array(np.random.normal(loc=logg[0],scale=logg[1],size=ntot))
g_teff = np.array(np.random.normal(loc=teff[0],scale=teff[1],size=ntot))
g_feh  = np.array(np.random.normal(loc= feh[0],scale= feh[1],size=ntot))

#------------------------------------------------------------

#Create Torres variables
#a
a1 = np.random.normal(1.5689,0.058,ntot)
a2 = np.random.normal(1.3787,0.029,ntot)
a3 = np.random.normal(0.4243,0.029,ntot)
a4 = np.random.normal(1.139,0.024,ntot)
a5 = np.random.normal(-0.1425,0.011,ntot)
a6 = np.random.normal(0.01969,0.00019,ntot)
a7 = np.random.normal(0.1010,0.014,ntot)
#b
b1 = np.random.normal(2.4427,0.038,ntot)
b2 = np.random.normal(0.6679,0.016,ntot)
b3 = np.random.normal(0.1771,0.027,ntot)
b4 = np.random.normal(0.705,0.13,ntot)
b5 = np.random.normal(-0.21415,0.0075,ntot)
b6 = np.random.normal(0.02306,0.0013,ntot)
b7 = np.random.normal(0.04173,0.0082,ntot)

X = np.log10(g_teff) - 4.1

logm = a1 + a2*X + a3*X**2 + a4*X**3 + a5*g_logg**2 + a6*g_logg**3 + a7*g_feh
logr = b1 + b2*X + b3*X**2 + b4*X**3 + b5*g_logg**2 + b6*g_logg**3 + b7*g_feh

minm, medm, maxm = np.percentile(logm,[16,50,84])
minr, medr, maxr = np.percentile(logr,[16,50,84])
print 'Mass   = ',10**medm,'-',10**medm-10**minm,'+',10**maxm-10**medm
print 'Radius = ',10**medr,'-',10**medr-10**minr,'+',10**maxr-10**medr


