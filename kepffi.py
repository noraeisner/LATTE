
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import numpy as numpy
import pandas as pandas
from lightkurve import search_tesscut
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button
from sklearn.decomposition import PCA


searchtic = 'TIC 229940491'  #349972099
kepid = 229940491
sec = 14

#----- functions

search_result = lk.search_tesscut(searchtic, sector=sec)
tpf = search_result.download(cutout_size=11)

def norm(a):
    '''
    function to normalise the data - used in download_data_tpf
    '''
    m = np.median(a)
    x = a-m
    s = np.std(x)
    x = x/s
    return x, m, s

def unnorm(x,m,s):
    '''
    function to un-normalise the data - used in download_data_tpf
    '''
    y = x * s
    a = y + m
    return a

# ----- perform PCA ------

X1 = tpf.flux

s = X1.shape
X1 = X1.reshape(s[0],s[1]*s[2])

lkeep = np.isfinite(X1.sum(axis=1)) * (X1.sum(axis=1)>0)
X1 = X1[lkeep,:]

X2 = np.zeros_like(X1)
M = np.zeros(X1.shape[1])
S = np.zeros(X1.shape[1])

for n in range(X1.shape[1]):
    a = X1[:,n]
    x, m, s = norm(a)
    X2[:,n]=x
    M[n] = m
    S[n] = s

ncomp = 5

pca = PCA(n_components=ncomp)
trends = pca.fit_transform(X2)
weights = pca.components_

X3 = np.copy(X2)
for n in range(X2.shape[1]):
    for m in range(ncomp):
        X3[:,n] -= trends[:,m] * weights[m,n]

X4 = np.zeros_like(X3)
for n in range(X2.shape[1]):
    x = X3[:,n]
    a = unnorm(x, M[n], S[n])
    X4[:,n] = a

t=tpf.time[lkeep]


def close(event):
    plt.close('all')


def extract_LC(aperture):
    ax[0].cla()
    flux = X4[:,aperture.flatten()].sum(axis=1)
    m = np.nanmedian(flux)
    return t, flux/m


global mask
global aperture
mask = []


def onclick(event):
    global mask
    global aperture

    [p.remove() for p in reversed(ax[1].patches)]

    events = ((int(event.xdata+0.5), int(event.ydata+0.5)))

    if (len(mask) > 0) and (events in list(mask)): 
        mask = [x for x in mask if x != events] 
    else:
        mask.append(events)

    sqcol = '#ffffee'
    alpha = 0.5

    for pixel in mask:
        m = int(pixel[0])
        n = int(pixel[1])
        r = Rectangle((float(m)-0.5, float(n)-0.5), 1., 1., edgecolor='white', facecolor=sqcol, alpha = 0.5)
        ax[1].add_patch(r)
    
    fig.canvas.draw()

    #update the extraction aperture
    aperture = np.array(np.zeros_like(tpf.flux.sum(axis=0)), dtype=bool)
    for coord in mask:
        aperture[coord] = True

    [LC] = ax[0].plot(extract_LC(aperture)[0], extract_LC(aperture)[1],marker='o',color = '#054950', alpha = 0.9, lw = 0, markersize = 3, markerfacecolor='#054950')
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Normalized Flux")

    fig.canvas.draw_idle()

    plt.draw()


fig, ax = plt.subplots(1,2, figsize=(10,3), gridspec_kw={'width_ratios': [3, 1]})

fig.subplots_adjust(left=0.55, bottom=0.35)


plt.tight_layout()

ax[1].set_axis_off()

showflux = (tpf.flux).mean(axis = 0)

im = ax[1].imshow((tpf.flux).mean(axis = 0))
ax[0].set_xlabel("Time")
ax[0].set_ylabel("Normalized Flux")

fig.canvas.mpl_connect('button_press_event', onclick)

ebx = plt.axes([0.81, 0.04, 0.13, 0.06])
exit = Button(ebx, 'Close', color='orange')
exit.on_clicked(close)

plt.show()

#print (aperture)
#
#
#        
#        # -------- flatten the normal lighcurve --------
#        print ("Flatten LC...", end =" ")
#
#        l = np.isfinite(flux/m)
#
#        fr_inj = flux/m
#        alltime  = t
#
#        T_dur = 0.01  #may need to change!!
#        
#        nmed = int(720*3*T_dur)
#        nmed = 2*int(nmed/2)+1 # make it an odd number 
#        ff = filters.NIF(np.array(fr_inj),nmed,10,fill=True,verbose=True)
#        # first number is three time transit durations, the second quite small (10,20 )
#
#        l = np.isfinite(ff)
#        
#        g = interp1d(alltime[l],ff[l],bounds_error=False,fill_value=np.nan)
#        ff = g(alltime)
#        fr = fr_inj / ff
#
#        if save == True:
#
#            plt.figure(figsize=(16,5))
#            plt.plot(alltime,fr_inj + 0.95,'.')
#            plt.plot(alltime,ff + 0.95,'.')
#            plt.plot(alltime,fr,'o', color = 'r')
#
#            plt.savefig('{}/{}/{}_fit_test.png'.format(indir, tic, tic), format='png')
#            plt.clf()
#            plt.close()
#
#



