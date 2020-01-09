'''
filter.py: some simple routines to filter time-series.
'''

import numpy as np
from scipy import signal, interpolate
from norm import *

def filt1d(array, nmed, nlin, fill = False, circ = False, kernel = 'box'):
    '''
    Nonlinear (median+linear) filter with edge reflection.
    The 'kernel' parameter can set to either 'box' (boxcar, default)
    or 'gauss' (Gaussian)
    '''

    nmed = 2 * int(nmed/2) + 1
    N = len(array)
    if N < (3*nmed):
        print('filt1d:')
        print('    Warning: array too short compared to nmed')
        print('    Returning without doing anything')
        return array

    # Check for NaNs
    lnan = np.isnan(array)
    nnan = lnan.sum()
#    print ("nnan")
#    print (nnan)

    if nnan != 0:
        lnotnan = np.where(lnan==0)[0]
        if len(lnotnan) == 0:
            print('filt1d:')
            print('    Warning: no non-NaN data')
            print('    Returning without doing anything')
            return array
        med = np.median(array[lnotnan])
        # Fill in any data gaps by interpolating over them ...
        work = np.zeros(N)
        il = lnotnan.min()
        ih = lnotnan.max()
        xnew = np.arange(ih-il+1) + il
        f = interpolate.interp1d(lnotnan, array[lnotnan])  
        work[il:ih+1] = f(xnew)
        # ... or extending slope of nearest data points if at the edges
        if (il!=0):
            slope = work[il+1] - work[il]
            for i in range(il): work[i] = work[il] - slope*(il-i)
        if (ih!=N-1):
            slope = work[ih] - work[ih-1]
            for i in range(N-ih-1)+ih+1: work[i] = work[ih] + slope*(i-ih)
    else:
#        print ('nonans')
        work = np.copy(array)
    # Edge reflection
    nr = min(20, nlin)
    sz = max([nmed, nlin])
    if sz >= (N-1): 
        sz = N-2
    if circ != False:
        wl = work[N-sz:]
        wr = work[:sz]
    else:
        wl = array[0:sz]
        pivot = np.median(array[:nr])
        wl = 2 * pivot - wl
        wl = wl[::-1] # reverse array
        wr = array[N-sz:N]
        pivot = np.median(array[N-nr:N])
        wr = 2 * pivot - wr
        wr = wr[::-1] # reverse array
    work2 = np.zeros(N + 2 * sz)
    work2[0:sz] = wl
    work2[sz:sz+N] = work
    work2[sz+N:2*sz+N] = wr
    # Filter
    if nmed > 1: work = signal.medfilt(work2, nmed) 
    else: work = work2
    if kernel == 'gauss':
        box = signal.gaussian(10 * nlin, nlin)
    else:
        box = signal.boxcar(2*nlin+1)
    work = signal.convolve(work, box) / float(2*nlin+1)
    padd = int((len(work) - N - 2 * sz) / 2)
    #print(padd,N,sz)
    work = work[padd:padd+N+2*sz]
    # return to orginal array bounds
    result = work[sz:N+sz]
    # replace bad data if present
    if (fill==False) * (nnan!=0): result[lnan] = np.nan
    #print ("reached the end")
    return result

def boxcare(array, box, fill = False, circ = False):
    """Boxcar filter with edge reflection"""
    nlin = int(box / 2)
    return filt1d(array, 1, nlin, fill = fill, circ = circ)

def gausse(array, box, fill = False, circ = False):
    """Gaussian filter with edge reflection"""
    nlin = int(box / 2)
    return filt1d(array, 1, nlin, fill = fill, circ = circ, kernel = Gauss)

def NIF(array, nmed, nlin, nsig = 3.0, prefilt = False, \
            fill = True, circ = False, verbose = False):
    '''
    Non-linear iterative filter with k-sigma clipping (Aigrain & Irwin
    2004). If fill = True, clipped values are filled in by linear
    interpolation.
    '''
    # Pre-filter if requested

    if prefilt != False:
        work = filt1d(np.copy(array), 7, 3) 
    else:
        work = np.copy(array)
    # Start iteration loop
    irejo = -1
    irej = 0
    for i in range(5):
        if i > 0:
            irejo = irej
            out = np.isnan(work) + (myabs(array - work) > (nsig*sigma))
            keep = (out == False)
            irej = sum(out)
            work = np.copy(array)
            if irej != 0: work[out] = np.nan

        keep = np.isfinite(work)
        work[keep] = filt1d(work[keep], nmed, nlin, fill = fill, circ = circ)
        sigma = 1.48 * np.median(myabs(work[keep] - array[keep]))
        #if verbose == True: print('NIF:', i, irej, irejo, sigma)
        if irej <= irejo: break
    return work
