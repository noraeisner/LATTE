'''
filter.py: some simple routines to filter time-series.
'''
'''
norm.py: NaN-robust wrappers for mean, median, min, max... + MAD
scatter estimates + simple normalisation routines
'''

import numpy as np
from scipy import signal, interpolate

def medsig(array):
    '''
    Return median and outlier-robust estimate of standard deviation
    (1.48 x median of absolute deviations).
    '''
    l = np.isfinite(array)
    if sum(l) == 0:
        return np.nan, np.nan
    if sum(l) == 1:
        return array[l], np.nan
    med = np.median(array[l])
    sig = 1.48 * np.median(abs(array[l] - med))
    return med, sig
  
def mysum(array, axis = None):
    '''
    Return sum of array along specified axis (ignoring NaNs)
    '''
    arr = np.copy(array)
    l = np.isfinite(arr)
    if l.any() == False: return np.nan
    arr[l==False] = 0.0
    res = np.sum(arr, axis = axis)
    if axis != None:
        l_sum = np.sum(l, axis = axis)
        res[l_sum==0] = np.nan
    return res

def mymean(array):
    '''
    Return mean of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return np.nan
    return(array[l].mean())

def mystd(array):
    '''
    Return standard deviation of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return np.nan
    return(array[l].std())

def mymin(array):
    '''
    Return minimum of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return np.nan
    return(array[l].min())

def mymax(array):
    '''
    Return maximum of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return np.nan
    return(array[l].max())

def myargmin(array):
    '''
    Return index of minimum of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return -1
    return(array[l].argmin())

def myargmax(array):
    '''
    Return index maximum of array (ignoring NaNs)
    '''
    l = np.isfinite(array)
    if l.any() == False: return -1
    return(array[l].argmax())

def myabs(array):
    l = np.isfinite(array)
    a = np.copy(array)
    a[l] = abs(array[l])
    return a

def norm01(array):
    '''
    Normalize array to zero mean and unit variance
    '''
    med, sig = medsig(array)
    a = np.copy(array)
    a -= med
    a /= sig
    return a, med, sig

def unnorm01(array, med, sig):
    '''
    Unnormalize array
    '''
    a = np.copy(array) 
    a *= sig
    a += med
    return a

def scale01(array):
    '''
    Scale array to range from 0 to 1
    '''
    nmeas, nobj = np.shape(array)
    arr = np.zeros((nmeas, nobj)) + np.nan
    mins = np.zeros(nobj)
    maxs = np.zeros(nobj)
    for i in np.arange(nobj):
        cur = array[:,i]
        mins[i] = mymin(cur)
        maxs[i] = mymax(cur)
        if np.isfinite(mins[i]) * np.isfinite(maxs[i]) * (maxs[i] > mins[i]):
            arr[:,i] = (cur - mins[i]) / (maxs[i] - mins[i])
    return arr, mins, maxs

def unscale(array):
    '''
    Unscale array from 0 to 1 back to original range
    '''
    nmeas, nobj = np.shape(array)
    arr = np.zeros((nmeas, nobj)) + np.nan
    for i in np.arange(nobj):
        arr[:,i] = (array[:,i] * (maxs[i] - mins[i])) + mins[i]
    return arr


# ---------

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
