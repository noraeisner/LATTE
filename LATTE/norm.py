'''
norm.py: NaN-robust wrappers for mean, median, min, max... + MAD
scatter estimates + simple normalisation routines
'''

import numpy as np

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

