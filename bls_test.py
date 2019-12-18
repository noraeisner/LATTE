#!/usr/bin/env python

try:
    import pylab
    from pylab import *
    from matplotlib import *
except:
    print ('ERROR -- Cannot find the python module matplotlib in your executable path')
    print ('         download and install from http://matplotlib.sourceforge.net\n')
    usage()

# numpy

try:
    import numpy
    from numpy import *
    seterr(all="ignore") 
except:
    print ('ERROR -- Cannot find the python module numpy in your executable path')
    print ('         download and install from http://numpy.scipy.org\n')
    usage()

# pyfits

try:
    import pyfits
except:
    print ('ERROR -- Cannot find the python module pyfits in your executable path')
    print ('         download and install from http://www.stsci.edu/resources/software_hardware/pyfits\n')
    usage()

import getopt, sys, urllib, time, re

# global variables

ffifile = False; aperfile = False; maskfile = 'KeplerFFI.txt'
plotfile = 'KeplerFFI.png'; pimg = None; mask = []
xmin = 0.0; xmax = 1000.0; ymin = 0.0; ymax = 1000.0; zmin = False; zmax = False
kepid = ''; ra = ''; dec = ''; kepmag = ''; season = ''; quarter = -1
skygroup = ''; channel = ''; module = ''; output = ''; column = ''; row = ''
cmap = 'jet'; aid = None; cid = None; did = None; eid = None; fid = None
pylab_latex = True

# -----------------------------------------------------------
# core code

def main():

    global pimg, zmin, zmax, xmin, xmax, ymin, ymax, quarter
    global kepid, ra, dec, kepmag, skygroup, season, channel
    global module, output, row, column, cmap, maskfile, plotfile
    global pylab_latex

# input arguments

    #kepid,ra,dec,ffifile,zmin,zmax,zscale,cmap,npix = parameters()
    #if kepid:
    #    maskfile = re.sub('.txt',kepid+'.txt',maskfile)
    #    plotfile = re.sub('.png',kepid+'.png',plotfile)
    #else:
    #    status = 0

# open FFI FITS file

#    ffi, status = openfits(ffifile,'readonly')
#    try:
#        quarter = ffi[0].header['QUARTER']
#    except:
#        try:
#            dateobs = ffi[0].header['DATE-OBS']
#
#        except:
#            txt  = 'ERROR -- cannot determine quarter when FFI was taken. Either a  '
#            txt += 'QUARTER or DATE-OBS keyword is expected in the primary header'
#            sys.exit(txt)
#    if quarter == 0: quarter = 1
#    if quarter < 0:
#        txt  = 'ERROR -- cannot determine quarter from FFI. Try downloading a new\n'
#        txt += 'version of KeplerFFI.py from http://keplergo.arc.nasa.gov'
#        sys.exit(txt)
#    if int(quarter) == 0:
#        season = 3
#    else:
#        season = (int(quarter) - 2) % 4
#
# locate target in MAST


    #if kepid:
    #    kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column \
    #        = MASTKepID(kepid,season)
    #else:
    #     kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column \
    #        = MASTRADec(ra,dec,120.0,season)
    #     
    #     ra = float(int(ra * 1.0e5)) / 1e5
    #     dec = float(int(dec * 1.0e5)) / 1e5
     
# read and close FFI FITS file
	
    searchtic = 'TIC 229940491'  #349972099
    kepid = 229940491
    sec = 14
    
    search_result = lk.search_tesscut(searchtic, sector=sec)
    img = search_result.download(cutout_size=11)
    
    #try:
    #    img, status = readimage(ffi,int(channel))
    #    status = closefits(ffi)
    #except:
    #    txt  = 'ERROR -- target is not on a legal channel during season ' + str(season)
    #    sys.exit(txt)
#
## print target data
#
#    print ''
#    print '      KepID:  %s' % kepid
#    print ' RA (J2000):  %s' % ra
#    print 'Dec (J2000): %s' % dec
#    print '     KepMag:  %s' % kepmag
#    print '   SkyGroup:    %2s' % skygroup
#    print '     Season:    %2s' % str(season)
#    print '    Channel:    %2s' % channel
#    print '     Module:    %2s' % module
#    print '     Output:     %1s' % output
#    print '     Column:  %4s' % column
#    print '        Row:  %4s' % row
#    print ''
#
## subimage of channel for plot
#
    ymin = max([int(row)-npix/2,0])
    ymax = min([int(row)+npix/2+1,img.shape[0]])
    xmin = max([int(column)-npix/2,0])
    xmax = min([int(column)+npix/2+1,img.shape[1]])

# intensity scale

    nstat = 2; pixels = []
    for i in range(ymin,ymax):
        for j in range(xmin,xmax):
            pixels.append(img[i,j])
    pixels = array(sort(pixels),dtype=float32)
    if int(float(len(pixels)) / 10 + 0.5) > nstat:
        nstat = int(float(len(pixels)) / 10 + 0.5)
    if not zmin:
        zmin = median(pixels[:nstat])
    if not zmax:
        zmax = median(pixels[-nstat:])
    if 'log' in zscale:
        img = log10(img)
        zmin = log10(zmin)
        zmax = log10(zmax)
    if 'sq' in zscale:
        img = sqrt(img)
        zmin = sqrt(zmin)
        zmax = sqrt(zmax)
    pimg = img[ymin:ymax,xmin:xmax]

# plot limits

    ymin = float(ymin) - 0.5
    ymax = float(ymax) - 0.5
    xmin = float(xmin) - 0.5
    xmax = float(xmax) - 0.5

# plot style

    try:
        params = {'backend': 'png',
                  'axes.linewidth': 2.5,
                  'axes.labelsize': 24,
                  'axes.font': 'sans-serif',
                  'axes.fontweight' : 'bold',
                  'text.fontsize': 12,
                  'legend.fontsize': 12,
                  'xtick.labelsize': 16,
                  'ytick.labelsize': 16}
        pylab.rcParams.update(params)
        pylab_latex = True
    except:
        pylab_latex = False
    pylab.figure(1,figsize=[10,7])
    plotimage(pylab_latex)

    return

# -----------------------------------------------------------
# plot channel image

def plotimage(pylab_latex):

    global aid, cid, did, eid, fid

# print FFI and source location data on plot

    pylab.clf()
    #pylab.axes([0.73,0.09,0.25,0.4])
    #if pylab_latex:
    #    pylab.text(0.1,1.0,'      KepID: %s' % kepid,fontsize=12)
    #    pylab.text(0.1,0.9,' RA (J2000): %s' % ra,fontsize=12)
    #    pylab.text(0.1,0.8,'Dec (J2000): %s' % dec,fontsize=12)
    #    pylab.text(0.1,0.7,'     KepMag: %s' % kepmag,fontsize=12)
    #    pylab.text(0.1,0.6,'   SkyGroup: %2s' % skygroup,fontsize=12)
    #    pylab.text(0.1,0.5,'     Season: %2s' % str(season),fontsize=12)
    #    pylab.text(0.1,0.4,'    Channel: %2s' % channel,fontsize=12)
    #    pylab.text(0.1,0.3,'     Module: %2s' % module,fontsize=12)
    #    pylab.text(0.1,0.2,'     Output: %1s' % output,fontsize=12)
    #    pylab.text(0.1,0.1,'     Column: %4s' % column,fontsize=12)
    #    pylab.text(0.1,0.0,'        Row: %4s' % row,fontsize=12)
    #else:
    #    pylab.text(0.1,1.0,'KepID: %s' % kepid,fontsize=12)
    #    pylab.text(0.1,0.9,'RA (J2000): %s' % ra,fontsize=12)
    #    pylab.text(0.1,0.8,'Dec (J2000): %s' % dec,fontsize=12)
    #    pylab.text(0.1,0.7,'KepMag: %s' % kepmag,fontsize=12)
    #    pylab.text(0.1,0.6,'SkyGroup: %2s' % skygroup,fontsize=12)
    #    pylab.text(0.1,0.5,'Season: %2s' % str(season),fontsize=12)
    #    pylab.text(0.1,0.4,'Channel: %2s' % channel,fontsize=12)
    #    pylab.text(0.1,0.3,'Module: %2s' % module,fontsize=12)
    #    pylab.text(0.1,0.2,'Output: %1s' % output,fontsize=12)
    #    pylab.text(0.1,0.1,'Column: %4s' % column,fontsize=12)
    #    pylab.text(0.1,0.0,'Row: %4s' % row,fontsize=12)
    #pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    xlim(0.0,1.0)
    ylim(-0.05,1.12)

# clear button

    pylab.axes([0.73,0.87,0.25,0.09])
    pylab.text(0.5,0.5,'CLEAR',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    aid = connect('button_press_event',clicker1)

# dump custom aperture to file button

    pylab.axes([0.73,0.77,0.25,0.09])
    pylab.text(0.5,0.5,'DUMP',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    cid = connect('button_press_event',clicker3)

# print window to png file button

    pylab.axes([0.73,0.67,0.25,0.09])
    pylab.text(0.5,0.5,'PRINT',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    did = connect('button_press_event',clicker4)

# print window to png file button

    pylab.axes([0.73,0.57,0.25,0.09])
    pylab.text(0.5,0.5,'CLOSE',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    eid = connect('button_press_event',clicker5)

# plot the image window

    ax = pylab.axes([0.08,0.09,0.63,0.88])
    pylab.subplots_adjust(0.06,0.1,0.93,0.88)
    pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    labels = ax.get_yticklabels()
    setp(labels, 'rotation', 90)
    imshow(pimg,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zmin,vmax=zmax,extent=(xmin,xmax,ymin,ymax),cmap=cmap)
    pylab.gca().set_autoscale_on(False)
    xlabel('Pixel Column Number', {'color' : 'k'})
    ylabel('Pixel Row Number', {'color' : 'k'})
    grid()

# plot the mask

    if cmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                'gray','pink','RdGy']:
        sqcol = 'g'
        alpha = 0.5
    else:
        sqcol = '#ffffee'
        alpha = 0.8
    for pixel in mask:
        m = int(pixel.split(',')[0])
        n = int(pixel.split(',')[1])
        x = [m-0.5,m+0.5,m+0.5,m-0.5,m-0.5]
        y = [n-0.5,n-0.5,n+0.5,n+0.5,n-0.5]
        pylab.fill(x,y,sqcol,alpha=alpha,ec=sqcol)
    fid = connect('key_press_event',clicker6)
    show()
    sys.exit()
    return

# -----------------------------------------------------------
# get input parameters

def parameters():

    global mask, cmap

    kepid = False
    ra = False
    dec = False
    quarter = False
    zmin = False
    zmax = False
    zscale = 'log'
    npix = 30
    decimal = True

    try:
		opts, args = getopt.getopt(sys.argv[1:],"h:krdfayzscn",
                                   ["help","kepid=","ra=","dec=","ffifile=",
                                    "aperfile=","zmin=","zmax=","zscale=",
                                    "cmap=","npix="])
    except getopt.GetoptError:
		usage()
    for o, a in opts:
		if o in ("-h", "--help"):
	    	usage()

# kepid

	if o in ("-k", "--kepid"):
            try:
                kepid = int(a)
                if kepid <= 0:
                    txt = 'ERROR -- kepid is not a positive integer'
                    sys.exit(txt)
            except:
                txt = 'ERROR -- kepid is not an integer'
                sys.exit(txt)
            kepid = str(kepid)

# ra coordinate

	if o in ("-r", "--ra"):
            try:
                ra = float(a)
                if ra < 0.0 or ra > 360.0:
                    txt = 'ERROR -- ra is not a sensible positive floating point number'
                    sys.exit(txt)
            except:
                decimal = False
                ra = str(a)

# dec coordinate

        if o in ("-d", "--dec"):
            try:
                dec = float(a)
                if dec < 0:
                    txt = 'ERROR -- dec is not a sensible floating point number'
                    sys.exit(txt)
            except:
                decimal = False
                dec = str(a)

# ffifile

        if o in ("-f", "--ffifile"):
            ffifile = str(a)

# aperfile: an aperture definition file

        if o in ("-a", "--aperfile"):
            aperfile = str(a)
            if os.path.isfile(aperfile):
                out = open(aperfile,'r')
                for line in out:
                    if 'HALO' in line and 'UNDERSHOOT' in line:
                        line = line.strip().split('|')
                        r = int(line[3])
                        c = int(line[4])
                        for pixel in line[5].split(';'):
                            y = int(pixel.split(',')[0])
                            x = int(pixel.split(',')[1])
                            mask.append(str(x+c) + ',' + str(y+r))
                out.close()
            else:
                txt = 'ERROR -- a custom aperture file ' + aperfile + ' doesn''t exist'
                sys.exit()

# zmin: minimum intensity level

        if o in ("-y", "--zmin"):
            try:
                zmin = float(a)
                if zmin < 0:
                    txt = 'ERROR -- zmin is not a positive floating point number'
                    sys.exit(txt)
                if zmin == 0:
                    zmin = 1 
            except:
                txt = 'ERROR -- zmin is not a floating point numnber'
                sys.exit(txt)

# zmax: maximum intensity level

        if o in ("-z", "--zmax"):
            try:
                zmax = float(a)
                if zmax <= 0:
                    txt = 'ERROR -- zmax is not a positive floating point number'
                    sys.exit(txt)
                if zmax <= zmin:
                    txt = 'ERROR -- zmax <= zmin'
                    sys.exit(txt)
            except:
                txt = 'ERROR -- kepid is not a floating point numnber'
                sys.exit(txt)

# zscale: form of image intensity scale

        zscales = ['lin','log','sq']
        if o in ("-s", "--zscale"):
            y = str(a)
            for z in zscales:
                if z in y:
                    zscale = z
            if zscale not in zscales:
                txt = 'ERROR -- not a recognized image intensity scale'
                sys.exit(txt)

# cmap: pylab colormap for image

        if o in ("-c", "--cmap"):
            cmap = str(a)
            if cmap == 'browse':
                cmap_plot()
            else:
                cmaps=[m for m in cm.datad if not m.endswith("_r")]
                if cmap not in cmaps:
                    txt = 'ERROR -- cmap is not a valid colormap. Try --cmap=browse' 
                    sys.exit(txt)

# npix: sub-image dimension

        if o in ("-n", "--npix"):
            try:
                npix = int(a)
                if npix <= 1:
                    txt = 'ERROR -- npix is too small'
                    sys.exit(txt)
            except:
                txt = 'ERROR -- npix is not a postive integer'
                sys.exit(txt)

# decimal coordinate conversion

    if not decimal:
        try:
            ra,dec = sex2dec(ra,dec)
        except:
            if not kepid: 
                txt = 'ERROR -- no sensible RA and Dec coordinates provided'
                sys.exit(txt)

# check we have all input parameters

    if not kepid and (not ra or not dec):
        txt = 'ERROR -- kepid or coordinates have not been supplied'
        sys.exit(txt)
    if not ffifile:
        txt = 'ERROR -- ffifile has not been supplied'
        sys.exit(txt)
    if not os.path.isfile(ffifile):
        txt = 'ERROR -- can''t find FFI file ' + ffifile
        sys.exit(txt)

    return kepid,ra,dec,ffifile,zmin,zmax,zscale,cmap,npix

# -----------------------------------------------------------
# target data retrieval from MAST based upon KepID

def MASTKepID(id,season):

    global skygroup, column, row

# build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&kic_kepler_id=' + id
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=CSV'

# retrieve results from MAST

    out = ''
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')
    if len(out) > 0:
        kepid = out[0]
        ra = out[4]
        dec = out[5]
        #kepmag = out[19]
        #skygroup = out[45]
        #channel = out[67 + season * 5]
        #module = out[68 + season * 5]
        #output = out[69 + season * 5]
        #row = out[70 + season * 5]
        #column = out[71 + season * 5]
    else:
        txt = 'ERROR -- no target found with KepID %s' % id
        sys.exit(txt)

    return kepid,ra,dec #,kepmag,skygroup,channel,module,output,row,column

# -------------------------------------
# detector location retrieval based upon RA and Dec

def MASTRADec(ra,dec,darcsec,season):

    global skygroup, column, row

# WCS data

    cd1_1 = 0.000702794927969
    cd1_2 = -0.000853190160515
    cd2_1 = -0.000853190160515
    cd2_2 = -0.000702794927969
    cd = array([[cd1_1,cd1_2],[cd2_1,cd2_2]])
    cd = linalg.inv(cd)

# coordinate limits

    x1 = 1.0e30
    x2 = x1
    darcsec /= 3600.0
    ra1 = ra - darcsec / 15.0 / cos(dec * pi / 180)
    ra2 = ra + darcsec / 15.0 / cos(dec * pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

# build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&kic_degree_ra=' + str(ra1) + '..' + str(ra2)
    url += '&kic_dec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=CSV'

# retrieve results from MAST: nearest KIC source to supplied coordinates

    z = ''
    x = 1.0e30
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')
            r = (float(out[4].split(' ')[0]) + \
                float(out[4].split(' ')[1]) / 60.0 + \
                float(out[4].split(' ')[2]) / 3600.0) * 15.0
            d = float(out[5].split(' ')[0]) + \
                float(out[5].split(' ')[1]) / 60.0 + \
                float(out[5].split(' ')[2]) / 3600.0
            a = sqrt((abs(r - ra) / 15.0 / cos(d * pi / 180))**2 + abs(d - dec)**2)
            if a < x:
                x = a
                z = line.split(',')

    if len(z) > 0:
        kepid = None
        kepmag = None
        skygroup = out[45]
        channel = out[67 + season * 5]
        module = out[68 + season * 5]
        output = out[69 + season * 5]
    else:
        txt = 'ERROR -- row and column could not be calculated. Is location on silicon?'
        sys.exit(txt)
            
# convert coordinates to decimal for the two targets, determine distance from input

    zra,zdec = sex2dec(z[4],z[5])
    dra = zra - ra
    ddec = zdec - dec
    drow = cd[0,0] * dra + cd[0,1] * ddec
    dcol = cd[1,0] * dra + cd[1,1] * ddec
   
# pixel coordinates of the nearest KIC target

    row = z[70 + season * 5]
    column = z[71 + season * 5]

# pixel coordinate of target

    row = str(int(float(row) + drow + 0.5))
    column = str(int(float(column) + dcol + 0.5))

    return kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column

# -----------------------------------
# convert sexadecimal hours to decimal degrees

def sex2dec(ra,dec):

    ra = re.sub('\s+','|',ra.strip())
    ra = re.sub(':','|',ra.strip())
    ra = re.sub(';','|',ra.strip())
    ra = re.sub(',','|',ra.strip())
    ra = re.sub('-','|',ra.strip())
    ra = ra.split('|')
    outra = (float(ra[0]) + float(ra[1]) / 60 + float(ra[2]) / 3600) * 15.0

    dec = re.sub('\s+','|',dec.strip())
    dec = re.sub(':','|',dec.strip())
    dec = re.sub(';','|',dec.strip())
    dec = re.sub(',','|',dec.strip())
    dec = re.sub('-','|',dec.strip())
    dec = dec.split('|')
    if float(dec[0]) > 0.0:
        outdec = float(dec[0]) + float(dec[1]) / 60 + float(dec[2]) / 3600
    else:
        outdec = float(dec[0]) - float(dec[1]) / 60 - float(dec[2]) / 3600

    return outra, outdec

# -----------------------------------------------------------
# open HDU structure

def openfits(file,mode):

    status = 0
    try:
        struct = pyfits.open(file,mode=mode)
    except:
	txt = 'ERROR -- cannot open ' + file + ' as a FITS file'
	sys.exit(txt)
    return struct, status

# -----------------------------------------------------------
# close HDU structure

def closefits(struct):

    status = 0
    try:
	struct.close()
    except:
	txt = 'ERROR -- cannot close HDU structure'
	sys.exit(txt)
    return status

# -----------------------------------------------------------
# read image from HDU structure

def readimage(struct,hdu):

    status = 0
    try:
        imagedata = struct[hdu].data
    except:
	txt = 'ERROR -- cannot read image data from HDU ' + str(hdu)
	sys.exit(txt)
    return imagedata, status

# -----------------------------------------------------------
# clear all pixels from pixel mask

def clicker1(event):

    global mask, aid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 488 and event.y < 537):
                disconnect(aid)
                disconnect(cid)
                disconnect(did)
                disconnect(eid)
                disconnect(fid)
                mask = []
                pylab.clf()
                plotimage(pylab_latex)

    return

# -----------------------------------------------------------
# dump custom aperture definition file

def clicker3(event):

    global aid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 432 and event.y < 480):
                masktxt  = 'NEW|'
                masktxt += skygroup + '|'
                masktxt += '{' + re.sub('\s+',':',str(ra))
                masktxt += ',' + re.sub('\s+',':',str(dec))
                masktxt += '},TAD_NO_HALO,TAD_NO_UNDERSHOOT_COLUMN|'
                masktxt += row + '|'
                masktxt += column + '|'
                for coord in sorted(set(mask)):
                    masktxt += str(int(coord.split(',')[1]) - int(row)) + ','
                    masktxt += str(int(coord.split(',')[0]) - int(column)) + ';'
                if (os.path.isfile(maskfile)):
                    os.remove(maskfile)
                out = open(maskfile,'a')
                out.write(masktxt[:-1]+'\n')
                out.close()
                print ('Wrote custom aperture definition file ' + maskfile + ' with ' + str(len(mask)) + ' pixels')
    return

# -----------------------------------------------------------
# print plot to png with left-mouse click

def clicker4(event):

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 377 and event.y < 425):
                pylab.savefig(plotfile)
                print ('Wrote plot hardcopy file ' + plotfile)
    return

# -----------------------------------------------------------
# close plot and exit program

def clicker5(event):

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 320 and event.y < 368):
                pylab.close('all')
    return

# -----------------------------------------------------------
# this function will be called with every click of the mouse

def clicker6(event):

    global mask, aid, cid, did, eid, fid

    if event.inaxes:
        if event.key == '' or \
                event.key == 'x':
            if cmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                        'gray','pink','RdGy']:
                sqcol = 'g'
                alpha = 0.5
            else:
                sqcol = '#ffffee'
                alpha = 0.8
            m = float(int(event.xdata + 0.5))
            n = float(int(event.ydata + 0.5))
            txt = str(int(m))+','+str(int(n))
            if txt in mask:
                tmpmask = []
                for pixel in mask:
                    if pixel != txt:
                        tmpmask.append(pixel)
                mask = tmpmask
            else:
                mask.append(txt)
            disconnect(aid)
            disconnect(cid)
            disconnect(did)
            disconnect(eid)
            disconnect(fid)
            plotimage(pylab_latex)

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot():

    pylab.figure(1,figsize=[5,10])
    a=outer(ones(10),arange(0,1,0.01))
    subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        subplot(l,1,i+1)
        pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        pylab.text(100.85,0.5,m,fontsize=10)
    show()
    sys.exit()


#-------------------------------
# main

if __name__ == "__main__":
    main()