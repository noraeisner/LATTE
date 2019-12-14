from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import numpy as numpy
import pandas as pandas
from lightkurve import search_tesscut
import astropy.io.fits as pf
import matplotlib.pyplot as plt



searchtic = 'TIC 229940491'  #349972099
kepid = 229940491
sec = 14

#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_title('click on points')
#

search_result = lk.search_tesscut(searchtic, sector=sec)
tpf = search_result.download(cutout_size=11)


mask = []

def onclick(event):
    mask.append("{},{}".format(int(event.xdata+0.5), int(event.ydata+0.5)))
    #print ("pixel {} {}".format(int(event.xdata+0.5), int(event.ydata+0.5)))
    print ("mask {}".format(mask))
    # plot the mask
    sqcol = '#ffffee'
    alpha = 0.5

    for pixel in mask:
        m = int(pixel.split(',')[0])
        n = int(pixel.split(',')[1])
        x = [m-0.5,m+0.5,m+0.5,m-0.5,m-0.5]
        y = [n-0.5,n-0.5,n+0.5,n+0.5,n-0.5] 
        plt.fill(x,y,sqcol,alpha=alpha,ec=sqcol)
        print (x,y)
    

fig, ax = plt.subplots(figsize=(7,7))
plt.imshow((tpf.flux).mean(axis = 0))
plt.fill(5,5,'#ffffee',alpha=1,ec='#ffffee')
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()









