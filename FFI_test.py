from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import numpy as numpy
import pandas as pandas
from lightkurve import search_tesscut
import astropy.io.fits as pf
import matplotlib.pyplot as plt


#
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


def clicker6(event):

	global mask, aid, cid, did, eid, fid

	if event.inaxes:
		if event.key == '' or \
				event.key == 'x':

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
			fig.canvas.mpl_disconnect(aid)
			fig.canvas.mpl_disconnect(cid)
			fig.canvas.mpl_disconnect(did)
			fig.canvas.mpl_disconnect(eid)
			fig.canvas.mpl_disconnect(fid)
			plotimage()




def plotimage():
	print ("plotting im")
	ax = plt.axes([0.08,0.09,0.63,0.88])
	plt.subplots_adjust(0.06,0.1,0.93,0.88)
	
	plt.imshow((tpf.flux).mean(axis = 0))
	plt.gca().set_autoscale_on(False)
	plt.xlabel('Pixel Column Number', {'color' : 'k'})
	plt.ylabel('Pixel Row Number', {'color' : 'k'})
	#plt.grid()
	
	mask = []

	# plot the mask
	sqcol = '#ffffee'
	alpha = 0.8
	for pixel in mask:
		m = int(pixel.split(',')[0])
		n = int(pixel.split(',')[1])
		x = [m-0.5,m+0.5,m+0.5,m-0.5,m-0.5]
		y = [n-0.5,n-0.5,n+0.5,n+0.5,n-0.5]
		plt.fill(x,y,sqcol,alpha=alpha,ec=sqcol)
	fid = fig.canvas.mpl_connect('key_press_event',clicker6)
	plt.show()
	
	return


def main():


	search_result = lk.search_tesscut(searchtic, sector=sec)
	tpf = search_result.download(cutout_size=11)

	fig = plt.figure(figsize=(3.5,3.5))
	
	plotimage()

	return

if __name__ == "__main__":
	main()
