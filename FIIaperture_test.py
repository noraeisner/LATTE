
import lightkurve as lk
import numpy as numpy
import pandas as pandas
from lightkurve import search_tesscut
import astropy.io.fits as pf
import matplotlib.pyplot as plt

searchtic = 'TIC 229940491'  #349972099
kepid = 229940491
sec = 14
plt_latex = True


# -----------------------------------------------------------
# clear all pixels from pixel mask

def clicker1(event):

	global mask, aid, cid, did, eid, fid

	if event.inaxes:
		if event.button == 1:
			if (event.x > 585 and event.x < 783 and
				event.y > 488 and event.y < 537):
				fig.canvas.mpl_disconnect(aid)
				fig.canvas.mpl_disconnect(cid)
				fig.canvas.mpl_disconnect(did)
				fig.canvas.mpl_disconnect(eid)
				fig.canvas.mpl_disconnect(fid)
				mask = []
				plt.clf()
				plotimage(plt_latex)

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
# close plot and exit program

def clicker5(event):

	if event.inaxes:
		if event.button == 1:
			if (event.x > 585 and event.x < 783 and
				event.y > 320 and event.y < 368):
				plt.close('all')
	return

# -----------------------------------------------------------
# this function will be called with every click of the mouse

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
			plotimage(plt_latex)

# -----------------------------------------------------------
# these are the choices for the image colormap


# plot the image

def main():

	global pimg, zmin, zmax, xmin, xmax, ymin, ymax, quarter
	global kepid, ra, dec, kepmag, skygroup, season, channel
	global module, output, row, column, cmap, maskfile, plotfile
	global pylab_latex, fig, tpf


	search_result = lk.search_tesscut(searchtic, sector=sec)
	tpf = search_result.download(cutout_size=11)

	fig = plt.figure(figsize=(3.5,3.5))
	plotimage(plt_latex)

	return


def plotimage(plt_latex):
	global aid, cid, did, eid, fid

	plt.clf()


	#plt.tight_layout()
	##plt.subplots_adjust(wspace=0.000001)
	
	#plt.subplot(111)
	#plt.axis('off')
	#plt.imshow((tpf.flux).mean(axis = 0))
	
	#plt.title("In Transit Flux (e-/candence)")
	
	
	# clear button
	plt.axes([0.73,0.87,0.25,0.09])
	plt.text(0.5,0.5,'CLEAR',fontsize=24,weight='heavy',
			   horizontalalignment='center',verticalalignment='center')
	plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
	plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
	#xlim(0.0,1.0)
	#ylim(0.0,1.0)
	aid = fig.canvas.mpl_connect('button_press_event',clicker1)
	
	# dump custom aperture to file button
	
	plt.axes([0.73,0.77,0.25,0.09])
	plt.text(0.5,0.5,'DUMP',fontsize=24,weight='heavy',
			   horizontalalignment='center',verticalalignment='center')
	plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
	plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
	#xlim(0.0,1.0)
	#ylim(0.0,1.0)
	cid = fig.canvas.mpl_connect('button_press_event',clicker3)
	

	# print window to png file button
	
	plt.axes([0.73,0.57,0.25,0.09])
	plt.text(0.5,0.5,'CLOSE',fontsize=24,weight='heavy',
			   horizontalalignment='center',verticalalignment='center')
	plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
	plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
	#xlim(0.0,1.0)
	#ylim(0.0,1.0)
	eid = fig.canvas.mpl_connect('button_press_event',clicker5)
	
	# plot the image window
	
	ax = plt.axes([0.08,0.09,0.63,0.88])
	plt.subplots_adjust(0.06,0.1,0.93,0.88)
	plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
	plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
	labels = ax.get_yticklabels()
	plt.setp(labels, 'rotation', 90)
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
#-------------------------------
# main

if __name__ == "__main__":
	main()


