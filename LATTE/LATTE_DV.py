import os 
import numpy as np
import pandas as pd

from os.path import exists
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch

# modules to generate the a pdf validation report
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.colors import red, black, grey, darkcyan
from reportlab.lib.enums import TA_RIGHT, TA_JUSTIFY, TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image,  Table, TableStyle, PageBreak, Flowable


def LATTE_DV(tic, indir, syspath, transit_list, sectors_all, target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, tpf_corrupt, astroquery_corrupt, FFI, bls = False, model = False, mpi = False, test = 'no'):

	'''
	funtion that makes compiles all of the information and figures into a comprehensive pdf summary document.

    Parameters
    ----------
    tic  :   str
        target TIC ID
    indir  :  str
        path to directory where all the plots and data will be saved. 
    sectors_all  :   list
        all the sectors in which the target has been/ will be observed
    target_ra   :  float
        the right ascension of the target
    target_dec   :  float
        the declination of the target
    tessmag  :  float
        TESS magnitude of the target star
    teff  :  float
        effective temperature of the tagret star (K)
    srad  :  float
        radius of the target star (solar radii)
	bls_stats1  :  list
		list of the returned stats of the initial bls search
	bls_stats2  :  list
		list of the returned stats of the second bls search
	FFI  :  bool
		whether the input data is from FFIs (True = from FFIs)
	bls  :  bool  (false)
		whether the bls search was run 
	model  :  bool  (false)
		whether the transit was modelled (only works if payenti has sucessfully been installed)
    Returns
    -------
    LATTE Data Validation report in PDF format.

	'''

	# ---- CHECK WHETHER THE TARGET IS A TCE OR A TOI ----
	print ("\n Start compiling the data validation report...")

	# TCE -----
	lc_dv = np.genfromtxt('{}/data/tesscurl_sector_all_dv.sh'.format(indir), dtype = str)

	TCE_links = []

	for i in lc_dv:
		if str(tic) in str(i[6]):
			TCE_links.append(i[6])

	if len(TCE_links) == 0:
		TCE = " - "
	
	else:
		TCE_links = np.sort(TCE_links)
		TCE_link = TCE_links[0]  # this link should allow you to acess the MAST DV report
		TCE = 'Yes **'
		TCE_link = '<link href="%s" color="blue">here</link>' % TCE_link
	

	# TOI -----
	TOI_planets = pd.read_csv('{}/data/TOI_list.txt'.format(indir), comment = "#")
	
	TOIpl = TOI_planets.loc[TOI_planets['TIC'] == float(tic)]
		
	if len(TOIpl) == 0:
		TOI = ' -  '
	else:
		TOI = (float(TOIpl["Full TOI ID"]))
	
	# ------ PARAMS ------
	ra = float(target_ra)
	dec = float(target_dec)
	
	
	def addPageNumber(canvas, doc):
		"""
		Add the page numbers to the document
		"""
		width, height = A4 # this is useful when defining where to plot something on the page
	
		page_num = canvas.getPageNumber()
		text = "%s" % page_num
		header = "TIC {}".format(tic)
		
		canvas.setFont('Helvetica',8)
		canvas.drawString(width*0.85, height * 0.95, header)
		canvas.drawRightString(200*mm, 10*mm, text)
 

	#------------------------------------------------------------------
	# Recall the names of all the plots that will be in the DV report
	#------------------------------------------------------------------
	
	# plot the full light curve, with marked sectors and marked transit - binned and unbinned

	full_LC_name = "{}/{}/{}_fullLC_md.png".format(indir,tic,tic)
	
	background_flux_name = '{}/{}/{}_background.png'.format(indir, tic, tic)
	
	centroid_positions_name = '{}/{}/{}_centroids.png'.format(indir, tic, tic)
	
	flux_aperture_name= '{}/{}/{}_aperture_size.png'.format(indir, tic, tic)
	
	tess_stars_name = '{}/{}/{}_star_field.png'.format(indir, tic, tic)
	
	#SDSS_stars_name = '{}/{}/{}_SDSSstar_field.png'.format(indir, tic, tic)

	nearest_neighbour_name = '{}/{}/{}_nearest_neighbours.png'.format(indir, tic, tic)
	
	pixel_LCs_name = '{}/{}/{}_individual_pixel_LCs_0.png'.format(indir, tic,tic)
	
	bls1 = '{}/{}/{}_bls_first.png'.format(indir, tic, tic)
	
	bls2 = '{}/{}/{}_bls_second.png'.format(indir, tic, tic)
	
	in_out_name = '{}/{}/{}_flux_comparison.png'.format(indir, tic, tic)
	
	model_name = '{}/{}/model_out/{}b_tr.png'.format(indir, tic, tic)

	phasefold_name = '{}/{}/{}_phase_folded.png'.format(indir, tic, tic)
	
	apertures_name = '{}/{}/{}_apertures_0.png'.format(indir, tic, tic)
	
	# ----- LOGOS ------
	# if this is a unittest run, find the files for the logos stored in the test folder
	if test != 'no':
		PHT_logo_name  =  '{}/LATTE_imgs/PHT_logo.jpg'.format(test)
		LATTE_logo_name = '{}/LATTE_imgs/LATTE_logo.png'.format(test)
		TESS_logo_name  = '{}/LATTE_imgs/TESS_logo.png'.format(test)

	# otherwise they're located in the place where the program is insatlled. 
	else:
		PHT_logo_name  =  '{}/LATTE_imgs/PHT_logo.jpg'.format(syspath)
		LATTE_logo_name = '{}/LATTE_imgs/LATTE_logo.png'.format(syspath)
		TESS_logo_name  = '{}/LATTE_imgs/TESS_logo.png'.format(syspath)


	# -------------------------------------------
	# Make a PDF summary file
	# -------------------------------------------

	doc = SimpleDocTemplate("{}/{}/DV_report_{}.pdf".format(indir,tic,tic) ,pagesize=A4,
						   rightMargin=72,leftMargin=72,
						   topMargin=40,bottomMargin=20)

	width, height = A4 # this is useful when defining where to plot something on the page
	
	Story=[]

	fig_count = 0
	table_count = 0
	
	# title
	title = "PHT Data Validation Report"
	subheading = "TIC {}".format(tic)

	styles=getSampleStyleSheet()
	styles.add(ParagraphStyle(name='centre', alignment=TA_CENTER))
	ptext = '<font size=12><b>%s</b></font>' % title
	subtext = '<font size=12><b>%s</b></font>' % subheading
	
	Story.append(Paragraph(ptext, styles["centre"]))
	Story.append(Paragraph(subtext, styles["centre"]))
	Story.append(Spacer(1, 30))
	
	# ----- ADD THE LOGOS -------
	PHT_logo =   Image(PHT_logo_name)
	PHT_logo.drawHeight = 0.5*inch*PHT_logo.drawHeight / PHT_logo.drawWidth
	PHT_logo.drawWidth = 0.5*inch

	LATTE_logo = Image(LATTE_logo_name)
	LATTE_logo.drawHeight = 0.5*inch*LATTE_logo.drawHeight / LATTE_logo.drawWidth
	LATTE_logo.drawWidth = 0.5*inch

	TESS_logo =  Image(TESS_logo_name)
	TESS_logo.drawHeight = 0.8*inch*TESS_logo.drawHeight / TESS_logo.drawWidth
	TESS_logo.drawWidth = 0.8*inch


	logo_table = (Table([[PHT_logo, LATTE_logo, TESS_logo]],
					colWidths=[width * 0.1],
					rowHeights=[1 * mm]))

	logo_table.setStyle(TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTRE'),('VALIGN', (0, 0), (-1,-1), 'MIDDLE')]))
	
	Story.append(logo_table)

	Story.append(Spacer(1, 30))
	# ----------------------------

	line = MCLine(width*0.77)
	Story.append(line)
	Story.append(Spacer(1, 2))
	line = MCLine_color(width*0.77)
	Story.append(line)
	Story.append(Spacer(1, 20))

	# --------------------------------------------
	# Full Image with momentum dumps
	# --------------------------------------------
	im = Image(full_LC_name)
	im._restrictSize(width*0.8, width*0.8)
	Story.append(im)
	
	fig_count += 1

	full_image_text = "Fig {}. Full lightcurve for target TIC {}. The solid red lines at the bottom of the figure indicated the \
	times of the reaction wheel momentum dumps and the dashed black line(s) show the time(s) of the marked transit event(s). Momentum dumps \
	occur around every 2 to 2.5 days and typically last around half an hour.".format(fig_count,tic)
	
	ptext = '<font size=8>%s</font>' % full_image_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	
	# --------------------------------------
	# ------ stellar parameters table ------
	# --------------------------------------
	try:
		srad = float(srad)
	except:
		srad = -999
	
	try:
		teff = float(teff)
	except:
		teff = -999
			
	data_stellar= [['Parameter',  "Value", "Unit"],
				   ['TIC ID',	 tic, ],
				   ['RA',		 ra, "degrees"],
				   ['Dec',		dec, "degrees"],
				   ['Radius',	 "{:.4f}".format(srad), "Solar Radii"],
				   ['Tess Mag',   "{:.4f}".format(tessmag), "Mag"],
				   ['Teff',	   "{:.0f}".format(teff), "Kelvin"],
				   ['Sectors',	  "{} *".format(str(sectors_all)[1:-1]), ],
				   ['TCE',	  TCE, ],
				   ['TOI',	  "{}".format(str(TOI)), ],
				   ]
	
	table_stellar=Table(data_stellar)
	table_stellar=Table(data_stellar,colWidths=width * 0.2, style=[
						('LINEABOVE',(0,1),(-1,1),1,colors.black),
						('LINEABOVE',(0,10),(-1,10),1,colors.black),
						('FONTSIZE', (0,0),(-1,9), 8),
						])

	data_len = len(data_stellar)
	
	for each in range(data_len):
		if each % 2 == 0:
			bg_color = colors.whitesmoke
		else:
			bg_color = colors.white
	
		table_stellar.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))
	

	# ------ ADD A LINE TO SEPERATE SECTIONS -----

	Story.append(Spacer(1, 20))
	line = MCLine(width*0.77)
	Story.append(line)

	# ------
	
	Story.append(Spacer(1, 20))
	ptext = '<font size=9>Target Properties</font>'
	Story.append(Paragraph(ptext, styles["Normal"]))
	Story.append(Spacer(1, 12))
	
	Story.append(table_stellar)
	Story.append(Spacer(1, 15))
	
	table_count += 1

	exofop_url = "https://exofop.ipac.caltech.edu/tess/target.php?id={}".format(tic)
	exofop_link = '<link href="%s" color="blue">TIC %s</link>' % (exofop_url, tic)


	if TCE == 'Yes **':
		Stellartable_text = "Table {}. Stellar properties of {}. \
			* List of the sectors in which the target will be has been \
			observed. ** Click {} for the TCE report.".format(table_count, exofop_link, TCE_link)


	else:
		Stellartable_text = "Table {}. Stellar properties of the {}. \
			* List of the sectors in which the target will be has been \
			observed.".format(table_count,exofop_link)

		


	ptext = '<font size=8>%s</font>' % Stellartable_text
	Story.append(Paragraph(ptext, styles["Normal"]))


	# --------------------------------------------
	# Background
	# --------------------------------------------
	Story.append(PageBreak()) # always start a new page for this analysis
	im2 = Image(background_flux_name)
	
	if len(transit_list) == 1:
		im2._restrictSize(width*0.55, width*0.55)
	
	else:
		im2._restrictSize(width*0.8, width*0.8)
	
	Story.append(im2)

	fig_count += 1
	Story.append(Spacer(1, 1))
	background_text = "Fig {}. Background flux vs. time around the time of each transit-like event. \
		The vertical orange line indicates the time of the transit-like event.".format(fig_count)
	
	ptext = '<font size=8>%s</font>' % background_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	

	# --------------------------------------------
	# Centroid Position
	# --------------------------------------------
	
	if FFI == False:
		Story.append(Spacer(1, 10))
		im3 = Image(centroid_positions_name)

		if len(transit_list) == 1:
			im3._restrictSize(width*0.52, width*0.52)
		else:
			im3._restrictSize(width*0.7, width*0.7)
	
		Story.append(im3)
		
		fig_count += 1
		centroid_text = "Fig {}. The x and y centroid positions around the time of each transit-like event. The black points shows the CCD column and row position of the target’s flux-weighted centroid. \
			The red shows the CCD column and row local motion differential velocity aberration (DVA), pointing drift, and thermal effects. \
			The vertical orange line indicates the time of the transit-like event".format(fig_count)
		
		ptext = '<font size=8>%s</font>' % centroid_text

		Story.append(Spacer(1, 5))

		Story.append(Paragraph(ptext, styles["Normal"]))
		
		Story.append(Spacer(1, 13))

		#Story.append(PageBreak()) # always start a new page for this analysis

	# the following plots will only exist if the TPF file is not corrupt - otherwise skip these.
	if tpf_corrupt == False:

		# --------------------------------------------
		# Flux Aperture
		# --------------------------------------------
		im4 = Image(flux_aperture_name)
	
		if len(transit_list) == 1:
			im4._restrictSize(width*0.55, width*0.55)
		else:
			im4._restrictSize(width*0.7, width*0.7)
		Story.append(im4)
		
		fig_count += 1
		Story.append(Spacer(1, 10))
		flux_aperture_text = "Fig {}. The lightcurve around the time of each transit-like event extracted with the SPOC pipeline \
			defined aperture (binned:blue, unbinned:grey) and the with an aperture that is 40 per cent smaller (red). The flux is extracted \
			from the target pixel files (TPFs) and has not been detrended or \
			corrected for systematics. The vertical orange line indicates the time of the transit-like event.".format(fig_count)
		
		ptext = '<font size=8>%s</font>' % flux_aperture_text
		Story.append(Paragraph(ptext, styles["Normal"]))
		
		# --------------------------------------------
		# Apertures Sizes
		# --------------------------------------------
		im45 = Image(apertures_name)
		
		im45._restrictSize(width*0.4, width*0.4)
		
		Story.append(im45)
		
		fig_count += 1

		Story.append(Spacer(1, 10))

		if FFI == False:
			aperture_text = "Fig {}. The apertures used to extract the lightcurves. The blue aperture on the right shows the \
			optimum aperture determined by the SPOC pipeline, which is used for the extraction of 2-minute cadence light curves shown in Figure 1. \
			The red outline on the left shows an aperture that is around 40 per cent smaller than the SPOC pipeline aperture which was used to extract the \
			red lightcurve shown in Figure {}.".format(fig_count, (fig_count-1))
		else:
			aperture_text = "Fig {}. The larger (right hand side, blue) and the smaller (left hamd side, red) apertures used to extract the lightcurves shown in Figure {}.".format(fig_count, (fig_count-1))

		ptext = '<font size=8>%s</font>' % aperture_text
		Story.append(Paragraph(ptext, styles["Normal"]))
			
		
		# --------------------------------------------
		# In and Out of Transit Comparison
		# --------------------------------------------
		
		Story.append(Spacer(1, 12))
		im5 = Image(in_out_name)
	
	
		im5._restrictSize(width*0.9, width*0.9)
	
		Story.append(im5)
	
		fig_count += 1
		Story.append(Spacer(1, 10))
		flux_aperture_text = "Fig {}. Difference images for target TIC {} for each transit like event. \
		Left: mean in-transit flux(left). Middle: mean out-of-transit flux. Right: difference between the mean out-of-transit and mean in-transit flux. \
		Ensure that the change in brightness occurs on target.".format(fig_count, tic)
		
		ptext = '<font size=8>%s</font>' % flux_aperture_text
		Story.append(Paragraph(ptext, styles["Normal"]))
	
	
		# --------------------------------------------
		# tess stars + SDSS star field
		# --------------------------------------------
		# can only put this in the report if astroquery is working. 
		if astroquery_corrupt == False:

			Story.append(Spacer(1, 12))
		
			im6 = Image(tess_stars_name)
			
			# if not with mpi (two star images)
			if mpi == False:
				im6._restrictSize(width*0.7, width*0.5)
			else:
				im6._restrictSize(width*0.35, width*0.35)
	
			Story.append(im6)
		
			fig_count += 1
			tess_stars_text = "Fig {}. Left: The locations of nearby GAIA DR2 stars with mag < 15 (orange circle) within the Tess \
			Cut Out around TIC {} (red star). Only shown for one sector. Right: SDSS image of the surrounding field.".format(fig_count, tic)
			
			ptext = '<font size=8>%s</font>' % tess_stars_text
			Story.append(Paragraph(ptext, styles["Normal"]))
		
	
	# --------------------------------------------
	# nearest neighbours
	# --------------------------------------------
	#Story.append(PageBreak()) # always start a new page for this analysis
	if FFI == False:
		im7 = Image(nearest_neighbour_name)
		
		im7._restrictSize(width*0.8, width*0.8)
	
		Story.append(im7)
		fig_count += 1
		Story.append(Spacer(1, 10))
		nn_text = "Fig {}. Lightcurves of the five closest stars to target {} (top pannel). \
			The distances to the target star and the TESS magnitudes are shown for each star. Only ever shown for one sector.".format(fig_count,tic)
		
		ptext = '<font size=8>%s</font>' % nn_text
		Story.append(Paragraph(ptext, styles["Normal"]))
	

	# this plot also only exists if the TPF downloaded sucessfully
	if tpf_corrupt == False:
		# --------------------------------------------
		# pixel_LCs_name
		# --------------------------------------------
		Story.append(Spacer(1, 10))
		im8 = Image(pixel_LCs_name)
		
	
		im8._restrictSize(width*0.65, width*0.65)
	
		Story.append(im8)
		fig_count += 1
		pixLC_text = "Fig {}. Normalised flux extracted for each pixel, using the SPOC pipeline mask, around the time of the transit-like event. \
		The orange/red data points show the in-transit data. The solid red lines show the SPOC pipeline mask. Only shown for one sector.".format(fig_count)
		
		ptext = '<font size=8>%s</font>' % pixLC_text
		Story.append(Paragraph(ptext, styles["Normal"]))
	
	# ------ Phase Folded LC ------


	if len(transit_list) > 1:

		# --------------------------------------------
		# Phase Folded
		# --------------------------------------------
		impf = Image(phasefold_name)

		impf._restrictSize(width*0.35, width*0.35)

		Story.append(impf)
		
		fig_count += 1
		Story.append(Spacer(1, 10))
		flux_aperture_text = "Fig {}. Phase folded lightcurve where the odd and the even transits are shown in different colours. Ensure that the odd and even transits have comparabel shapes and depths.".format(fig_count)
		
		ptext = '<font size=8>%s</font>' % flux_aperture_text
		Story.append(Paragraph(ptext, styles["Normal"]))
		

	# ------ BLS -------
	Story.append(PageBreak()) # always start a new page for this analysis
	# ------

	if bls == True:
		
		Story.append(Spacer(1, 12))
		blsim1 = Image(bls1)
		blsim2 = Image(bls2)

		blsim1._restrictSize(width*0.6, width*0.6)
		blsim2._restrictSize(width*0.6, width*0.6)

		bls_table = (Table([[blsim1, blsim2]],
						colWidths=[width * 0.45], rowHeights=[width * 0.6]))
	
		bls_table.setStyle(TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTRE'),('VALIGN', (0, 0), (-1,-1), 'MIDDLE')]))
		
		Story.append(bls_table)

		fig_count += 1

		if FFI == False:

			bls1_text = "Fig {}. Box Least Square fitting (BLS) for whole lightcurve binned to 10 minutes. Top left panel: log liklihood periodogram. \
							The solid red line indicates the peak period and the dashed orange lines show the integer \
							harmonics of this period. Middle left panel: Full light curve, unbinned (orange) and binned to 10 minutes (black). \
							The peak period is highlighted by the solid red lines. Bottom left Panel: Phase folded light curve where the found transit-event is fit \
							with a simple box (red line). The pannels on the right show the same diagnostics, however the diagnostic \
							was run with the highest detected signal-to-noise transits, from the initial BLS search, removed. ".format(fig_count)
		
		else:
			bls1_text = "Fig {}. Box Least Square fitting (BLS) for whole lightcurve. Top left panel: log liklihood periodogram. \
							The solid blue line indicates the peak period and the dashed red lines show the integer \
							harmonics of this period. Middle left panel: Full light curve, unbinned LC (orange) . \
							The peak period is highlighted by the solid blue lines. Bottom left Panel: Phase folded light curve where the found transit-event is fit \
							with a simple box (blue line). The pannels on the right show the same diagnostics, however the diagnostic \
							was run with the highest detected signal-to-noise transits, from the initial BLS search, removed. ".format(fig_count)
		

		ptext = '<font size=8>%s</font>' % bls1_text
		Story.append(Paragraph(ptext, styles["Normal"]))
		


		# --------------------
		# ---- BLS TABLE -----

		data_bls= [['Parameter',				  "bls1",														"bls2"],
					   ['period',			    "{:.3f}".format(bls_stats1[0]),	         						 "{:.3f}".format(bls_stats2[0])],
					   ['t0',			        "{:.2f}".format(bls_stats1[1]),	         						 "{:.2f}".format(bls_stats2[1])],
					   ['depth',			    "{:.5f} ± {:.5f}".format(bls_stats1[2][0],bls_stats1[2][1]),	 "{:.5f} ± {:.5f}".format(bls_stats2[2][0],bls_stats2[2][1]) ],
					   ['depth phased',			"{:.5f} ± {:.5f}".format(bls_stats1[3][0],bls_stats1[3][1]),	 "{:.5f} ± {:.5f}".format(bls_stats2[3][0],bls_stats2[3][1]) ],
					   ['depth half',			"{:.5f} ± {:.5f}".format(bls_stats1[4][0],bls_stats1[4][1]),	 "{:.5f} ± {:.5f}".format(bls_stats2[4][0],bls_stats2[4][1]) ],
					   ['depth odd',			"{:.5f} ± {:.5f}".format(bls_stats1[5][0],bls_stats1[5][1]),	 "{:.5f} ± {:.5f}".format(bls_stats2[5][0],bls_stats2[5][1]) ],
					   ['depth even',			"{:.5f} ± {:.5f}".format(bls_stats1[6][0],bls_stats1[6][1]),	 "{:.5f} ± {:.5f}".format(bls_stats2[6][0],bls_stats2[6][1]) ],
					   ]
		
		table_bls=Table(data_bls)
		table_bls=Table(data_bls,colWidths=width * 0.2, style=[
							('LINEABOVE',(0,1),(-1,1),1,colors.black),
							('LINEABOVE',(0,8),(-1,8),1,colors.black),
							('FONTSIZE', (0,0),(-1,7), 8),
							])
	

		# ------ ADD A LINE TO SEPERATE SECTIONS -------
		
		Story.append(Spacer(1, 20))
		line = MCLine(width*0.77)
		Story.append(line)
		
		# ------

		Story.append(Spacer(1, 16))
		ptext = '<font size=10>BLS parameters</font>'
		Story.append(Paragraph(ptext, styles["Normal"]))
		Story.append(Spacer(1, 15))
	
		data_len = len(data_bls)
		
		for each in range(data_len):
			if each % 2 == 0:
				bg_color = colors.whitesmoke
			else:
				bg_color = colors.white
		
			table_bls.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))
		
		
		Story.append(table_bls)
		Story.append(Spacer(1, 15))

		table_count += 1
		Stellartable_text = "Table {}. Summary of the two BLS fits. Fit one is run with the whole lightcurve and fit two is run with the highest detected signal-to-noise transits removed.".format(table_count)
		ptext = '<font size=8>%s</font>' % Stellartable_text
		Story.append(Paragraph(ptext, styles["Normal"]))

		Story.append(PageBreak())
		# -----------


	if model == True:

		#Story.append(PageBreak()) # always start a new page for this analysis
		pyaneti_url = 'https://academic.oup.com/mnras/article/482/1/1017/5094600'

		pyaneti_link = '<link href="%s" color="blue">Pyaneti</link>' % pyaneti_url


		model_title = "Modeling"
		model_text = "The modeling of target TIC {} using the open source {} package.".format(tic, pyaneti_link)
		
		ptext = '<font size=11><b>%s</b></font>' % model_title
		Story.append(Paragraph(ptext, styles["centre"]))
		Story.append(Spacer(1, 10))
		ptext = '<font size=11>%s</font>' % model_text
		Story.append(Paragraph(ptext, styles["centre"]))
		Story.append(Spacer(1, 15))
		
		line = MCLine(width*0.77)
		Story.append(line)

		# --------------------------------------------
		# Pyaneti Modeling results
		# --------------------------------------------
		
		Story.append(Spacer(1, 12))
		im_model = Image(model_name)
			
		im_model._restrictSize(width*0.5, width*0.5)
			
		Story.append(im_model)
		fig_count += 1

		model_text = "Fig {}. The phase folded lightcurve from the the Pyaneti modeling. The solid black line shows the best fit model. See Table 2 for model parameters.".format(fig_count)
			
		ptext = '<font size=8>%s</font>' % model_text
		Story.append(Paragraph(ptext, styles["Normal"]))
		

		# ---------------
		# pyaneti modeling table
	
		# ----------------------------------------------
		# print a table of the model/pyaneti parameters
		# ----------------------------------------------
		
		# this can only be done if the pyaneti model has been run
		# as this takes quite a while this might not always be the case.


		if os.path.exists("{}/{}/model_out/{}_parameters.csv".format(indir, tic, tic)):
			
			Story.append(Spacer(1, 12))
			line = MCLine(width*0.77)
			Story.append(line)
			Story.append(Spacer(1, 12))

			ptext = '<font size=10>Candidate Model Parameters</font>'
			Story.append(Paragraph(ptext, styles["Normal"]))
			Story.append(Spacer(1, 8))
		
			manifest_table = pd.read_csv('{}/{}/model_out/{}_parameters.csv'.format(indir, tic, tic))
		
			#manifest_table
			params = ['T0', 'P', 'e', 'w', 'b', 'a/R*', 'rp/R*', 'Rp', 'Tperi', 'i', 'a', 'Insolation', 'rho*', 'g_p', 'Teq', 'T_tot', 'T_full']
			elements = []
		
			param_vals = []
			param_errs = []
			param_units= []
		
		
			for st in params:
				try:
					s = (manifest_table[manifest_table['Var']==st])
					param_vals.append("{:.3f}".format(s['Val'].values[0] ))
					param_errs.append("+{:.4f}  -{:.4f}".format(s['Pos'].values[0], s['Neg'].values[0]))
					param_units.append("{}".format(s['Unit'].values[0] ))
				except:
					param_vals.append(-999)
					param_errs.append(-999)
					param_units.append(-999)
		
			data_params = [['Parameters', 'Value', 'Uncertainty', 'Unit']]
			for p,v,e,u in zip(params, param_vals, param_errs, param_units):
				data_params.append([p,v,e,u])
		
			table_model=Table(data_params)
			table_model=Table(data_params,colWidths=width * 0.2, style=[
								('LINEABOVE',(0,1),(-1,1),1,colors.black),
								('LINEABOVE',(0,8),(-1,8),1,colors.black),
								('LINEABOVE',(0,len(params)+1),(-1,len(params)+1),1,colors.black),
								('FONTSIZE', (0,0),(-1,len(params)), 8),
								])
		

			data_len = len(data_params)
			
			for each in range(data_len):
				if each % 2 == 0:
					bg_color = colors.whitesmoke
				else:
					bg_color = colors.white
			
				table_model.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))
			

			Story.append(table_model)

			table_count += 1
			Story.append(Spacer(1, 10))
			Stellartable_text = "Table {}. The candidate paramaters from the Pyaneti modeling.".format(table_count)
			ptext = '<font size=8>%s</font>' % Stellartable_text
			Story.append(Paragraph(ptext, styles["Normal"]))
		


	doc.build(Story, onFirstPage=addPageNumber, onLaterPages=addPageNumber)


# ----- ADDITIONAL MODULES ------

class MCLine(Flowable):
	"""
	Line flowable --- draws a line in a flowable
	http://two.pairlist.net/pipermail/reportlab-users/2005-February/003695.html
	"""
 
	#----------------------------------------------------------------------
	def __init__(self, width, height=0):
		Flowable.__init__(self)
		self.width = width
		self.height = height
 
	#----------------------------------------------------------------------
	def __repr__(self):
		return "Line(w=%s)" % self.width
 
	#----------------------------------------------------------------------
	def draw(self):
		"""
		draw the line
		"""
		self.canv.setStrokeColor(grey)
		self.canv.line(0, self.height, self.width, self.height)
	  
class MCLine_color(Flowable):
	"""
	Line flowable --- draws a line in a flowable in COLOUR
	http://two.pairlist.net/pipermail/reportlab-users/2005-February/003695.html
	"""
 
	#----------------------------------------------------------------------
	def __init__(self, width, height=0):
		Flowable.__init__(self)
		self.width = width
		self.height = height
 
	#----------------------------------------------------------------------
	def __repr__(self):
		return "Line(w=%s)" % self.width
 
	#----------------------------------------------------------------------
	def draw(self):
		"""
		draw the line
		"""
		self.canv.setStrokeColor(darkcyan)
		self.canv.line(0, self.height, self.width, self.height)

