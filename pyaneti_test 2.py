import os 
import pandas as pd

from os.path import exists
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_RIGHT, TA_JUSTIFY, TA_CENTER
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image,  Table, TableStyle, PageBreak, Flowable
from reportlab.lib.units import mm
from reportlab.lib.colors import red, black, grey, darkcyan

tic = 55525572
indir = "/Users/Nora/Documents/research/TESS/planethunters/LATTE"
peak_list = [1622.49] 
sectors_all = [11]
target_ra = -99
target_dec = -99
tessmag = -99
teff = -99
srad = -99
bls_stats1 = [-99]
bls_stats2 = [-99]
bls = True
model = True



# ---------------- PARAMS
ra = float(target_ra)
dec = float(target_dec)

print ("RA: {}".format(ra))
print ("DEC: {}".format(dec))
print ("rad: {}".format(srad))
print ("Tmag: {}".format(tessmag))
print ("Teff: {}".format(teff))

print ("tess-point sectors {}".format(sectors_all))


#------------------------------------------------------------------
# Run all the plotting routines that save the figures that we care about.
#------------------------------------------------------------------

# plot the full light curve, with marked sectors and marked transit - binned and unbinned
print ("Creating plots...")

full_LC_name = "{}/{}/{}_fullLC_md.png".format(indir,tic,tic)

background_flux_name = '{}/{}/{}_background.png'.format(indir, tic, tic)

centroid_positions_name = '{}/{}/{}_centroids.png'.format(indir, tic, tic)

flux_aperture_name= '{}/{}/{}_aperture_size.png'.format(indir, tic, tic)

tess_stars_name = '{}/{}/{}_star_field.png'.format(indir, tic, tic)

nearest_neighbour_name = '{}/{}/{}_nearest_neighbours.png'.format(indir, tic, tic)

pixel_LCs_name = '{}/{}/{}_individual_pixel_LCs_0.png'.format(indir, tic,tic)

bls1 = '{}/{}/{}_bls_first.png'.format(indir, tic, tic)

bls2 = '{}/{}/{}_bls_second.png'.format(indir, tic, tic)

in_out_name = '{}/{}/{}_flux_comparison.png'.format(indir, tic, tic)

model_name = '{}/{}/model_out/{}b_tr.png'.format(indir, tic, tic)

# ----- LOGOS ------
PHT_logo_name  =  './LATTE_imgs/PHT_logo.jpg'
LATTE_logo_name = './LATTE_imgs/LATTE_logo.jpg'
TESS_logo_name  = './LATTE_imgs/TESS_logo.png'


class MCLine(Flowable):
    """
    Line flowable --- draws a line in a flowable
    http://two.pairlist.net/pipermail/reportlab-users/2005-February/003695.html
    """
 
    #----------------------------------------------------------------------
    def __init__(self, width, color, height=0):
        Flowable.__init__(self)
        self.color = color
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

        self.canv.setStrokeColor(color)
        self.canv.line(0, self.height, self.width, self.color, self.height)
        


def addPageNumber(canvas, doc):
    """
    Add the page number
    """
    width, height = A4 # this is useful when defining where to plot something on the page

    page_num = canvas.getPageNumber()
    text = "%s" % page_num
    header = "TIC {}".format(tic)
    
    canvas.setFont('Helvetica',11)
    canvas.drawString(width*0.85, height * 0.95, header)
    canvas.drawRightString(200*mm, 20*mm, text)
 

def LATTE_DV(tic, indir, peak_list, sectors_all,target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, bls = False, model = False):


	# -------------------------------------------
	# Make a PDF summary file
	# -------------------------------------------

	doc = SimpleDocTemplate("{}/{}/DV_report_{}.pdf".format(indir,tic,tic) ,pagesize=A4,
						   rightMargin=72,leftMargin=72,
						   topMargin=25,bottomMargin=18)


	width, height = A4 # this is useful when defining where to plot something on the page
	
	Story=[]

	fig_count = 0

	# title
	title = "PHT Data Validation Report"
	subheading = "TIC {}".format(tic)


	styles=getSampleStyleSheet()
	styles.add(ParagraphStyle(name='centre', alignment=TA_CENTER))
	ptext = '<font size=12><b>%s</b></font>' % title
	subtext = '<font size=12><b>%s</b></font>' % subheading
	
	Story.append(Paragraph(ptext, styles["centre"]))
	Story.append(Paragraph(subtext, styles["centre"]))
	Story.append(Spacer(1, 60))
	

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

	chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTRE'),
                          ('VALIGN', (0, 0), (-1,-1), 'MIDDLE')])


	Story.append(Table([[PHT_logo, LATTE_logo, TESS_logo],],
					colWidths=[width * 0.07],
                    rowHeights=[1 * mm]))

	Story.append(Spacer(1, 10))
	# ----------------------------

	line = MCLine(width*0.75, grey)
	Story.append(line)
	Story.append(Spacer(1, 2))
	line = MCLine(width*0.75, grey)
	Story.append(line)
	Story.append(Spacer(1, 10))

	# --------------------------------------------
	# Full Image with momentum dumps
	# --------------------------------------------
	im = Image(full_LC_name)
	im._restrictSize(width*0.8, width*0.8)
	Story.append(im)
	
	fig_count += 1
	full_image_text = "Fig {}. Full lightcurve for target TIC {}. The solid red lines at the bottom of the figure indicated the \
	times of the momentum dumps and the dashed black line(s) show the time(s) of the marked transit event(s).".format(fig_count,tic)
	

	ptext = '<font size=8>%s</font>' % full_image_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	
	# --------------------------------------
	# ------ stellar parameters table ------
	# --------------------------------------

	data_stellar= [['Parameter',  "Value", "Unit"],
				   ['TIC ID',	 tic, ],
				   ['RA',		 ra, "degrees"],
				   ['Dec',		dec, "degrees"],
				   ['Radius',	 "{:.4f}".format(srad), "Solar Radii"],
				   ['Tess Mag',   "{:.4f}".format(tessmag), "Mag"],
				   ['Teff',	   "{:.0f}".format(teff), "Kelvin"],
				   ['Sectors *',	  "{}".format(str(sectors_all)[1:-1]), ],
				   ]
	
	table_stellar=Table(data_stellar)
	table_stellar=Table(data_stellar,colWidths=width * 0.2, style=[
						('LINEABOVE',(0,1),(-1,1),1,colors.black),
						('LINEABOVE',(0,8),(-1,8),1,colors.black),
						('FONTSIZE', (0,0),(-1,7), 8),
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
	line = MCLine(width*0.75, grey)
	Story.append(line)

	# ------
	
	Story.append(Spacer(1, 20))
	ptext = '<font size=9>Target Star Properties</font>'
	Story.append(Paragraph(ptext, styles["Normal"]))
	Story.append(Spacer(1, 12))
	
	Story.append(table_stellar)
	Story.append(Spacer(1, 15))
	Stellartable_text = "Table 1. Stellar properties of the target. * List of the sectors in which the target will be/ has been observed."
	ptext = '<font size=8>%s</font>' % Stellartable_text
	Story.append(Paragraph(ptext, styles["Normal"]))


	
	# --------------------------------------------
	# tess stars 
	# --------------------------------------------
	
	Story.append(Spacer(1, 12))
	im6 = Image(tess_stars_name)
	
	im6._restrictSize(width*0.3, width*0.3)

	Story.append(im6)
	fig_count += 1
	tess_stars_text = "Fig {}. The locations of stars with mag < 15 (orange circle) within the Tess Cut Out around TIC {} (red star). Only shown for one Sector. Data from Gaia DR2.".format(fig_count, tic)
	
	ptext = '<font size=8>%s</font>' % tess_stars_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	

	# --------------------------------------------
	# Background
	# --------------------------------------------
	Story.append(Spacer(1, 12))
	im2 = Image(background_flux_name)
	
	if len(peak_list) == 1:
		im2._restrictSize(width*0.5, width*0.5)
	
	else:
		im2._restrictSize(width*0.9, width*0.9)
	
	Story.append(im2)

	fig_count += 1
	Story.append(Spacer(1, 10))
	background_text = "Fig {}. Background flux vs. time around the time of each transit-like event. \
		The vertical orange line indicates the time of the transit-like event.".format(fig_count)
	
	ptext = '<font size=8>%s</font>' % background_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	

	# --------------------------------------------
	# Centroid Position
	# --------------------------------------------
	Story.append(Spacer(1, 12))
	im3 = Image(centroid_positions_name)
	print (len(peak_list))
	if len(peak_list) == 1:
		im3._restrictSize(width*0.5, width*0.5)
	else:
		im3._restrictSize(width*0.8, width*0.8)

	Story.append(im3)
	
	fig_count += 1
	centroid_text = "Fig {}. The x and y centroid positions around the time of each transit-like event. The black points shows the CCD column and row position of the target’s flux-weighted centroid. \
		The red shows the CCD column and row local motion differential velocity aberration (DVA), pointing drift, and thermal effects. \
		The vertical orange line indicates the time of the transit-like event".format(fig_count)
	
	ptext = '<font size=8>%s</font>' % centroid_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	

	# --------------------------------------------
	# Flux Aperture
	# --------------------------------------------
	Story.append(Spacer(1, 12))
	im4 = Image(flux_aperture_name)
	print (peak_list)
	if len(peak_list) == 1:
		im4._restrictSize(width*0.5, width*0.5)
	else:
		im4._restrictSize(width*0.9, width*0.9)
	Story.append(im4)
	
	fig_count += 1
	Story.append(Spacer(1, 10))
	flux_aperture_text = "Fig {}. The lighcurve around the time of each transit-like event extracted with the SPOC pipeline \
		defined aperture (binned:blue, unbinned:grey) and the with an aperture that is 50% smaller (red). The flux is extracted \
		from the target pixel files (TPFs) and has not been detrended or \
		corrected for systematics. The vertical orange line indicates the time of the transit-like event.".format(fig_count)
	
	ptext = '<font size=8>%s</font>' % flux_aperture_text
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
	flux_aperture_text = "Fig {}. Diffrence images for target TIC {} for each transit like event. \
	Left: mean in-transit flux(left). Right: mean out-of-transit flux. Right: difference between the mean out-of-transit and mean in-transit flux.".format(fig_count, tic)
	
	ptext = '<font size=8>%s</font>' % flux_aperture_text
	Story.append(Paragraph(ptext, styles["Normal"]))

	# --------------------------------------------
	# nearest neighbours
	# --------------------------------------------
	
	Story.append(Spacer(1, 12))
	im7 = Image(nearest_neighbour_name)
	

	im7._restrictSize(width*0.8, width*0.8)

	Story.append(im7)
	fig_count += 1
	nn_text = "Fig {}. Lightcurves of the six closest stars to target {}".format(fig_count,tic)
	
	ptext = '<font size=8>%s</font>' % nn_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	
	
	# --------------------------------------------
	# pixel_LCs_name
	# --------------------------------------------
	
	Story.append(Spacer(1, 12))
	im8 = Image(pixel_LCs_name)
	

	im8._restrictSize(width*0.8, width*0.8)

	Story.append(im8)
	fig_count += 1
	pixLC_text = "Fig {}. Ligthcurve of each pixel. The red shows the SPOC pipeline mask.".format(fig_count)
	
	ptext = '<font size=8>%s</font>' % pixLC_text
	Story.append(Paragraph(ptext, styles["Normal"]))
	

	if bls == True:
		
		try:
			Story.append(Spacer(1, 12))
			im9 = Image(bls1)
			
			im9._restrictSize(width*0.4, width*0.4)
			
			Story.append(im9)
			fig_count += 1
			bls1_text = "Fig {}. The initial BLS excecuted on the whole light curve binned to 10 minutes. Top panel: log liklihood periodogram. \
							The solid red line indicates the peak period and the dashed orange lines show the integer \
							harmonics of this period. Middle panel: Full light curve, unbinned (orange) and binned to 10 minutes (black). \
							The peak period is highlighted by the solid red line. Bottom Panel: Phase folded light curve where the found transit-event is fit \
							with a simple box (red line).".format(fig_count)
			
			ptext = '<font size=8>%s</font>' % bls1_text
			Story.append(Paragraph(ptext, styles["Normal"]))
		
		
			Story.append(Spacer(1, 12))
			im10 = Image(bls2)
			
			im10._restrictSize(width*0.4, width*0.4)
			
			Story.append(im9)
			
			bls2_text = "Fig {}. Second iteration of the BLS where the highest detected signal-to-noise transits have been removed. Top panel: log liklihood periodogram. \
							The solid light blue line indicates the peak period and the dashed blue lines show the integer \
							harmonics of this period. Middle panel: Full light curve, unbinned (blue) and binned to 10 minutes (black). \
							The peak period is highlighted by the solid light blue line. Bottom Panel: Phase folded light curve where the found transit-event is fit \
							with a simple box (light blue line).".format(fig_count)
			
			ptext = '<font size=8>%s</font>' % bls1_text
			Story.append(Paragraph(ptext, styles["Normal"]))
	
			# make a BLS table
			# --------
	
			data_bls= [['Parameter',  				"bls1",                                                        "bls2"],
						   ['depth',	 		"{:.5f} ± {:.5f}".format(bls_stats1[0][0],bls_stats1[0][1]),     "{:.5f} ± {:.5f}".format(bls_stats2[0][0],bls_stats2[0][1]) ],
						   ['depth phased',		"{:.5f} ± {:.5f}".format(bls_stats1[1][0],bls_stats1[1][1]),     "{:.5f} ± {:.5f}".format(bls_stats2[1][0],bls_stats2[1][1]) ],
						   ['depth half',		"{:.5f} ± {:.5f}".format(bls_stats1[2][0],bls_stats1[2][1]),     "{:.5f} ± {:.5f}".format(bls_stats2[2][0],bls_stats2[2][1]) ],
						   ['depth off',	    "{:.5f} ± {:.5f}".format(bls_stats1[3][0],bls_stats1[3][1]),     "{:.5f} ± {:.5f}".format(bls_stats2[3][0],bls_stats2[3][1]) ],
						   ['depth even',       "{:.5f} ± {:.5f}".format(bls_stats1[4][0],bls_stats1[4][1]),     "{:.5f} ± {:.5f}".format(bls_stats2[4][0],bls_stats2[4][1]) ],
						   ['',''],
						   ]
			
			table_bls=Table(data_bls)
			table_bls=Table(data_bls,style=[
								('LINEABOVE',(0,1),(-1,1),1,colors.black),
								('LINEABOVE',(0,6),(-1,7),1,colors.black),
								('FONTSIZE', (0,0),(-1,6), 8),
								])
			
			
			# ------
			
			Story.append(Spacer(1, 12))
			ptext = '<font size=10>Stellar Parameters</font>'
			Story.append(Paragraph(ptext, styles["Normal"]))
			Story.append(Spacer(1, 12))
			
			Story.append(Spacer(1, 15))
			Story.append(table_bls)
	
		except:
			print ("We tried to run the BLS but it didn't work. You'll have to look into this - try running it in a Jupyter notebook to find out why it didnt' work.")


	if model == True:
		
		# --------------------------------------------
		# Pyaneti Modeling results
		# --------------------------------------------
		
		Story.append(Spacer(1, 12))
		im_model = Image(model_name)
			
		im_model._restrictSize(width*0.5, width*0.5)
			
		Story.append(im_model)
		fig_count += 1
		model_text = "Fig {}. Ligthcurve of each pixel. The red shows the SPOC pipeline mask.".format(fig_count)
			
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
		
		    ptext = '<font size=10>Candidate Parameters</font>'
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
		            param_errs.append("(+{:.4f}  -{:.4f})".format(s['Pos'].values[0], s['Neg'].values[0]))
		            param_units.append("{}".format(s['Unit'].values[0] ))
		        except:
		            param_vals.append(-999)
		            param_errs.append(-999)
		            param_units.append(-999)
		
		    data_params = [['Parameters', 'Value']]
		    for p,v,e,u in zip(params, param_vals, param_errs, param_units):
		        data_params.append([p,v,e,u])
		
		    table=Table(data_params)
		    table=Table(data_params,style=[
		                        ('LINEABOVE',(0,1),(-1,1),1,colors.black),
		                        ('LINEABOVE',(0,8),(-1,8),1,colors.black),
		                        ('LINEABOVE',(0,len(params)+1),(-1,len(params)+1),1,colors.black),
		                        ('FONTSIZE', (0,0),(-1,len(params)), 8),
		                        ])
		
		    Story.append(table)
		

		#except: 
		#	print ("Something didn't work with Pyaneti - maye try adjusting the priors.")
	


	doc.build(Story, onFirstPage=addPageNumber, onLaterPages=addPageNumber)



LATTE_DV(tic, indir, peak_list, sectors_all,target_ra, target_dec, tessmag, teff, srad, bls_stats1, bls_stats2, bls = False, model = True)