# LATTE

Written by Nora L. Eisner

email: *nora.eisner@new.ox.ac.uk*

### THE CODE

--------

*The aim of this code is to provide a series of diagnostic tests which are used in order to determine the nature of events found in TESS lightcurves.*


--------
--------
How to run it? 


**Normal Mode**

LATTE is simply run through the command line with: (don't forget to run it with --new-data the first time - see below)

			python3 LATTE.py         

This will open up a box for you that prompts you to enter a TIC ID. 

Once a TIC ID has been entered, the program will tell you in what sectors that target has been observed. If you want to look at all of the sectors, either enter them of simply press enter with nothing in the box. Alternatively, enter the sectors that you are interested in and enter them separated by commas. Remember that LATTE will have to download all the data for each sector so you might not always want to look at all of the sectors. 

Next, you will see a screen that has the full lighcurve as well as a zoom in of the lightcurve. The solid red line across the entire full lightcurve lets you know where on the lighcurtver you are zooming in on. Move the 'zoomed bar' by sliding the teal coloured slider across with your mouse. This allows you to locate transit-like events in the data. Once you have found an event that you want to explore further, enter the time of it in the etxt box below (the time is indicated to the right of the slider). You can enter as many times as you like, separated by commas. I note that the slider currently becomes more difficult to navigate when you are looking at multiple sectors. 

Options:
- change the binning of full lightcurve (not the zoomed in one at the moment) on the left.
- change the y scale (flux) with the sliding bar underneath the plot.
- In this widnow you have to indicate whether you want to run a BLS search, model the transit, save the output files (default is yes) and generate summary pdf - Data Validation Report (DVR) - at the end. 

All of the options and the listed times of the events are entered when the orange 'close' button is pressed.

Once the window has been closed, all of the plots will be generated automatically with no further input. If the prgram is not run with --noshow (see the additional arguments section below), all of the figures have to be closed manually before the program continues to run. 


**FFI Mode**

In FFI mode the data is downloaded using TESScut and the data detrended using PCA, a moving average filter and k-sigma clipping. 
Unlike the TESS 2-minute cadence targets, the FFIs do not come with a pre-chosen optimal aperture. By default, the user is given the option to manually select the pixels within each aperture for a 'large' and a 'small' aperture. This window pops up and the apertures are chosen but clicking on the pixels in the interface. The two (large and small aperture) lightcurves are simultaneously displayed. When the FFI mode is run with the comman '--auto', the aperture size is chosen manually using a threshold flux values and places at the centre of the image. 


### Arguments

!!!!!!  **--new-data**     The code needs a couple of texts files stored on your computer in order to run - these are downloaded from the internet automatically for you. The first time you run the program, and any time that there is new data available run add this to the command line when runnign the program. Don't worry, the code won't re download anything that already exists - there are checks to see what is already downloaded and won't download it again.

**--targetlist***=path/to/the/csv/input/file* instead of entering the information manually everytime, you can provide LATTE with an input file which lists all the TIC IDs and the times of the transit events. Look at the example file to see the required format of the information.

**--noshow** if you do not want the figures to pop up on your screen you can run the program with this command in the command line. The program will run significantly faster with this run. If this option is chosen the files are always saved. 

**--o** If the input target list option is chosen, the program will check whether each target has already been analysed, in which case it will skip said target. If you do not wish to skip targets that have already been assessed use this in order to specify the 'overwrite' option. When the program is run interactively (not with an input file) this option has no effect.

**--nickname** In order to keep track of all of the candidates, it can be useful to asign them a nickname. This can be entered here which will simply change the name of the folder at the end. 


**--tic** You can skip the box to enter the tic ID by entering it in the command line with e.g. --tic=55525572. 

**--sector** You can skip entering the sectors by entering them in the command line with e.g. --sector=2,5. You will need to know in what sectors this target was observed.

**--FFI** If you want to look at a FFI write '--FFI' in the command line. 

**--auto** When looking at the FFIs, the default option is that you choose both the large and small apertures interactivelty. In order for the system to choose them run the command with '--auto'. 


### Example

	python3 LATTE.py --nickname=example     


**Figures:**

- Full lightcurve with the times of the momentum dumps marked. 

![Full LC](https://github.com/noraeisner/LATTE/blob/master/example_output/94986319_fullLC_md.pngs=200)

- Background flux around the times of the marked transit event(s).
- Centroid positions around the time of the marked transit event(s).
- The lightcurve around the time of the marked event(s) extracted in two different aperture sizes. 
- The average flux in and out of the marked event(s) and the differences between the two.
- The average flux of the target pixel file with the locatiosn of nearby stars (magnitude < 15) indicated (GAIA DR2 queried).
- The lightcurves of the 6 closest stars that were also observed by TESS (TP files).
- A lightcurve around the time of the marked event(s) extracted for every pixel in the target pixel file.
- (optional) Two simple BLS plots. The second with the highest detected signal-to-noise transits from the initial run removed.
- (optional) Modelling of the signal using a Bayesian approach with an MCMC sampling. This makes use of the Pyaneti code (Barragan et al. 2017). 
- **FFI Mode**The FFI mode currently does not plot the nearby stars lightcurves (will be implemented soon)
- Saves the extracted apertures used
- Saves the corrected and the uncorrected lightcurves to verify that the detrending is workign correctly - there's are nt stored in the DV reports. 


**Tables:**

- Stellar parameters summarized, as well as information to whether the target has been flagged as a TCE or a TOI. The table links to the relevant reports (if applicable) as well as to the exofop page on the star.
- (optional) Summary of the BLS results. 
- (optional) Fitting parameters with uncertainties from the modelling. 













