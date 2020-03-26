# LATTE

Written by Nora L. Eisner

email: *nora.eisner@new.ox.ac.uk*

### THE CODE

--------

*The aim of this code is to provide a series of diagnostic tests which are used in order to determine the nature of events found in *TESS* lightcurves.*

LATTE is an open source Python package that performs standard diagnostic tests to help identify, vet and characterise signals in the *TESS* lightcurves in order to weed out instrumental and astrophysical false positives. The program is designed to allow for a fast, in depth analysis of targets that have already been identified as promising candidates by the main *TESS* pipelines or via alternative methods such as citizen science. The code automatically downloads the data products for any chosen TIC ID (short or long cadence TESS data) and produces a number of diagnostic plots (outlined below) that are compiled in a concise report. 

The implemented diagnostic tests are similar to the ‘Data Validation’ steps performed by the Science Processing Operations Center (SPOC)pipeline for the short-cadence data, and by the quick-look pipeline (QLP) for the long-cadence data. However, while the SPOC pipeline tests are performed on a selection of two-minute cadence *TESS* targets only, LATTE allows for the analysis of any *TESS* target with a valid TIC ID, including both the two-minute and thirty-minute cadence targets.

The code is designed to be fully interactive and does not require any prior knowledge of python coding. 

--------
--------

### Installation

LATTE requires python3 to be installed on your computer, which can be download from https://www.python.org/downloads/. Once you have python3, you can download the code directly from the LATTE github page. Alternatively you can install LATTE using pip (https://pypi.org/project/tessLATTE/) via your command line with:

	pip3 install tessLATTE      

In order for LATTE to work you will need to have the right versions of certain modules installed, so downloading it in a virtual environemt. **Note: ensure that the matplotlib version that you are using is v3.2.0rc1 (pip install matplotlib==3.2.0rc1). You will also need a module called TKinter installed. If this is not already installed please use: sudo apt-get install python3-tk.**

The first time that the program is run you will be prompted to enter a path to a file on your computer where all the output data will be stored (this path can be changed later using --new-path). The first time that the code is run it will also have to download the text data files from MAST, which are around 325M. This download will only have to run in full once but may take a couple of minutes to complete. Please note that this step does not download any of the actual TESS data, but only text files that contain curl scipts that are needed to download individual light curves when LATTE is run. This one time bulk download incerases the speed of LATTE by ensuring that these curl scipts do not have be be downloaded everytime that the code is executed.

If the code doesn't run when you first install it, this could be due to an SSL issue. To check, please try executing the following command from your command line: **python -m tess_stars2px -t 261136679**.

If this command returns the error: *"ssl.SSLError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed"*, you will need to change your SSL settings, which can be done following these [instructions](https://timonweb.com/tutorials/fixing-certificate_verify_failed-error-when-trying-requests_html-out-on-mac/). 


### How to run it?

LATTE is simply run through the command line with:

	python3 -m LATTE      or     python3 -m LATTE --new-data     (if run for the first time or when new data is released)

This will open up a box for you that prompts you to enter a TIC ID and indicate whether you would like to extract the information from the 2-minute cadence ('Standard mode') or 30-minute candence Full Frame Image (FFI) data ('FFI mode').

Once a TIC ID has been entered, the program will tell you in what sectors that target has been observed. If you want to look at all of the sectors, either enter them or simply press enter with nothing in the box. Alternatively, enter the sectors that you are interested in and enter them separated by commas. Remember that LATTE will have to download all the data for each sector so you might not always want to look at all of the sectors. 

**Normal Mode**

The '*normal mode*' looks at the short-cadence *TESS* data which has already been detrended and extracted by the SPOC pipeline. Optimal aperture lightcurve extraction aperture sizes have therefore already been identified and do not need to be selected manually.

**FFI Mode**

In *FFI mode* the data is downloaded using TESScut and the data detrended using PCA, a moving average filter and k-sigma clipping. Unlike the *TESS* 2-minute cadence targets, the FFIs do not come with a pre-chosen optimal aperture. By default, the user is given the option to manually select the pixels within each aperture for a 'large' and a 'small' aperture. This GUI opens automatically and the two apertures are selected by clicking on the desired pixels in the interface (see image below). The two (large and small aperture) lightcurves are simultaneously displayed. When the FFI mode is run with the command '--auto', the aperture size is chosen manually using a threshold flux value centred at the midpoint of the postage stamp.

Please see the LATTE github page for an example of what the interface looks like.

## Input target list

In order to efficiently generate diagnostic plots for multiple targets without having to interactively enter the target TIC IDs, LATTE can be executed with an input file that lists the TIC IDs, the times of the transit-like events, and the observational sectors to consider (if known). See example. In the cases where the times of the transit events or the sectors have not been entered, the user is prompted to enter these manually via the GUI, as descibed above. For longer target lists the code can also be parallelized (see example). If you want to run it with this mode, enter: 

	python3 -m LATTE --targetlist=path/to/the/csv/input/file

The program will skip targets that have already been analysed by you, so if you want to overwrite data add '--o' to the command.

### Transit time selection

Once you have identified the TIC ID, the observational sectors and the aperture sizes (*FFI mode* only), you will see a screen that has the full lighcurve as well as a zoom in of the lightcurve (see figure below). The solid red line across the full lightcurve lets you know where on the lighcurve you are zooming in on. Click on the top (full) or bottom (zoomed in) plots to change the location of the zoom in until the red vertical line is centred on the mid-point of a transit-event. When you are happy with the location press the 'Add Time' button below the plots in order to record this time (in TBJD). You can delete wrongly entered times with the 'Delete time'. The saved times will be shown on the screen. The position of the red line can also be changed by dragging the teal coloured 'Transit' slider with your mouse. The y-scale of the plots can be changed with the grey coloured slider.

Please see the LATTE github page for an example of what the interface looks like.

Additional options are displayed to the left of the plot.

Sector Selection: 
- The sectors that you have chosen to look at are listed in the 'Sector' bar below the plot. You can change the sector to look at by clicking on the corresponding number in the bar, or on 'all' too look at all of them. Alternatively you can use the arrows to the right of the bar to scroll through them. Only the top panel will change. 

Binning Factor: 
- change the binning of the top plot (only available in the normal and not FFI mode)

Settings:
- Simple: only run the most basica diagnostic tests (not suing TPF). This is for a very quick analysis. 
- Show Plots: Display the plots after they are made. This slows down the code but allows you to see and analyse them as they are being made.
- North: Align all the images due North (this slows down the code).
- BLS: Run a Box-Least-Squares algorithm that searches for periodic signals in the lightcurve. 
- Save: Save all the images (default this is already checked)
- Report: Generate a pdf file at the end of the code to summaise the findings (default this is already checked).

Finally, there is an optional box to enter a 'Memorable name' for the candidate. This name is used to store the data at the end in order to make identifying certain targets easier.

Once all the options have been chosen and the transit times stored (you have to enter at least one transit time), press the orange 'Done' button to continue.

The code will then generate download and process all of the data. Note that all the data is downloaded but not stored. The more sectors you analyse in one go the longer the code will take to run.


### Output

For more information on how to interpret the results shown in the reports look at the 'How to make a LATTE' document in the 'example' folder. Please see the LATTE github page for example plots.

**Figures:**

- Full lightcurve with the times of the momentum dumps marked. 


The solid red lines at the bottom of the figure indicated the times of the reaction wheel momentum dumps and the dashed black line(s) show the times of the marked transit events. Momentum dumps occur around every 2 to 2.5 days and typically last around half an hour. This can affect the satellite pointing, resulting in spurious transit-like signals. 
	

- Background flux around the times of the marked transit event(s).

Enhanced scattered light in the telescope optics can cause the background flux to rise sharply each time TESS approaches the perigee of its eccentric orbit around the Earth. Additionally, solar system objects, such as asteroids passing through the field of view will show up as a spike in the background flux. Both of these can result in a transit-like signal in the light curve. 


- Centroid positions around the time of the marked transit event(s).

The black points shows the CCD column and row position of the target’s flux-weighted centroid. The red shows the CCD column and row local motion differential velocity aberration (DVA), pointing drift, and thermal effects. Significant spikes in either of these at the time of the transit-like event suggest that the signal is not due to a transiting planet. 


- The lightcurve around the time of the marked event(s) extracted in two different aperture sizes. 

For a transiting planet, the transit shape and depth are expected to remain constant when extracted with different aperture sizes. Conversely, they are expected to change if the signal is due to a blended eclipsing binary. This is because for blended binaries, the signal is not centred on the target and instead originates from elsewhere in the field of view, thus changing the shape and depth of the signal with different apertures. I caution that if the selected aperture is not centred on the target, this test can lead to the mis-interpretation of a moving signal. I, therefire, recomend to always check the figure showing the used apertures before interpreting these results. These data are extracted from the Target Pixel Files (TPFs) and are neither detrended nor corrected for systematics.


- The apertures used for the extraction. Please note that for very bright and very faint stars the aperture selection for the smaller aperture may not work correctly so these should be checked in order to correctly interperet the results.

- The average flux in and out of the marked event(s) and the differences between the two.

A plot showing the average in-transit and average out of-transit flux as well as the difference between the two. A change in spatial distribution of the flux suggests that the transit-like event is caused by a blend with a background eclipsing binary. The data are extracted from the TPFs and corrected for systematic effects using Principal Component Analysis. 

- The average flux of the target pixel file with the locations of nearby stars (magnitude < 15) indicated (GAIA DR2 queried).

A plot showing the location of nearby stars brighter than *TESS* magnitude 17 as queried from the Gaia Data Release 2 catalog and the DSS2 red field of view around the target star. The sizes of the markers indicate the positions of the Gaia stars relate to the magnitude of each companion. The *Astropy reproject* module is used to align the plots such that north points towards the top of the page.

- The lightcurves of the 5 closest stars that were also observed by *TESS* (TP files).

Lightcurves of the five two-minute cadence stars closest to the selected target. The occurrence of similar transit-like events in nearby lightcurves may suggest that the signal is caused by a background event or a blended eclipsing binary.

- A lightcurve around the time of the marked event(s) extracted for every pixel in the target pixel file.

A plot showing lightcurves extracted for each individual pixel around the target. The data are extracted from the TPFs and detrended using a third order spline fit. If the transit-like signal is more prominent in pixels that are not centred on the target star, this suggests that the signal is caused by a blended background eclipsing binary. 

- (optional) Two simple BLS plots that look for periodic signals in the data. The second with the highest detected signal-to-noise transits from the initial run removed.
- (in progress, will be available in next release of LATTE) Modelling of the signal using a Bayesian approach with an MCMC sampling. This makes use of the Pyaneti code (Barragan et al. 2017). 

**FFI Mode**

- The FFI mode currently does not plot the nearby stars lightcurves (will be implemented soon).
- Saves the extracted apertures used
- Saves the corrected and the uncorrected lightcurves to verify that the detrending is workign correctly - there's are nt stored in the DV reports. 


**Tables:**

- Stellar parameters summarized, as well as information to whether the target has been flagged as a TCE or a TOI. The table links to the relevant reports (if applicable) as well as to the exofop page on the star.
- (optional) Summary of the BLS results. 
- (optional) Fitting parameters with uncertainties from the modelling. 

### Arguments

NOTE: all of these arguments (except new-path, auto and targetlist, new-data) can be changed as option in the GUI. They are arguments in case the same options wish to be run multiple times and the user therefore wishes to identify them in the command line when the program is executed.

**--new-data**  The code requires multiple text files to be stored on your computer in order to run - these are downloaded automatically from the MAST server. The first time the proghram is run, and any time that there is new data available, add **--new-data** to the command line when running the program. The code checks what data has already been downloaded and doesn't re-download anything that already exists.

**--auto** When looking at the FFIs, the default option is that you choose both the large and small apertures interactivelty. In order for the system to choose them run the command with '--auto'.

**--targetlist***=path/to/the/csv/input/file* instead of entering the information manually everytime, you can provide LATTE with an input file which lists all the TIC IDs and the times of the transit events. Look at the example file to see the required format of the information.

**--new-path** If you want to define a new path to store the data.

**--o** If the input target list option is chosen, the program will check whether each target has already been analysed, in which case it will skip this target. If you do not wish to skip targets that have already been assessed use this in order to specify the 'overwrite' option. When the program is run interactively (not with an input file) this option has no effect.

**--noshow** if you do not want the figures to pop up on your screen you can run the program with this command in the command line. The program will run significantly faster with this run. If this option is chosen the files are always saved (also option in the GUI)

**--nickname** In order to keep track of all of the candidates, it can be useful to asign them a nickname. This can be entered here which will simply change the name of the folder at the end (also option in the GUI)

**--FFI** If you want to look at an FFI write '--FFI' in the command line (also option in the GUI) 

**--north** If you want all the images to be oriented due north (this is not the default as it takes longer to run)

**--mpi** If the code is parallelized (see example), this needs to be entered in the command line. This is because the module that reprojects the TPFs, astroplan, cannot be parallelized due to problems with multiple processes accessing python's shelve storage at the same time.

**--tic** You can skip the box to enter the tic ID by entering it in the command line with e.g. --tic=55525572. 

**--sector** You can skip entering the sectors by entering them in the command line with e.g. --sector=2,5. You will need to know in what sectors this target was observed (option in the GUI)


### Example

Please see the LATTE github page for a short video that shows how to run LATTE for example target TIC 94986319 (TOI 421). The the output files that this example produces can be found in the folder called 'example_output'. The summary DV report, which can be found in this folder, can be used to verify that this these atrget stars pass all of the initial vetting stages and inform decisions as to whether to follow it up using ground based telescopes. 


If you want to set or change the input/output path, alter the terminal command to: 

	python -m LATTE --new-path

and if you want to download new data (or check for new data) use the command: 

	python -m LATTE --new-data


#### MPI example

If you have a computer with multiple cores you can generate multiple reports in one go by paralizing the code. There is an example scipt for how this can be done on the github page, although I caution that the code may have to be adapted to work with your operating system. 

The MPI code requires an input file (see 'example_input_file.csv'). In order to keep track of the targets that have been processed, the code generates a manifest file listing some of the information of the processed targets such as the transit-event times, the RA, Dec, sectors, TMag, Teff, whether it's a TOI and the link to the TCE report if it's available. 

An additional input parameter for this code is --number. This sets the number of targets that you want to analyse from the list. For example, if you have a list of 1000 targets but you only want to look at the first 10, you can set --number=10. 

On a mac, the example parallelized code can be executed using: 

	mpiexec -np 2 python example_parallelized_code.py --targetlist=/path/to/input/code/example/example_input_file.csv --number=10  --mpi

It's import that you add the --mpi, otherwise the code cannot be parallelized due to memory conflicts. 


### LATTE workflow - what is being downloaded, when and why?

Below is an outline of how the code works and when it downloads the various different data files. 

1) When the code is first exectuted, a number of text files are downloaded which are required for the code run run more quickly for all subsequent runs. These include: 

	a) Two text files for each TESS observational sector, one for the short cadence light curves, and one for the target pixel files, which contain curl commands that link to the download of each individual TESS light curve. These are simply text files and no actual TESS data and will only have to be downlaoded once for each sector. This is in order to increase the speed of the code later on. These files are around 3.9 MB each.

	b) One text file for each sector containing the RA, Dec and TESS magnitude for each target (760 KB per file). This information is used to identify tge TIC IDs of the 5 closests TESS short-cadence targets.

	c) One file that lists the times of the momentum dumps across all of the sectors (23 KB). In order to generate this file one lightcurve is downloaded from each sector and the data extracted. These lightcurves are openend but never saved. This is only needed when looking at teh FFIs.

	d) A list of all of the confirmed TOIs (274 KB) and all of the targets that were flagged as TCEs by the SPOC pipeline (20.2 MB). 

2) When you run LATTE, you will first be prompted to enter the TIC ID. Once the TIC ID has been entered, LATTE calls a program known as 'tess-point', which tells us in which TESS osbervational sectors this target was osberved in. This information is stored in a temporary text file so that it can be accessed by the code later on. 

3) You are then prompted to enter the observational sectors that you want to look at in more detail. At this point LATTE opens the short cadence text files for each observational sector that you entered, and uses the curl scipts for the relevant TIC ID to download a short-cadence light curve (2MB per light curve). The more sectors that you stated to look at, the longer this step will take. If you chose the 'simple' option within the TESS GUI, this will be the only data that is downloaded. 

4) Once you have identified the times of the transit like events, the code will generate a plot showing the background flux, and one showing the mean centroid positions. These plots use the alreadt downloaded lightcurve file. 

5) Next, the code will download the target pixel file using the lightkurve package (unless you run the 'simple' mode). This file is downlaoded once the centroid plot has been closed and used to extract lighcurves with two different aperture sizes. Due to the size of the file and the time taken to determine the optimal aperture size, this step may take a couple of seconds to complete. 

6) The next figure shows the brightests stars surrounding the target, both on the average TPF data and on an SDSS image. The target pixel files are downloaded directly from the MAST server at this stage (48 MB per light curve), using the previously downlaoded curl scripts. This file is larger than the light curve  data and will therfore take longer to download. This file will be used for multiple plots. If the 'north' option has been selected, the image will be reprojected such that both images are orentied such that North is directed upwards. This option takes longer.

7) The next figure compares the flux in-transit to out-of-transit for each marked transit-like event. This uses the already downlaoded TPF and should therefore not take long to load. 

8) The next figure shows a lightcurve extracted for each pixel surrounding the target around the time of the marked transit-like event using the already downlaoded TPF. . This is done individually for each amrked transit and will therefore take longer if more transit have been identified. Your terminal window will indicate how much longer this step takes to run. 

9) If more than one transit-like event has been identified, the lightcurve will be phase folded. This uses the previously downlaoded lightcurve file and should be completed very quickly. 

10) The next figure shows the lightcurves of the five closest short-cadence TESS stars to the target. The data for each star is downloaded using the curl scipts. Your terminal window will indicate once each lighcurve has been downloaded sucessfully. 

11) The final step is to compile and save the LATTE data validation report. 


### Overview of the main LATTE scipts

__main__.py      : Intitialises the parameters: TIC ID, Sector, FFI or not, run with input file or not, check that the data has been downloaded etc. In brief, this code is designed check that all of the input parameters are consistent with one another and informs the rest of the code how to run. 

LATTEutils.py    : All the functions needed to download data and text files, function to run the interactive GUI and all of the plotting and data handling. Functions in this scipt are called both from 'main' and from 'brew'.

LATTE_DV.py      : Scipt to combine all of the results from LATTEbrew in order to generate a pdf data validation report.

LATTEbrew.py     : Scipt that calls the LATTEutils functions in turn in order to store the results. Calls the LATTE_DV.py function to collate all the results at the end.


### Perform automated tests

The code has been tested using the python unittest module. There are a number of tests including testing the connections to the servers where data is downloaded as well as tests that process pre-downloaded data in order to ensure that the outputs are as expected. In order to execute the tests please clone the GitHub LATTE repository. You must then be in the parent LATTE directory and from there you can run all the tests using: 

	python3 -m unittest

or individual tests using, for example: 
	
	python3 -m unittest LATTE/tests/test_DVreport.py

Before running the tests makes sure that you have configured an output path (*python3 -m LATTE --new-path*) and that the data files are downloaded (*python3 -m LATTE --new-data*). 

### Contributing

If you believe that you have found any bugs in the code or if you need help or have any questions or suggestions, please feel free to file a new Github issue. You can also raise an issue or a pull request which fixes the issue. Please see CONTRIBUTING.md for more details. 


### Future Work

The next release will allow for the option to model the transit-like events using the open source package *Pyaneti* [@pyaneti] which uses a Bayesian approach with an MCMC sampling to determine the parameters of the best-fit.

Future releases will also allow for modeling and removal of periodic trends in the light curve. This is useful for the detection and circumbinary planets and for the analysis of multi-star systems. 


### License 

LATTE can be used under GNU Lesser General Public License v3.0


### Acknowledgements

I thank all of the Planet Hunters *TESS* volunteers whose dedication to the project encouraged me to write this analysis tool. I am also extremely grateful to all of the support provided by the Zooniverse team and the Oxford exoplanet group.














