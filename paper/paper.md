---
title: 'LATTE: Lightcurve Analysis Tool for Transiting Exoplanets'

tags:
  - Python
  - astronomy
  - TESS
  - exoplanets
  - transit

authors:
  - name: Nora L. Eisner
    orcid: 0000-0002-9138-9028
    affiliation: 1

  - name: Chris J. Lintott
    affiliation: 1

  - name: Suzanne Aigrain
    affiliation: 1


affiliations:
 - name: Oxford Astrophysics, University of Oxford, Denys Wilkinson Building, Keble Rd, Oxford OX1 3RH, UK
   index: 1

date: 9 January 2020
bibliography: paper.bib

---

# Summary

NASA's Transiting Exoplanet Survey Satellite [*TESS*, @ricker15] is the first nearly all-sky space-based exoplanet search mission. During the nominal two year primary mission alone, *TESS* will monitor over 80 per cent of the sky, split up into 26 observational sectors. The prime objective of the mission is to discover transiting exoplanets around 200,000 pre-selected bright nearby stars which are being monitored every two minutes. In addition to these short-cadence observations, *TESS* is obtaining Full Frame Images (FFIs) every thirty minutes.

Unlike its predecessor *Kepler*, *TESS* uses very large (20 arcsecond) pixels, meaning that astrophysical false positives (e.g. blends, where the signal of a target is contaminated by a faint eclipsing binary) are expected to be common. Many false positive scenarios can be ruled out from detailed examination of the *TESS* data (lightcurves and target pixel files) using standard diagnostic tests. 

``LATTE`` is an open source Python package that performs these standard diagnostic tests to help identify, vet and characterise signals in the *TESS* lightcurves in order to weed out instumental and astrophysical false positives. The program is designed to allow for a fast, in depth analysis of targets that have already been identified as promising candidates by the main *TESS* pipelines or via alternative methods such as citizen science [e.g., @fischer12; @Christiansen2018; @eisner19]. The implemented diagnostic tests are similar to the ‘Data Validation’ steps performed by the Science Processing Operations Center [SPOC; @jenkins16] pipeline. However, while the SPOC pipeline tests are performed on a selection of two-minute cadence *TESS* targets only, ``LATTE`` allows for the analysis of any *TESS* target with a valid TIC ID [@Stassun19], including both the two-minute and thirty-minute cadence targets.


# ``LATTE`` 

The standard vetting tests performed by ``LATTE`` provide a useful tool for both researchers and students. Furthermore, ``LATTE`` was developed with the growing exoplanet citizen science community [e.g., @fischer12; @Christiansen2018; @eisner19] in mind and thus all interactions with the code are done via a user-friendly interface, requiring minimal experience in python programming.

Upon execution, ``LATTE`` prompts the user enter the TIC ID [@Stassun19] of a target star, as well as the observational sector(s) to analyse. The corresponding data are then downloaded directly from the Mikulski Archive for Space Telescopes (MAST, http://archive.stsci.edu/tess/) and displayed in a Graphical User Interface (GUI), as shown in Figure 1. The user is asked to identify the times of the transit-like events by directly clicking on the lightcurve plots and recording the transit-times using the 'Add time' button. Additional settings, such as whether to run a Box-Least-Squares (BLS) algorithm that searches for periodic signals in the lightcurve, can be selected within the GUI (see left side of Figure 1). Once the times of one or more transit-like events have been identified the program automatically downloads and handles the *TESS* data to provide standard exoplanet diagnostic plots and statistics.

![User-interface displaying the whole lightcurve (top pannel) and a zoom-in of the data around the selected transit-time (bottom pannel) as indicated by the vertical red line. The selected time can be changed by clicking on either one of the two lightcurve plots. Times are stored and removed with the 'Add time' and 'Remove time' buttons respectively. Further program options are shown to the left of the lightcurve plots: '*Simple*' for a quick analysis that doesn't require the Target Pixel Files (TPFs); '*Hide plots*' to prevent the figures from being displayed as they are being produced; '*North*' to reproject the TPF plots due north; '*BLS*' to run a Box-Least-Squares algorithm to search for periodic signals; '*Save*' to store the files locally (selected by default); and '*Report*' to create a PDF summary report (selected by default). There is also an option to change the binning of the displayed full lighcurve.](Fig1.png)


## Different Modes

### Standard mode 

By default ``LATTE`` runs in *standard mode* using the pre-selected short-cadence *TESS* data and lightcurves that were extracted and detrended by the Science Processing Operations Center [SPOC, @jenkins16] pipeline. The *Lightkurve* package [@2018lightkurve] is used to extract raw lightcurves from the Target Pixel Files.

### FFI mode

The *FFI mode*, which can be selected with a check-box at the time when the TIC ID is entered, downloads the FFI *TESS* data using *TESScut* [@brasseur2019]. The data are then detrended using the *Astropy PCA* module [@astropy] and an iterative non-linear filter [@Aigrain04] to estimate and subtract residual systematics. Unlike the two-minute cadence *TESS* targets, the FFI targets have no pre-identified optimal apertures. The user is, therefore, prompted to manually identifying two different sized photometric extraction apertures by interacting with a GUI, which displays the average flux per pixel around the target (see Figure 2). The apertures are chosen by clicking on the desired pixels, where the lighter coloured pixels represent higher flux. The two lightcurves, extracted using the chosen 'large' and 'small' apertures, are simultaneously displayed. Alternatively, the user can opt to have the two apertures selected automatically using a threshold flux value where the apertures are centred on the coordinates of the target star.

![The FFI interface which allows the user to manually select the pixels for two seperate apertures (right), where the lighter coloured pixels correspond to a higher flux. The resultant lightcurves, extracted with the larger (teal green circles) and a smaller (pink crosses) apertures, are displayed on the left.](Fig2.png)

#  Diagnostic tests

Summary of the ``LATTE`` diagnostic tests:

- **Momentum Dumps.** A plot highlighting the times of the *TESS* reaction wheel momentum dumps to ensure that the times of the transit-like events do not coincide with times of increased systematics due to the satellite unloading built up photon pressure.

- **Target Properties Table.** A table that summarized the properties of the target star, including the stellar radius, effective temperature, TESS magnitude, and the galactic coordinates. The table also lists all of the *TESS* sectors in which the target will be osberved. Additionaly, it indicates whether the target has been flagged as a Target Threshold Crossing Event (TCE) or as a Tess Object of Interest (TOI) by the SPOC pipeline, with links to the corresponding validation reports as well as a link to the ExoFOP page that contains additional stellar information [@exofop2019AAS].

- **Background Flux.** A plot to verify that the transit-like dips in the lighcurve are not caused by background events such as asteroids passing through the field of view.

- **x and y centroid positions.**  A plot showing the CCD column and row local position of the target’s flux-weighted centroid, and the CCD column and row motion which considers differential velocity aberration (DVA), pointing drift, and thermal effects.

- **Aperture size test.** A plot showing the lightcurves extracted with two different aperture sizes. In the *FFI mode* the two different aperture sizes can be selected manually (see above), while in the *standard mode* the SPOC pipeline defined aperture, and an aperture that is 50% smaller, are used. Differences in the transit depth and shape of the two different lightcurves may be indicative of the transit-like events resulting from a blended eclipsing binary. However, we caution that a sub-optimal mask selection can lead to the mis-interpretation of a moving signal. The data are extracted from the Target Pixel Files (TPFs) and are neither detrended nor corrected for systematics.

- **Pixel-level centroid analysis.** A plot showing the average in-transit and average out of-transit flux as well as the difference between the two. A change in spatial distribution of the flux suggests that the transit-like event is caused by a blend with a background eclipsing binary. The data are extracted from the TPFs and corrected for systematic effects using Principal Component Analysis [e.g., @Stumpe2012].

- **Nearby companion stars.** A plot showing the location of nearby stars brighter than *TESS* magnitude 17 as queried from the Gaia Data Release 2 catalog [@Gaia2018] and the DSS2 red field of view around the target star. The sizes of the markers that indicate the positions of the Gaia stars relate to the magnitude of each companion. The *Astropy reproject* module [@astropy] was used to align the plots such that north points towards the top of the page.

- **Nearest neighbour lightcurves** (not available in *FFI mode*). Lightcurves of the five two-minute cadence stars closest to the selected target. The occurence of similar transit-like events in nearby lightcurves may suggest that the signal is caused by a background event or a blended eclipsing binary.

- **Pixel level lighcurves.** A plot showing lightcurves extracted for each individual pixel around the target. The data are extracted from the TPFs and detrended using a third order spline fit.

- **Box-Least-Squares fit.** A plot showing the results from two BLS searches. The data are corrected for residual systematics using an iterative non-linear filter [@Aigrain04] prior to being searched for periodic signals using the *Astropy BLS* module [@astropy]. The intitial BLS search identifies the times of the highest detected signal-to-noise event and removes these data from the lightcurve. A second search is then carried out in order to search for additional periodic signals.


All of the results are summarised in a pdf report. 

## Input target list

In order to efficiently generate diagnostic plots for multiple targets without having to interactively enter the target TIC IDs, ``LATTE`` can be executed with an input file that lists the TIC IDs, the times of the transit-like events, and the observational sectors to consider (if known). In the cases where the times of the transit events or the sectors have not been entered, the user is prompted to enter these manually via the GUI, as descibed above. For longer target lists the code can also can be parallelized (see https://github.com/noraeisner/LATTE for an example.)

# Installation and command line execution

``LATTE`` is developed publicly through a GitHub repository, which provides the user with example input and output files, as well as the opportunity to highlight arising issues with the code. The main aspects of the code have been tested using the *unittest* package. ``LATTE`` can be installed either directly from the GitHub repository or via the pip command on PyPi (https://pypi.org/project/tessLATTE/).

Once installed, the program is easily excecuted via the command line (e.g. python3 -m LATTE). Further command line options, such as to run the code with an input target list, can be entered in order execute the code with more advanced settings. For further information on these options please see the documentations provided on the GitHub page. The first time that the code is executed, the user is prompted to enter a local path that directs to a location on the computer where the data and output files will be stored. Additionaly, on the initial run, ``LATTE``  automatically downloads text files that are required for the program run and may, therefore, take a couple of minutes to set-up. 


# Future

The next release will allow for the option to model the transit-like events using the open source package *Pyaneti* [@pyaneti] which uses a Bayesian approach with an MCMC sampling to determine the parameters of the best-fit.


# Acknowledgements

I thank all of the Planet Hunters *TESS* volunteers whose dedication to the project encouraged me to write this analysis tool. I am also extremely grateful to all of the support provided by the Zooniverse team and the Oxford exoplanet group.


# References


