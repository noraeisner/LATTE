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

affiliations:
 - name: Department of Physics, University of Oxford, Oxford
   index: 1

date: 9 January 2020
bibliography: paper.bib

---

# Summary

NASA's Transiting Exoplanet Survey Satellite (*TESS*, [@ricker15]) is the first nearly all-sky space-based transit search mission. During the nominal two year primary mission alone, TESS will monitor over 80 per cent of the sky, split up into 26 observational sectors. The prime objective of the mission is to discover transiting exoplanets around 200,000 pre-selected bright nearby stars which are being monitored every two minutes. In addition to the 2-minute cadence observations, TESS is obtaining Full Frame Images (FFIs) at a 30-minute cadence.

Unlike its predecessor *Kepler*, *TESS* uses very large (20 arcsecond) pixels, meaning that astrophysical false positives (e.g. blends, where the photometric aperture of a bright target is contaminated by a faint eclipsing binary) are expected to be very common. Many of these false positive scenarios can be ruled out from a detailed examination of the TESS data, (lightcurves and target pixel files), using standard diagnostic tests. 

``LATTE`` is an open source Python package that performs these standard vetting tests in order to help weed out astronomical and astrophysical false positive scenarios. The implemented tests are similar to the ‘Data Validation’ steps performed by the Science Processing Operations Center (SPOC; @jenkins16) pipeline, however, while the SPOC pipeline tests are only carried on a selection of 2-minute cadence TESS targets, ``LATTE`` allows for the analysis of any *TESS* target with a valid TIC ID [@Stassun19], include both the pre-selected 2-minute cadence and 30-minute cadence targets.

# ``LATTE`` 

The standard vetting test performed by ``LATTE`` provide a useful tool for both astronomical researchers and students. Furthermore, ``LATTE`` was developed with the growing exoplanet citizen science community [e.g., @fischer12; @Christiansen2018; @eisner19] in mind and thus all interactions with the code are done via a user-interface, requiring minimal knowledge of python coding.

Upon execution, ``LATTE`` prompts the user enter a TIC ID [@Stassun19] of a target star, as well as the observational sector(s) to analyse. The corresponding data is then downloaded directly from the Mikulski Archive for Space Telescopes (MAST, http://archive.stsci.edu/tess/) and displayed in a Graphical User Interface (GUI) as shown in Figure 1. The user is then asked to identify the times of the transit-like events by clicking on the lightcurve plots and storing the times of the transit like events with the 'Add Time' button. Additional options for the running of the program, such as whether to run a Box-Least-Squares (BLS) algorithm that searches for periodic signals in the lightcurve, can be selected within the GUI (see Figure 1). Once one or more transits have been identified and the ‘Done’ button is pressed, the script downloads and handles the *TESS* data to provide essential exoplanet diagnostic plots and figures. 

![Figure 1.](figure.png)


# Diagnostic Tests

The ``LATTE`` diagnostic test:

- Momentum Dumps. A plot of the times of the satellite momentum dumps allows the user to ensure that the times of the transit-like events do not coincide with times of increased systematics caused by the satellite.

- Background Flux. A plot to verify that the transit-like dips in the lighcurve are not caused by background events such as an asteroid passing through the field of view.
- The x and y centroid positions.  The plot shows both the CCD column and row local position of the target’s flux-weighted centroid and the motion differential velocity aberration (DVA), pointing drift, and thermal effects, allowing the user to identify times of enhanced systematic effects introduced by the satellite.
- Light curve extraction with different aperture sizes. Difference in the transit depth and shape of the two different light curves may be indicative of the transit-like events being the result of a blended eclipsing binary. The data is extracted from the Target Pixel Files (TPFs) and is not detrended or corrected for systematics.
- Pixel-level centroid analysis. This plot shows the average in-transit and average out of-transit flux as well as the difference between the two. A change in spatial distribution of the flux suggests that the change in flux at the time of the transit-like event are caused by a blend with a background eclipsing binary. The data is extracted from the TPFs and are corrected for systematic effects using Principal Component Analysis [e.g. @Stumpe2012].
- Nearby companion stars. A plot showing the location of nearby stars brighter than TESS magnitude 17, queried from the Gaia Data Release 2 catalog [Gaia Collaoration 2018], as well as a the DSS2 red field of view around the target star. The plots are reprojected so that north aligns with the top of the image [@astropy].
- Nearest neighbour lightcurves (not available in FFI mode). Lightcurves of the five closest pre-selected 2-minute cadence stars to the selected target. The manifestation of the transit-like event in other nearby lightcurves would suggest that the event is caused by a background event or a blended eclipsing binary.

- Pixel level lighcurves. A plot showing a lightcurve extracted from each individual pixel. The data is extracted from the TPFs and detrended using a third order spline fit.
- (optional) Box-Least-Squares (BLS) fit. The data is corrected for residual systematics using an iterative non-linear filter [@Aigrain04] prior to being searched for periodic signals using astropy’s BLS function [@astropy]. The initial BLS search identifies the times of the highest detected signal-to-noise event and removes these data from the LC. The search is then carried out again in order to search for additional periodic signals.

# FFI mode

The default setting for ``LATTE`` is to download and analyse the 2-minute cadence pre-selected data where the lightcurves were extracted and detrended by the SPOC pipeline [SPOC, @jenkins16]. The FFI mode, which can be selected with a check-box at the time when the TIC ID is entered, downloads the TESS data using TESScut [@brasseur2019]. The data is then detrended using the ``Astropy`` PCA module [@astropy] and an iterative non-linear filter [@Aigrain04] to estimate and subtract residual systematics. Unlike the TESS 2-minute cadence targets, the FFIs do not come with a pre-identified optimal aperture. By default, the user is prompted to manually select the pixels to be used for both a 'large' and a 'small' photometric extraction aperture by clicking on the desired pixels in the user interface. The two lightcurves (from the 'large' and 'small' apertures) are simultaneously displayed. Alternative, the aperture size can be defined automatically by the program using a threshold flux value and centred on the coordinated of the target star.

![Figure 1.](figure.png)

# Future 

The next release will allow for the option to model the transit-like events using the open source package Pyaneti [@pyaneti] which uses a Bayesian approach with an MCMC sampling.


# Acknowledgements

I thank the all of the Planet Hunters TESS volunteers who prompted me to write a program that may be helpful in identifying the signal of true exoplanets in the TESS data and the support from the entire Zooniverse team and Oxford exoplanet group.

# References


