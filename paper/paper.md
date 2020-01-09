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

date: 1 January 2020
bibliography: paper.bib

---

# Summary


NASA's Transiting Exoplanet Survey Satellite (*TESS*, [@ricker15]) is the first nearly all-sky space-based transit search mission which is monitoring over 80 per cent of the sky, split up into 26 observational sectors, during the nominal two year primary mission alone. TESS's objective is to discover transiting exoplanets around 200,000 pre-selected bright nearby stars which are being monitored every two minutes. Furthermore, TESS is obtaining Full Frame Images (FFIs) at a thirty minute cadence, providing brightness measurements of millions of additional field stars. 

Unlike its predecessor *Kepler*, *TESS* uses very large (20 arcsecond) pixels, meaning that astrophysical false positives (e.g. blends, where the photometric aperture of a bright target is contaminated by a faint eclipsing binary) are expected to be very common. Many false positive scenarios can be ruled out from a detailed examination of the TESS data, (lightcurves and target pixel files), using standard diagnostic tests. 

``LATTE`` is an open source Python package that performs standard vetting tests in order to help weed out astronomical and astrophysical false positive scenarios. The implemented tests are similar to the ‘Data Validation’ steps performed by of the Science Processing Operations Center (SPOC; @jenkins16) pipeline, however, the SPOC pipeline tests are only carried on a selection of TESS targets (Threshold Crossing Events, TCE) and, to date, not on FFI targets. ``LATTE`` allows for the analysis of any *TESS* target with a valid TIC ID [@Stassun19], which include both the pre-selected 2 minute cadence and 30 minute cadence targets.

These standard vetting test can provide a useful tool for astronomical researchers and students. However, ``LATTE`` was also developed with the growing exoplanet citizen science community [e.g., @fischer12; @Christiansen2018; @eisner19] in mind. All interactions with the code are done via a user-interface and the generation of the standardexoplanet  diagnostic tests requires no prior knowledge of Python coding. Upon execution, ``LATTE`` prompts the user enter a TIC ID [@Stassun19] of a target star, as well as the observational sector(s) to analyse (upon being prompted in which Sectors the target was observed) and the location of the transit-like event, which can be identified by zooming in on the plot using the mouse button. All further options, such as whether to create a result summary report, are indicated using 'tick-boxes' in the Graphical User Interface (GUI). Once this information has been entered, the script downloads and handles the data without further user-interference in order to provide essential exoplanet diagnostic plots and figures. 

The default setting for ``LATTE`` is to download and analyse the 2-minute cadence pre-selected data where the lightcurves were extracted and detrended by the SPOC pipeline [SPOC, @jenkins16]. However, run the program in the FFI mode (specified in the command line when the program is executed). In FFI mode the data is downloaded using TESScut [@brasseur2019] and the data detrended using the ``Astropy`` PCA module [@astropy] and an interactive non-linear filter [@Aigrain04] to estimate and subtract residual systematics. Unlike the TESS 2-minute cadence targets, the FFIs do not come with a pre-identified optimal aperture. By default, the user is prompted to manually select the pixels to be used for both a 'large' and a 'small' photometric extraction aperture by clicking on the desired pixels in the user interface. The two lightcurves (from the 'large' and 'small' apertures) are simultaneously displayed. Alternative, the aperture size can be defined automatically by the program using a threshold flux value and centred on the coordinated of the target star.



# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

I thank the all of the Planet Hunters TESS volunteers who prompted me to write a program that may be helpful in identifying the signal of true exoplanets in the TESS data and the support from the entire Zooniverse team and Oxford exoplanet group.

# References


