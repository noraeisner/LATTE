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


NASA's Transiting Exoplanet Survey Satellite (TESS, [@ricker15]) is the first nearly all-sky space-based transit search mission which is monitoring over 80 per cent of the sky, split up into 26 osbervational sectors, during the nominal two year primary mission alone. TESS's objective is to discover transiting exoplanets around 200,000 pre-selected bright nearby stars whcih will be monitored every two minutes. Furthemore, TESS is obtaining Full Frame Images (FFIs) at a thirty minute cadence, providing brightness measurements of millions of additional field stars. 

TESS uses very large (20 arcsecond) pixels, so astrophysical false positives (e.g. blends, where the photometric aperture of a bright target is contaminated by a faint eclipsing binary) are expected to be common. Many false positive scenarios can be ruled out from a detailed examination of the TESS data (lightcurve and target pixel files) using standard diagnostic tests. 

``LATTE`` is an open source Python package that performs standard vetting tests, similar to the ‘Data Validation’ step of the Science Processing Operations Center (SPOC; @jenkins16) pipeline, using the TESS lightcurves and target pixel files for both the pre-selected 2 minute cadence data and 30 minutes FFI targets. The tests are designed to help the user identify false positive scienarios with either instrumental or astrophysical. 

The program was developed with the growing exoplanet citizen science community in mind [e.g., @fischer12; @Christiansen2018; @eisner19], and thus requires no prior knowledge of Python coding. Additinally, the code can provide a useful tool for astronomical researchers and students.  ``LATTE`` is has a user-friendly interface that prompts the user to enter the TIC ID [@Stassun19] of the target star, the observatonal sector(s) to analyse and the location of the transit-like event, which can be identified using an interactive slider. All further options, such as whether to create a PDF file and whether to model the data, are indicated using 'tick-boxes' in the Graphical User Interface (GUI). Once this minimal information has been entered, the scipt downloads and handles the data without further user-interference in oder to provide essential exoplanet diagnostic plots and figures, which are summarised in a PDF document. 

The default setting for ``LATTE`` is to download and analyse the 2-minute cadence pre-selected data where the lightcurves were extracted and detrended by the SPOC pipeline [SPOC, @jenkins16]. However, run the program in the FFI mode (specified in the command line when the program is executed). In FFI mode the data is downloaded using TESScut [@brasseur2019] and the data detrended using the ``Astropy`` PCA module [@astropy], an inerative non-linear filter [@Aigrain04] to estimate and subtract residual systematics. Unlike the TESS 2-minute cadence targets, the FFIs do not come with a pre-identified optimal aperture. By default, the user is given the option to manually select the pixels within each aperture for a 'large' and a 'small' aperture. This GUI opens automatically and the two apertures are selected by clicking on the desired pixels in the interface. The two (large and small aperture) lightcurves are simultaneously displayed. Alternative, the aperture size can be defined automatically by the program using a threshold flux value and centred on the given target coordinates.



# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

I thank the all of the Planet Hunters TESS volunteers who prompted me to write a program that may be helpful in identifying the signal of true exoplanets in the TESS data and the support from the entire Zooniverse team and Oxford exoplanet group.

# References


