# ASTEP 400 data

## 2020 reduction

The files contain the following columns

    BJD: Barycentric Julian Date
    FCAL0: Beta Pic flux divided by its daily median 
    FCAL1: Beta Pic flux with optimal calibration
    FCAL2: Beta Pic flux with a simple calibration using reference star #3 to calculate the daily median and reference star #2 for variations on shorter timescales
    FCAL3: Beta Pic flux with a simple calibration using reference star #5 to calculate the daily median and reference star #2 for variations on shorter timescales
    SUNELEV: Sun elevation (in degrees)
    AIRMASS: Beta Pic airmass
    SKY: Level of the sky (in ADU)
   FLAG: 0 if the stars were not correctly detected or the guiding was poor. 1 otherwise. 

I believe that for the paper, the “optimal calibration” FCAL1 should be used. It is calculated as follows: 
For BJD<57970, we use reference star #3 to calculate the daily median
Afterwards, we use  reference star #5 to calculate the daily median. Something changed, probably resulting from manipulations of the telescope and we noticed that reference star #3 was giving poorer stability.
Between BJD>57891 & BJD<57907, we use reference star #3 to calculate the daily median but apply a 0.985 coefficient: Ice was stuck on the window of the camera box, partially obscuring the light and affecting stars differentially. 

As you know we had issues with the presence of snow on the camera entrance window. BJD 57891-57907 occurred after a large snow storm and we later realized that the window had not been cleaned. Then, using 2 different reference stars (#3 and #5) was the simplest way of removing the trends. We can see that some snow deposits still affected the long-time trends. 

Djamel’s new analysis has decreased the influence of moonshine, but it is still present. For the paper, it is probably better to focus on the data when SUNELEV < -13°, (even though the oscillations are visible, the global level is very hard to estimate when the sky is too bright). The pulsations have *not* been removed from those lightcurves. 

    https://mail.google.com/mail/u/0/#search/tristan/FMfcgxwKjKnfjqrwmmBzBdZGnzjLTFDq

## 2018 reduction

`bjd_bp_flux_2017.csv` : 2017 ASTEP aperture flux of the beta Pictoris field (csv file)
`bjd_bp_flux_2018.csv` : 2018 ASTEP aperture flux of the beta Pictoris field (csv file)
`betapic_field.pdf`    : Field (1 deg x 1deg) of beta Pic, showing positions of the 18 stars for which aperture photometry was calculated


-----
Header of each csv file

BJD,EXPTIME,SUNELEV,AIRMASS,TRCKQLT,SKY,RPHOT,XS1,YS1,FLUX1,FSKY1,....   ....       ....,XS18,YS18,FLUX18,FSKY18

BJD: Barycentric Julian Day
EXPTIME: Exposure time (s)
SUNELEV: Sun elevation (degrees)
AIRMASS: Air mass
TRCKQLT: Tracking Quality : 0 = bad guiding, 1 = good guiding
SKY: Sky luminosity (ADU) 
RPHOT: Aperture radius (Pixels)
XS1,YS1: Position of the star 1 ( Beta Pic) on the CCD
FLUX1: Aperture flux (ADU) of the star 1 (Beta Pic)
FSKY1: Aperture Flux (ADU) of the Sky close to the star 1 (aperture of the sky = 0.5 x aperture of the star)
To calculate the flux of the star 1 : FLUX = FLUX1-4.*FSKY1
....

....
XS18,YS18: Position of the star 18 (see betapic_field.pdf) on the CCD
FLUX18: Aperture flux (ADU) of the star 18
FSKY18: Aperture Flux (ADU) of the Sky close to the star 18 (aperture of the sky = 0.5 x aperture of the star)
To calcule the flux of the star 18 : FLUX = FLUX18-4.*FSKY18
