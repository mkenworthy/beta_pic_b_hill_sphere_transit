# Beta Pictoris b Hill sphere transit

Data and figures for the paper on the 2017 Beta Pictoris b Hill sphere transit.

Data reduction and analysis are carried out with Python scripts with numerical prefixes indicating the order that the scripts are executed.

The telescopes are: ASTEP 400, bRing, BRITE, HST.

## Data

Data from the telescopes are stored in:

    data/astep
    data/bring
    data/brite
    data/hst

### Processing photometric data to common flux scale

The raw data is contained in `data/` and processed with the `01_*` scripts into a standard epoch and magnitude format.

    python 01_write_out_astep_all.py
    python 01_write_out_bring_all.py
    python 01_write_out_brite_all.py

Output is:

       astep_all.fits
       bring_all.fits
       brite_all.fits


### Remove outliers and rebin data

     python 02_write_out_binned_flux.py

Output is:

       binned_flux_bring.dat
       binned_flux_astep.dat
       binned_flux_brite.dat

### Fit a CPD disk to artificial data

Generate example of fitting CPD orientation to artificial data set:

     python 03_fit_tilted_thin_disk_model_chisq.py

### Search for CPD disk in data

     python 04_cpdfit_combined.py 0.30
     python 04_cpdfit_combined.py 0.60

### Fit data with model of 1981 eclipse event

     python 05_fit_1981_all_data_with_offset.py

