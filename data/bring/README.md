# BRING photometry

Columns are:

JD == Julian Date

Raw == with initial reduction (spatio-temporal calibration, interpixel calibration), but no further detrending (moon/daily/lst/…). Units are V-magnitude

Reduced == with initial reduction (spatio-temporal calibration, interpixel calibration), detrended (moon/daily/lst/…). Units are delta-V-magnitude (median should be approximately zero).

ReducedHF == same as reduced, but with a 3-day moving median removed, so only High-Frequency components remain. This is for Delta Scuti pulsation analysis.
