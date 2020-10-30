from astropy.io import ascii
from astropy.table import Table
import numpy as np

## astep data
fname_astep17 = 'data/astep/bjd_bp_flux_2017.csv'
fname_astep18 = 'data/astep/bjd_bp_flux_2018.csv'

fname_out = 'astep_all.fits'

data_astep17 = ascii.read(fname_astep17)
data_astep18 = ascii.read(fname_astep18)

astep17flux = (data_astep17['FLUX1']-4*data_astep17['FSKY1'])/(data_astep17['FLUX3']-4*data_astep17['FSKY3'])
astep18flux = (data_astep18['FLUX1']-4*data_astep18['FSKY1'])/(data_astep18['FLUX3']-4*data_astep18['FSKY3'])
astep17nflux = astep17flux/np.nanmedian(astep17flux)-1
astep18nflux = astep18flux/np.nanmedian(astep18flux)-1

x_astep_all = np.append(data_astep17['BJD'], data_astep18['BJD'])
y_astep_all = np.append(astep17nflux, astep18nflux)

t = Table([x_astep_all, y_astep_all], names=('time','flux'))

t.write(fname_out)
print('read in {} and wrote out {}'.format(fname_astep17,fname_out))
