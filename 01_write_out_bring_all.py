from astropy.io import ascii
from astropy.table import Table

fname_bring = 'data/bring/Reduced_2110269.fits'
fname_bring = 'data/bring/Reduced_2110269_v20190625.fits'

fname_out = 'bring_all.fits'

data_bring = Table.read(fname_bring)

time = data_bring['jd']-2400000
y_bring = -data_bring['reduced']

t = Table([time, y_bring], names=('time','flux'))

t.write(fname_out)
print('read in {} and wrote out {}'.format(fname_bring,fname_out))
