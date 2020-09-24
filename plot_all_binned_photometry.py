import numpy as np
import matplotlib.pyplot as plt
import betapic as bp
from astropy.table import Table


import os, sys
import datetime
runtime = os.path.abspath((sys.argv[0])) + " run at " + datetime.datetime.now().strftime("%c")
tyb = dict(color='black', fontsize=20)

from kepler3 import *
from exorings3 import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

t_mid = 58009. * u.day # time when star has closest approach to planet

f_brite = Table.read('binned_flux_brite.dat', format='ascii.ecsv')
f_astep = Table.read('binned_flux_astep.dat', format='ascii.ecsv')
f_bring = Table.read('binned_flux_bring.dat', format='ascii.ecsv')

good=(f_bring['ferr']>0.00001)

f_bring['ferr'][good] = 0.00001

figc, (a1, a2, a3) = plt.subplots(3,1,figsize=(10,6), constrained_layout=True, sharex=True, sharey=True)

a1.errorbar(f_brite['time'],  f_brite['flux'], yerr=f_brite['ferr'], fmt='.', alpha=0.2)
a2.errorbar(f_astep['time'],  f_astep['flux'], yerr=f_astep['ferr'], fmt='.', alpha=0.2)
a3.errorbar(f_bring['time'],  f_bring['flux'], yerr=f_bring['ferr'], fmt='.', alpha=0.2)

a1.set_xlim(57600, 58350)
a1.set_ylim(-0.06,0.06)
a1.text(0.05, 0.80, 'BRITE', ha='left', va='bottom', transform=a1.transAxes, **tyb)
a2.text(0.05, 0.80, 'ASTEP', ha='left', va='bottom', transform=a2.transAxes, **tyb)
a3.text(0.05, 0.80, 'BRING', ha='left', va='bottom', transform=a3.transAxes, **tyb)

a1.text(0.05, 0.80, 'BRITE', ha='left', va='bottom', transform=a1.transAxes, **tyb)
a1.text(0.98, 0.90, runtime, ha='right', va='bottom', transform=a1.transAxes, fontsize=8)
plt.draw()
plt.savefig('all_binned_photometry.pdf')

