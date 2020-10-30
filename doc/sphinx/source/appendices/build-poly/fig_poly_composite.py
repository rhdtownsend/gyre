#!/usr/bin/env python
#
# Composite polytrope plot

# Imports

import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np

from read_poly import *
from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Read plot data

p = read_poly('poly.composite.h5')

rho = p['t_z']*p['theta']**p['n_poly_z']

gamma_0 = (p['n_poly'][0]+1)/p['n_poly'][0]
gamma_z = (p['n_poly_z']+1)/p['n_poly_z']

P = 1/p['B_z']*(p['n_poly'][0]+1)/(p['n_poly_z']+1)*p['t_z']**(2-gamma_z)*rho**gamma_z

m = p['mu_z']/p['mu_z'][-1]

# Do the plot

fig, ax = plt.subplots()

ax.plot(p['z']/p['z_s'], p['theta'], color=ORANGE, label=r'$\theta$')
ax.plot(p['z']/p['z_s'], rho, color=SKY_BLUE, label=r'$\rho/\rho_{\rm c}$')
ax.plot(p['z']/p['z_s'], P, color=BLUE_GREEN, label=r'$P/P_{\rm c}$')
ax.plot(p['z']/p['z_s'], m, color=BLACK, label=r'$m/M$')

ax.set_xlabel(r'$z/z_{\rm s}$')

ax.legend()

# Write out the figure

fig.tight_layout()
fig.savefig('fig_poly_composite.svg')
