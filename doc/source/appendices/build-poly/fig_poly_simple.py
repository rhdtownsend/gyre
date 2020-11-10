#!/usr/bin/env python
#
# Simple polytrope plot

# Imports

import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np

from pygyre import *
from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Read plot data

tab = read_model('poly.simple.h5')

# Do the plot

fig, ax = plt.subplots()

ax.plot(tab['x'], tab['theta'], color=ORANGE, label=r'$\theta$')
ax.plot(tab['x'], tab['rho/rho_c'], color=SKY_BLUE, label=r'$\rho/\rho_{\rm c}$')
ax.plot(tab['x'], tab['P/P_c'], color=BLUE_GREEN, label=r'$P/P_{\rm c}$')
ax.plot(tab['x'], tab['M_r/M'], color=BLACK, label=r'$M_{r}/M$')

ax.set_xlabel(r'$z/z_{\rm s}$')

ax.legend()

# Write out the figure

fig.tight_layout()
fig.savefig('fig_poly_simple.svg')
