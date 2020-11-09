#!/usr/bin/env python
#
# Frequency vs radial order

# Imports

import matplotlib.pyplot as plt
import pygyre as pg

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Read data

s = pg.read_output('summary.h5')

# Do the plot

fig, ax = plt.subplots()

ax.plot(s['n_pg'], s['freq'].real, color=SKY_BLUE)

ax.set_xlabel(r'$n_{\rm pg}$')
ax.set_ylabel(r'$\nu\,[{\rm cyc\,d^{-1}}]$')

ax.xaxis.get_major_locator().set_params(integer=True)

fig.tight_layout()
fig.savefig('fig_freq.svg')
